module io
#ifdef MPI
  use mpi
#endif
  use sink_mod, only: sink_prop
  implicit none

  interface read_var
    module procedure read_int4
    module procedure read_real8
    module procedure read_array_1d_char
    module procedure read_array_3d_real8
    module procedure read_array_4d_real8
    module procedure read_array_1d_sink
  end interface read_var

  interface write_var
    module procedure write_int4
    module procedure write_real8
    module procedure write_array_1d_char
    module procedure write_array_3d_real8
    module procedure write_array_4d_real8
    module procedure write_array_1d_sink
  end interface write_var

  public :: read_var, read_dummy_recordmarker
  public :: write_var, write_dummy_recordmarker
  public :: open_file_read, close_file, open_file_write
  private

#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: offset = 0
  integer :: ierr
#endif

contains

#ifdef MPI
subroutine get_file_end(fh, end_bytes)
  integer, intent(in) :: fh
  integer :: ierr
  integer(kind=MPI_OFFSET_KIND), intent(out) :: end_bytes
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer :: method = 1

  if (method == 1) then
    ! This may or may not be better than the method below... idk yet
    call mpi_file_sync(fh, ierr)
    call MPI_File_get_size(fh, end_bytes, ierr)
  else
! --- Alternative method------------------------------------
    ! Reset view to default ("absolute" view in bytes)
    call mpi_file_set_view(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
    ! Move individual file positions to absolute end of file
    call mpi_file_seek(fh, 0_MPI_OFFSET_KIND, MPI_SEEK_END, ierr)
    call mpi_file_get_position(fh, disp, ierr)
    call mpi_file_get_byte_offset(fh, disp, end_bytes, ierr)
! ----------------------------------------------------------
  endif

end subroutine get_file_end
#endif

subroutine open_file_write(filename, fh)
  use mpi_utils, only: myrank
  character(len=*), intent(in) :: filename
  integer, intent(out) :: fh
#ifdef MPI
  integer :: ierr
#endif

  if (myrank == 0) then
    ! Create new file
    open(newunit=fh, file=filename, status='replace', form='unformatted', action='write', access='stream')
#ifdef MPI
    ! Close the file so that other ranks can open it in MPI mode
    close(fh)
#endif
  end if

#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR + MPI_MODE_APPEND, MPI_INFO_NULL, fh, ierr)
  call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

end subroutine open_file_write

subroutine open_file_read(filename, fh)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: fh
#ifndef MPI
  integer :: istat
#endif

#ifdef MPI
  offset = 0   ! reset the offset counter
  call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  if (ierr /= 0) then
    print*,'Error opening binary dump file'
    print'(3a)','File name = "',trim(filename),'"'
    call mpi_abort(MPI_COMM_WORLD, 1, ierr)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)
#else
  open(newunit=fh,file=filename,status='old',form='unformatted',iostat=istat, access='stream')
  if(istat/=0)then
    print*,'Error opening binary dump file'
    print'(3a)','File name = "',trim(filename),'"'
    stop
  end if
#endif

end subroutine open_file_read

subroutine close_file(fh)
  integer, intent(inout) :: fh

#ifdef MPI
  call mpi_file_close(fh, ierr)
  call mpi_barrier(MPI_COMM_WORLD, ierr)
#else
  close(fh)
#endif

end subroutine close_file

#ifdef MPI
subroutine update_offset(fh, dtype)
  integer, intent(in) :: fh, dtype
  integer(kind=MPI_OFFSET_KIND) :: bytes

  ! get the number of bytes read
  call mpi_file_get_type_extent(fh, dtype, bytes, ierr)
  offset = offset + bytes
  call mpi_barrier(MPI_COMM_WORLD, ierr)

end subroutine update_offset
#endif

subroutine read_int4(fh, var)
  integer, intent(in) :: fh
  integer, intent(out) :: var

#ifdef MPI
  call mpi_file_set_view(fh, offset, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, var, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, MPI_INTEGER4)
#else
  read(fh) var
#endif

end subroutine read_int4

subroutine read_real8(fh, var)
  integer, intent(in) :: fh
  real(8), intent(out) :: var

#ifdef MPI
  call mpi_file_set_view(fh, offset, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, var, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, MPI_REAL8)
#else
  read(fh) var
#endif

end subroutine read_real8

subroutine read_array_1d_char(fh, arr, istart, iend)
  integer, intent(in) :: fh
  character(len=10), allocatable, intent(inout) :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend
#ifdef MPI
  integer :: i, nbuff

  nbuff = (iend - istart + 1) * 10

  call mpi_file_set_view(fh, offset, MPI_CHARACTER, MPI_CHARACTER, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend), nbuff, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  do i = 1,nbuff
    call update_offset(fh, MPI_CHARACTER)
  end do
#else
  read(fh) arr(istart:iend)
#endif

end subroutine read_array_1d_char

subroutine read_array_3d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend)
  use mpi_utils, only: mpitype_array3d_real8
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend
#ifdef MPI
  integer :: nbuff
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)

  call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, mpitype_array3d_real8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpitype_array3d_real8)
#else
  read(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine read_array_3d_real8

subroutine read_array_4d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend, lstart, lend)
  use mpi_utils, only: mpitype_array4d_real8
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend, lstart, lend
#ifdef MPI
  integer :: nbuff
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)*(lend-lstart+1)

  call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, mpitype_array4d_real8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend, jstart:jend, kstart:kend, lstart:lend), nbuff, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpitype_array4d_real8)
#else
  read(fh) arr(istart:iend, jstart:jend, kstart:kend, lstart:lend)
#endif

end subroutine read_array_4d_real8

subroutine read_array_1d_sink(fh, arr, istart, iend)
  use mpi_utils, only: mpitype_sink_prop
  integer, intent(in) :: fh
  type(sink_prop), allocatable, intent(inout) :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend
#ifdef MPI
  integer :: nbuff
  nbuff = (iend-istart+1)

  call mpi_file_set_view(fh, offset, mpitype_sink_prop, mpitype_sink_prop, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend), nbuff, mpitype_sink_prop, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpitype_sink_prop)
#else
  read(fh) arr(istart:iend)
#endif

end subroutine read_array_1d_sink

subroutine read_dummy_recordmarker(fh, legacy)
  integer, intent(in) :: fh
  logical, intent(in) :: legacy
  integer :: dummy

  if (legacy) then
    call read_int4(fh, dummy)
  endif

end subroutine read_dummy_recordmarker

subroutine write_dummy_recordmarker(fh, legacy)
  integer, intent(in) :: fh
  logical, intent(in) :: legacy
  integer, parameter :: dummy = -1 ! Assmuing record marker is the size of an integer

  if (legacy) then
    call write_int4(fh, dummy)
  endif

end subroutine write_dummy_recordmarker

subroutine write_int4(fh, var)
  integer, intent(in) :: fh
  integer, intent(in) :: var
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: ierr

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, var, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
#else
  write(fh) var
#endif

end subroutine write_int4

subroutine write_real8(fh, var)
  integer, intent(in) :: fh
  real(8), intent(in) :: var
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: ierr

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, var, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
  write(fh) var
#endif

end subroutine write_real8

subroutine write_array_3d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend)
  use mpi_utils, only: mpitype_array3d_real8
  integer, intent(in) :: fh
  real(8), intent(in), allocatable :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: ierr, nbuff

  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, mpitype_array3d_real8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
  write(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine write_array_3d_real8

subroutine write_array_4d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend, lstart, lend)
  use mpi_utils, only: mpitype_array4d_real8
  integer, intent(in) :: fh
  real(8), intent(in), allocatable :: arr(:,:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend, lstart, lend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: ierr, nbuff

  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)*(lend-lstart+1)

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, mpitype_array4d_real8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend, jstart:jend, kstart:kend, lstart:lend), nbuff, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
  write(fh) arr(istart:iend, jstart:jend, kstart:kend, lstart:lend)
#endif

end subroutine write_array_4d_real8

subroutine write_array_1d_sink(fh, arr, istart, iend)
  use mpi_utils, only: mpitype_sink_prop
  integer, intent(in) :: fh
  type(sink_prop), intent(in), allocatable :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: ierr, nbuff

  nbuff = (iend-istart+1)

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, mpitype_sink_prop, mpitype_sink_prop, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend), nbuff, mpitype_sink_prop, MPI_STATUS_IGNORE, ierr)
#else
  write(fh) arr(istart:iend)
#endif

end subroutine write_array_1d_sink

subroutine write_array_1d_char(fh, arr, istart, iend)
  integer, intent(in) :: fh
  character(len=10), intent(in), allocatable :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: nbuff

  nbuff = (iend - istart + 1) * 10

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_CHARACTER, MPI_CHARACTER, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend), nbuff, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
#else
  write(fh) arr(istart:iend)
#endif

end subroutine write_array_1d_char

end module io
