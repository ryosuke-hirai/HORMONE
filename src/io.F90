module io
#ifdef MPI
  use mpi
#endif
  implicit none

  interface read_var
    module procedure read_int4
    module procedure read_real8
    module procedure read_array_1d_char
    module procedure read_array_3d_real8
    module procedure read_array_4d_real8
    module procedure read_array_1d_sink
  end interface read_var

  public :: read_var, read_dummy_recordmarker, open_file_read, close_file
  private

#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: offset = 0
  integer :: ierr
#endif

contains

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

  ! TODO MPI
#ifdef MPI
  print*, "ERROR: read_array_1d_char not implemented for MPI"
  call mpi_abort(MPI_COMM_WORLD, 1, ierr)
#else
  read(fh) arr(istart:iend)
#endif

end subroutine read_array_1d_char

subroutine read_array_3d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend)
  use mpi_utils, only: type_mpi_array
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend
#ifdef MPI
  integer :: nbuff
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)

  call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, type_mpi_array, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, type_mpi_array)
#else
  read(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine read_array_3d_real8

subroutine read_array_4d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend, lstart, lend)
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend, lstart, lend

  ! TODO MPI
#ifdef MPI
  print*, "ERROR: read_array_4d_real8 not implemented for MPI"
  call mpi_abort(MPI_COMM_WORLD, 1, ierr)
#else
  read(fh) arr(istart:iend, jstart:jend, kstart:kend, lstart:lend)
#endif

end subroutine read_array_4d_real8

subroutine read_array_1d_sink(fh, arr, istart, iend)
  use sink_mod, only: sink_prop
  integer, intent(in) :: fh
  type(sink_prop), allocatable, intent(inout) :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend

  ! TODO MPI
#ifdef MPI
  print*, "ERROR: read_array_1d_sink not implemented for MPI"
  call mpi_abort(MPI_COMM_WORLD, 1, ierr)
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

end module io
