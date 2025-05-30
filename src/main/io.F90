module io
#ifdef MPI
  use mpi
  use mpi_utils, only: mpi_subarray_default, mpi_subarray_gravity, mpi_subarray_spc, mpi_type_sink_prop, mpi_subarray_extgrtv, mpi_type_sink_prop_array
#endif
  use derived_types, only: sink_prop
  implicit none

  logical, public :: legacy = .true.

  interface read_var
    module procedure read_int4
    module procedure read_real8
    module procedure read_array_1d_char
    module procedure read_array_3d_real8
    module procedure read_array_spc
    module procedure read_array_1d_sink
  end interface read_var

  interface write_var
    module procedure write_int4
    module procedure write_real8
    module procedure write_array_1d_char
    module procedure write_array_3d_real8
    module procedure write_array_spc
    module procedure write_array_1d_sink
  end interface write_var

  public :: read_var, read_dummy_recordmarker, read_extgrv_array
  public :: write_var, write_dummy_recordmarker, write_extgrv_array
  public :: write_string, write_string_master
  public :: open_file_read, close_file, open_file_write, open_file_write_ascii
  private

#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: offset = 0
  integer :: ierr
#endif

contains

#ifdef MPI
subroutine get_file_end(fh, end_bytes)
  integer, intent(in) :: fh
  integer(kind=MPI_OFFSET_KIND), intent(out) :: end_bytes
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer :: method = 2

  if (method == 1) then
    ! MPI_File_sync is very slow on some file systems, e.g. Lustre
    call mpi_file_sync(fh, ierr)
    call MPI_File_get_size(fh, end_bytes, ierr)
  elseif(method == 2) then
    end_bytes = offset

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
  character(len=*), intent(in) :: filename
  integer, intent(out) :: fh
#ifdef MPI
  offset = 0 ! reset the offset counter
  call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
#else
  open(newunit=fh, file=filename, status='replace', form='unformatted', action='write', access='stream')
#endif
end subroutine open_file_write

subroutine open_file_write_ascii(filename, fh)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: fh
#ifdef MPI
  call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
#else
  open(newunit=fh, file=filename, status='replace', form='formatted', action='write')
#endif
end subroutine open_file_write_ascii

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

subroutine read_array_3d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend, grav)
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend
  logical, optional, intent(in) :: grav
  logical :: gravity
#ifdef MPI
  integer :: nbuff, itype
#endif

  if (present(grav)) then
    gravity = grav
  else
    gravity = .false.
  end if

#ifdef MPI
  if (gravity) then
    itype = mpi_subarray_gravity
  else
    itype = mpi_subarray_default
  end if
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)

  call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, itype, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, itype)
#else
  read(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine read_array_3d_real8

subroutine read_array_spc(fh, arr, istart, iend, jstart, jend, kstart, kend, lstart, lend)
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend, lstart, lend
#ifdef MPI
  integer :: nbuff
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)*(lend-lstart+1)

  call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, mpi_subarray_spc, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend, jstart:jend, kstart:kend, lstart:lend), nbuff, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpi_subarray_spc)
#else
  read(fh) arr(istart:iend, jstart:jend, kstart:kend, lstart:lend)
#endif

end subroutine read_array_spc

subroutine read_array_1d_sink(fh, arr, istart, iend)
  integer, intent(in) :: fh
  type(sink_prop), allocatable, intent(inout) :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend
#ifdef MPI
  integer :: nbuff
  nbuff = (iend-istart+1)

  call mpi_file_set_view(fh, offset, mpi_type_sink_prop, mpi_type_sink_prop_array, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend), nbuff, mpi_type_sink_prop, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpi_type_sink_prop_array)

#else
  read(fh) arr(istart:iend)
#endif

end subroutine read_array_1d_sink

subroutine read_extgrv_array(fh, arr)
  use grid
  integer, intent(in) :: fh
  real(8), allocatable, intent(inout) :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer :: istart, iend, jstart, jend, kstart, kend
#ifdef MPI
  integer :: nbuff
#endif

  istart = gis; iend = gie
  jstart = gjs; jend = gje
  kstart = gks; kend = gke
  if (gis==gis_global) istart = gis-2
  if (gie==gie_global) iend = gie+2
  if (gjs==gjs_global) jstart = gjs-2
  if (gje==gje_global) jend = gje+2
  if (gks==gks_global) kstart = gks-2
  if (gke==gke_global) kend = gke+2

#ifdef MPI
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)

  call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, mpi_subarray_extgrtv, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpi_subarray_extgrtv)
#else
  read(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine read_extgrv_array

subroutine read_dummy_recordmarker(fh)
  integer, intent(in) :: fh
  integer :: dummy

  if (legacy) then
    call read_int4(fh, dummy)
  endif

end subroutine read_dummy_recordmarker

subroutine write_dummy_recordmarker(fh)
  integer, intent(in) :: fh
  integer, parameter :: dummy = -1 ! Assuming record marker is the size of an integer

  if (legacy) then
    call write_int4(fh, dummy)
  endif

end subroutine write_dummy_recordmarker

subroutine write_int4(fh, var)
  integer, intent(in) :: fh
  integer, intent(in) :: var
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, var, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, MPI_INTEGER4)
#else
  write(fh) var
#endif

end subroutine write_int4

subroutine write_real8(fh, var)
  integer, intent(in) :: fh
  real(8), intent(in) :: var
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, var, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, MPI_REAL8)
#else
  write(fh) var
#endif

end subroutine write_real8

subroutine write_array_3d_real8(fh, arr, istart, iend, jstart, jend, kstart, kend, grav)
  integer, intent(in) :: fh
  real(8), intent(in), allocatable :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend
  logical, optional, intent(in) :: grav
  logical :: gravity
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: nbuff, itype
#endif

  if (present(grav)) then
    gravity = grav
  else
    gravity = .false.
  end if

#ifdef MPI
  if (gravity) then
    itype = mpi_subarray_gravity
  else
    itype = mpi_subarray_default
  end if

  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, itype, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, itype)
#else
  write(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine write_array_3d_real8

subroutine write_array_spc(fh, arr, istart, iend, jstart, jend, kstart, kend, lstart, lend)
  integer, intent(in) :: fh
  real(8), intent(in), allocatable :: arr(:,:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend, lstart, lend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: nbuff

  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)*(lend-lstart+1)

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, mpi_subarray_spc, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend, jstart:jend, kstart:kend, lstart:lend), nbuff, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpi_subarray_spc)
#else
  write(fh) arr(istart:iend, jstart:jend, kstart:kend, lstart:lend)
#endif

end subroutine write_array_spc

subroutine write_array_1d_sink(fh, arr, istart, iend)
  integer, intent(in) :: fh
  type(sink_prop), intent(in), allocatable :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: nbuff

  nbuff = (iend-istart+1)

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_BYTE, mpi_type_sink_prop_array, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend), nbuff, mpi_type_sink_prop, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpi_type_sink_prop_array)
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
  integer :: nbuff, n

  nbuff = (iend - istart + 1) * 10

  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_CHARACTER, MPI_CHARACTER, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend), nbuff, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  do n = 1, nbuff
   call update_offset(fh, MPI_CHARACTER)
  end do

#else
  write(fh) arr(istart:iend)
#endif

end subroutine write_array_1d_char

subroutine write_extgrv_array(fh, arr)
  use grid
  integer, intent(in) :: fh
  real(8), intent(in), allocatable :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer :: istart, iend, jstart, jend, kstart, kend
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND) :: end_bytes
  integer :: nbuff
#endif

  istart = gis; iend = gie
  jstart = gjs; jend = gje
  kstart = gks; kend = gke
  if (gis==gis_global) istart = gis-2
  if (gie==gie_global) iend = gie+2
  if (gjs==gjs_global) jstart = gjs-2
  if (gje==gje_global) jend = gje+2
  if (gks==gks_global) kstart = gks-2
  if (gke==gke_global) kend = gke+2

#ifdef MPI
  nbuff = (iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)
  call get_file_end(fh, end_bytes)
  call mpi_file_set_view(fh, end_bytes, MPI_REAL8, mpi_subarray_extgrtv, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_write_all(fh, arr(istart:iend,jstart:jend,kstart:kend), nbuff, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
  call update_offset(fh, mpi_subarray_extgrtv)
#else
  write(fh) arr(istart:iend,jstart:jend,kstart:kend)
#endif

end subroutine write_extgrv_array

subroutine write_string(fh, str, advance)
  integer, intent(in) :: fh
  character(len=*), intent(in) :: str
  logical, optional, intent(in) :: advance
  logical :: adv

  if (present(advance)) then
    adv = advance
  else
    adv = .true.
  end if

#ifdef MPI
  call MPI_File_write_shared(fh, str, len_trim(str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  if (adv) call MPI_File_write_shared(fh, new_line(''), 1, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
#else
  if (adv) then
    write(fh, '(a)') trim(str)
  else
    write(fh, '(a)', advance='no') trim(str)
  end if
#endif
end subroutine write_string

subroutine write_string_master(fh, str, advance)
  use mpi_utils, only: barrier_mpi, myrank
  integer, intent(in) :: fh
  character(len=*), intent(in) :: str
  logical, optional, intent(in) :: advance
  logical :: adv

  if (present(advance)) then
    adv = advance
  else
    adv = .true.
  end if

  call barrier_mpi
  if (myrank==0) then
    call write_string(fh, str, advance=adv)
  end if
  call barrier_mpi

end subroutine write_string_master

end module io
