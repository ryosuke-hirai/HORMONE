module io
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

contains

subroutine open_file_read(filename, un)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: un
  integer :: istat

  open(newunit=un,file=filename,status='old',form='unformatted',iostat=istat, access='stream')
  if(istat/=0)then
    print*,'Binary dump file not found'
    print'(3a)','File name = "',trim(filename),'"'
    stop
  end if

end subroutine open_file_read

subroutine close_file(un)
  integer, intent(in) :: un

  close(un)

end subroutine close_file

subroutine read_int4(un, var)
  integer, intent(in) :: un
  integer, intent(out) :: var

  read(un) var

end subroutine read_int4

subroutine read_real8(un, var)
  integer, intent(in) :: un
  real(8), intent(out) :: var

  read(un) var

end subroutine read_real8

subroutine read_array_1d_char(un, arr, istart, iend)
  integer, intent(in) :: un
  character(len=10), allocatable, intent(inout) :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend

  read(un) arr(istart:iend)

end subroutine read_array_1d_char

subroutine read_array_3d_real8(un, arr, istart, iend, jstart, jend, kstart, kend)
  integer, intent(in) :: un
  real(8), allocatable, intent(inout) :: arr(:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend

  read(un) arr(istart:iend,jstart:jend,kstart:kend)

end subroutine read_array_3d_real8

subroutine read_array_4d_real8(un, arr, istart, iend, jstart, jend, kstart, kend, lstart, lend)
  integer, intent(in) :: un
  real(8), allocatable, intent(inout) :: arr(:,:,:,:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend, jstart, jend, kstart, kend, lstart, lend

  read(un) arr(istart:iend, jstart:jend, kstart:kend, lstart:lend)

end subroutine read_array_4d_real8

subroutine read_array_1d_sink(un, arr, istart, iend)
  use sink_mod, only: sink_prop
  integer, intent(in) :: un
  type(sink_prop), allocatable, intent(inout) :: arr(:) ! use allocatable attribute to preserve lower and upper bound indices
  integer, intent(in) :: istart, iend

  read(un) arr(istart:iend)

end subroutine read_array_1d_sink

subroutine read_dummy_recordmarker(un, legacy)
  integer, intent(in) :: un
  logical, intent(in) :: legacy
  integer :: dummy

  if (legacy) read(un) dummy

end subroutine read_dummy_recordmarker

end module io
