program extract

! purpose: To shrink the ejecta data.

  implicit none

  integer count,i,ejdatnum,every
  real*8 psmass,ejectadistance
  real*8,allocatable,dimension(:)::t_ej,d_ej,p_ej,e_ej,v_ej
  character*50 ejtbinfile

  every = 6

 ejtbinfile = 'ejtbin_8.7_1.0B_dis4.80E+13cm.dat'
 
 open(unit=90, file=ejtbinfile, status='old',form='unformatted')
 count = -2
 do
  read(90,end=311)
  count = count+1
 end do
311 close(90)!counting lines of ejtfile

 ejdatnum = count
 allocate( t_ej(count),d_ej(count),p_ej(count),e_ej(count),v_ej(count) )

 open(unit=90, file=ejtbinfile, status='old',form='unformatted')

 read(90)psmass
 read(90)ejectadistance

 do i = 1, count
  read(90)t_ej(i), d_ej(i), p_ej(i), e_ej(i), v_ej(i)
 end do
 close(90) !reading ejecta data


 open(unit=10,file='hogebin.dat',status='replace',form='unformatted')
 write(10)psmass
 write(10)ejectadistance

 do i = 1, count, every
  write(10)t_ej(i), d_ej(i), p_ej(i), e_ej(i), v_ej(i)
 end do
 close(10)

end program
