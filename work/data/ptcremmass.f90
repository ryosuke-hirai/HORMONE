program premmass

! purpose: To calculate remaining masses using tracer particle data.

  implicit none

  integer ej,ub,n, label, tn
  real*8 time, Macc, Mstr, Msta, Mini, x1,x3,mass
  real*8,parameter:: Msun = 1.989d33
  character*50 ptcfile

  open(unit=1,file='remmassptc.dat',status='replace')
  write(1,'(a10,3a16)')'time','Mstar','Mstrip','Maccrete'

  open(unit=10,file='ptc00000286s.dat',status='old')
  read(10,'()'); read(10,'()')

  Macc = 0d0 ; Mstr = 0d0 ; Msta = 0d0
  do n = 1, 1000000
   read(10,*,end=221)label, ej, ub, mass, x1, x3
   if(ej==1.and.ub==0)then
    Msta = Msta + mass
   elseif(ej==1.and.ub==1)then
    Mstr = Mstr + mass
   elseif(ej==2.and.ub==0)then
    Macc = Macc + mass
   end if
  end do
221 close(10)

  write(1,'(i10,3(1PE16.8e2))')tn,Msta/Msun, Mstr/Msun, Macc/Msun


 do tn = 300, 11900, 100
  write(ptcfile,'(a3,i8.8,a5)')'ptc',tn,'s.dat'
  open(unit=10,file=ptcfile,status='old')
  read(10,'()'); read(10,'()')

  Macc = 0d0 ; Mstr = 0d0 ; Msta = 0d0
  do n = 1, 1000000
   read(10,*,end=222)label, ej, ub, mass, x1, x3
   if(ej==1.and.ub==0)then
    Msta = Msta + mass
   elseif(ej==1.and.ub==1)then
    Mstr = Mstr + mass
   elseif(ej==2.and.ub==0)then
    Macc = Macc + mass
   end if
  end do
222 close(10)

  write(1,'(i10,3(1PE16.8e2))')tn,Msta/Msun, Mstr/Msun, Macc/Msun

 end do



end program premmass
