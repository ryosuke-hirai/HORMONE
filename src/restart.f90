!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE RESTART
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To restart calculation from arbitrary timestep.
!          from INITIALCONDITION

subroutine restart

 use settings!,only:start,outstyle,gravswitch,dt_out,include_particles
 use grid
 use physval
 use gravmod
 use pressure_mod
 use ejtfilemod,only:inimass
 use particle_mod
 use merger_mod
 use constants

 implicit none

 character*30 startfile, bptfile
real*8 dfac, Mdot, vinf, mmm
!-----------------------------------------------------------------------------

 if(outstyle==1)then
  write(startfile,'(a8,i11.11,a5)')'data/bin',start,'s.dat'
 elseif(outstyle==2)then
  write(startfile,'(a8,i8.8,a4)')'data/bin',start,'.dat'
 end if

 open(unit=11,file=startfile,status='old',form='unformatted')

  read(11)tn,time,iniEtot,inimass,de_dt,domega_dt
  read(11) d (is:ie,js:je,ks:ke), &
           v1(is:ie,js:je,ks:ke), &
           v2(is:ie,js:je,ks:ke), &
           v3(is:ie,js:je,ks:ke), &
!!$           b1(is:ie,js:je,ks:ke), &
!!$           b2(is:ie,js:je,ks:ke), &
!!$           b3(is:ie,js:je,ks:ke), &
           e (is:ie,js:je,ks:ke), &
!           phi(is:ie,js:je,ks:ke), &
           grvphi(gis:gie,gjs:gje,gks:gke)
  if(gravswitch==3)then
   read(11)grvphiold(gis:gie,gjs:gje,gks:gke), &
           dt_old
  end if
  if(compswitch>=2)then
   read(11)spc(1:spn,is:ie,js:je,ks:ke)
  end if

 close(11)

!temp
 if(start==int(inifile2))then 
  ie=1200;gie=ie
  in = ie + 2 ; jn = je + 2 ; kn = ke + 2
  gin = gie + 2 ; gjn = gje + 2 ; gkn = gke + 2
! lmax = ie*je*ke
  lmax = (gie-gis+1)*(gje-gjs+1)*(gke-gks+1)

  dfac = dxi1(is+1)/dxi1(is)
  do i = 901, ie+2
   dxi1(i) = dxi1(i-1)*dfac
   idxi1(i) = 1d0/dxi1(i)
   xi1(i) = xi1(i-1) + dxi1(i)
   x1(i) = 0.75d0*(xi1(i)+xi1(i-1))*(xi1(i)*xi1(i)+xi1(i-1)*xi1(i-1)) &
                 /(xi1(i)*xi1(i)+xi1(i)*xi1(i-1)+xi1(i-1)*xi1(i-1))
   dx1(i) = x1(i) - x1(i-1)
   idx1(i) = 1d0/dx1(i)
  end do
  xi1e = xi1(ie)

  Mdot = 1d-5*msun/year
  vinf = 4d7
  mmm  = sum(d(is:900,js:je,ks:ke)*dvol(is:900,js:je,ks:ke))*2d0

  k = ks
  do j = js, je
   do i = 901, ie
    d (i,j,k) = Mdot/(4d0*pi*x1(i)*x1(i)*vinf)
    v1(i,j,k) = vinf
    v2(i,j,k) = 0d0
    v3(i,j,k) = 0d0
    grvphi(i,j,k) = -G*mmm/x1(i)
    grvphiold(i,j,k) = -G*mmm/x1(i)
    p (i,j,k) = -d(i,j,k)*grvphi(i,j,k)
    spc(1:spn,i,j,k) = spc(1:spn,900,j,k)
    imu(i,j,k) = 0.25d0*(6d0*spc(1,i,j,k)+spc(2,i,j,k)+2d0)
    eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    e (i,j,k) = eint(i,j,k) + 0.5d0*d(i,j,k)*vinf*vinf
   end do
  end do
  deallocate(spin_coeffr,spin_coefft)
  call metric
  deallocate( hg11,hg12,hg21,hg22,hg123 )
  call gravsetup

  open(unit=40,file='data/gridfile2.bin',status='replace',form='unformatted')

  write(40) x1  (gis-2:gie+2), xi1(gis-2:gie+2), &
       dxi1(gis-2:gie+2), dx1(gis-2:gie+2), &
       x2  (gjs-2:gje+2), xi2(gjs-2:gje+2), &
       dxi2(gjs-2:gje+2), dx2(gjs-2:gje+2), &
       x3  (gks-2:gke+2), xi3(gks-2:gke+2), &
       dxi3(gks-2:gke+2), dx3(gks-2:gke+2)
  close(40)

  open(unit=50,file='data/gridfile2.dat',status='replace')
  write(50,'()')
  write(50,'(2(a5),3(a15))')'#  i','j','x1(3)','x2(4)','dvol(5)'

  if(crdnt==2)then
   k=ks;j=js-1
   do i = is, ie,2
    write(50,'(2(i5),3(1x,1PE14.6e2))')i,j,x1(i),0d0,dvol(i,j,k)
   end do
   write(50,*)
  end if
  do k = ks,ke
   do j = js,je
    do i = is,ie,2
     write(50,'(2(i5),3(1x,1PE14.6e2))')i,j,x1(i),x2(j),dvol(i,j,k)
    end do
    if(ie/=1)write(50,*)
   end do
  !   if(ie==1)write(30,*)
  end do
  close(50)
 end if


!bptfile----------------------------------------------------------------
 if(include_particles.and.time>inifile)then
  if(outstyle==1)then
   write(bptfile,'(a8,i11.11,a5)')'data/bpt',start,'s.dat'
  elseif(outstyle==2)then
   write(bptfile,'(a8,i11.11,a4)')'data/bpt',start,'.dat'
  end if

  open(unit=81,file = bptfile, status='old',form='unformatted')
  allocate( ptcx(0:dim,1:maxptc), ptci(0:2,1:maxptc) )

  read(81)np,npl
  read(81)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)

  close(81)
 end if

 t_out = time + dt_out
 if(gravswitch==3)grvtime = time
! call readejecta

end subroutine restart
