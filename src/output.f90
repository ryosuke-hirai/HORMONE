module ejtfilemod

  implicit none

  real*8 t_ejt_out, dt_ejt_out, remmass, inimass, belmass, totmass, totangmom, totangmoma, totangmomb, totenergy, totenergya, totenergyb, ToW, ToWb, ToWt, totintene, totinteneb, totintenet, totgrvene,totgrveneb,totgrvenet,totkinene,totkineneb,totkinenet

end module ejtfilemod
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output data

subroutine output

  use settings,only:outstyle,gravswitch,compswitch,spn
  use grid
  use physval
  use constants
  use gravmod
  use ejtfilemod
  use particle_mod
  use merger_mod

  implicit none

  logical extflag
  character*35 pltfile, binfile, ptcfile, bptfile
  real*8 phih, vsq, pr, pt
  real*8,allocatable:: lapphi(:,:,:)

!----------------------------------------------------------------------------


!gridfile---------------------------------------------------------------
  if(tn==0.or.time==inifile)then
   open(unit=40,file='data/gridfile.bin',status='replace',form='unformatted')

   write(40) x1  (gis-2:gie+2), xi1(gis-2:gie+2), &
             dxi1(gis-2:gie+2), dx1(gis-2:gie+2), &
             x2  (gjs-2:gje+2), xi2(gjs-2:gje+2), &
             dxi2(gjs-2:gje+2), dx2(gjs-2:gje+2), &
             x3  (gks-2:gke+2), xi3(gks-2:gke+2), &
             dxi3(gks-2:gke+2), dx3(gks-2:gke+2)
   close(40)

   open(unit=50,file='data/gridfile.dat',status='replace')
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
  end if

  close(50)


!pltfile---------------------------------------------------------------------
  if(outstyle==1)then
   write(pltfile,'(a8,i11.11,a5)')'data/plt',int(time),'s.dat'
  elseif(outstyle==2)then
   write(pltfile,'(a8,i8.8,a4)')'data/plt',tn,'.dat'
  end if

  open(unit=30,file = pltfile, status='replace')

  write(30,'(a,i7,2(a,1PE12.4e2))') '#tn =',tn,'  time= ',time
  write(30,'(18(a15))')'d(6)','e(7)','p(8)','v1(9)','v2(10)','v3(11)','phi(12)','T(13)','mu(14)','H1(15)','He4(16)','He3(17)','C12(18)','N14(19)','O16(20)','Ne20(21)','Mg24(22)'

  if(gravswitch==1)then
   j = js ; k = ks
   do i = is, ie-1
    grvphi(i,j,k) = 0d0
    do n = i+1, ie
     grvphi(i,j,k) = grvphi(i,j,k) - d(n,j,k)*x1(n)*dxi1(n)
    end do
    grvphi(i,j,k) = G*(-mc(i)/x1(i)+4d0*pi*grvphi(i,j,k))
   end do
   grvphi(ie,j,k) = -G*mc(ie)/x1(ie)
  end if

!!$  allocate(lapphi(is:ie,js:je,ks:ke))
!!$  do k = ks, ke
!!$   do j = js, je
!!$    do i = is, ie
!!$     h = min(dx1(i),dx1(i+1),dx3(k),dx3(k+1))
!!$     lapphi(i,j,k) =((sum(lag11(-1:1,i,j,k)*grvphi(i-1:i+1,j,k)) + &
!!$                     sum(lag12(-1:1,i,j,k)*grvphi(i-1:i+1,j,k)) + &
!!$                     sum(lag31(-1:1,i,j,k)*grvphi(i,j,k-1:k+1)) + &
!!$                     sum(lag32(-1:1,i,j,k)*grvphi(i,j,k-1:k+1)) + &
!!$                     sum(lag21(-1:1,i,j,k)*grvphi(i-1:i+1,j,k))*2d0 - &
!!$                     6d0*grvphi(i,j,k))/(h*h) - &
!!$                     4d0*pi*G*d(i,j,k))/ &
!!$                     ((sum(lag11(-1:1,i,j,k)*grvphi(i-1:i+1,j,k)) + &
!!$                     sum(lag12(-1:1,i,j,k)*grvphi(i-1:i+1,j,k)) + &
!!$                     sum(lag31(-1:1,i,j,k)*grvphi(i,j,k-1:k+1)) + &
!!$                     sum(lag32(-1:1,i,j,k)*grvphi(i,j,k-1:k+1)) + &
!!$                     sum(lag21(-1:1,i,j,k)*grvphi(i-1:i+1,j,k))*2d0 - &
!!$                     6d0*grvphi(i,j,k))/(h*h) + &
!!$                     4d0*pi*G*d(i,j,k))
!!$    end do
!!$   end do
!!$  end do
  if(crdnt==2)then
   j =js;k=ks
   do i = is, ie,2
    write(30,'(18(1x,1PE14.6e2))')&
         d(i,j,k),e(i,j,k),p(i,j,k),v1(i,j,k),v2(i,j,k),v3(i,j,k), &
         grvphi(i,j,k),T(i,j,k),1d0/imu(i,j,k),spc(1:8,i,j,k)
   end do
   write(30,'()')
  end if

  do k = ks,ke
   do j = js,je
    do i = is,ie,2
     write(30,'(18(1x,1PE14.6e2))')&
          d(i,j,k),e(i,j,k),p(i,j,k),v1(i,j,k),v2(i,j,k),v3(i,j,k), &
          grvphi(i,j,k),T(i,j,k),1d0/imu(i,j,k),spc(1:8,i,j,k)
    end do
    if(ie/=1)write(30,*)
   end do
!   if(ie==1)write(30,*)
  end do

  close(30)


!binfile----------------------------------------------------------------
if(time>inifile)then
 if(outstyle==1)then
   write(binfile,'(a8,i11.11,a5)')'data/bin',int(time),'s.dat'
  elseif(outstyle==2)then
   write(binfile,'(a8,i8.8,a4)')'data/bin',tn,'.dat'   
  end if

  open(unit=10,file=binfile,status='replace',form='unformatted')

  write(10)tn,time,iniEtot,inimass,de_dt,domega_dt
  write(10) d (is:ie,js:je,ks:ke), &
            v1(is:ie,js:je,ks:ke), &
            v2(is:ie,js:je,ks:ke), &
            v3(is:ie,js:je,ks:ke), &
!            b1(is:ie,js:je,ks:ke), &
!            b2(is:ie,js:je,ks:ke), &
!            b3(is:ie,js:je,ks:ke), &
            e (is:ie,js:je,ks:ke), &
!            phi(is:ie,js:je,ks:ke), &
            grvphi(gis:gie,gjs:gje,gks:gke)
  if(gravswitch==3)then
   write(10)grvphiold(gis:gie,gjs:gje,gks:gke), &
            dt_old
  end if
  if(compswitch>=2)then
   write(10)spc(1:spn,is:ie,js:je,ks:ke)
  end if

  close(10)
 end if
  if(include_particles)then
!ptcfile----------------------------------------------------------------
  if(outstyle==1)then
   write(ptcfile,'(a8,i11.11,a5)')'data/ptc',int(time),'s.dat'
  elseif(outstyle==2)then
   write(ptcfile,'(a8,i8.8,a4)')'data/ptc',tn,'.dat'
  end if

  if(crdnt==1.and.je==1)then
   do n = 1, np
    extflag = .false.
    do k = ks, ke
     if(ptcx(2,n)>=xi3(k-1))then;if(ptcx(2,n)<xi3(k))then
      do i = is, ie
       if(ptcx(1,n)>=xi1(i-1))then;if(ptcx(1,n)<xi1(i))then
        if(e(i,js,k)+grvphi(i,js,k)*d(i,js,k)>0d0)then
         ptci(2,n) = 1
        else
         ptci(2,n) = 0
        end if
        extflag=.true.
        exit
       end if;end if
      end do
      if(extflag)exit
     end if;end if
    end do
   end do
  elseif(crdnt==2.and.ke==1)then
   do n = 1, np
    extflag = .false.
    pr = sqrt(ptcx(1,n)**2d0+ptcx(2,n)**2d0)
    pt = acos(ptcx(2,n)/pr)
    do j = js, je
     if(pt>=xi2(j-1))then;if(pt<xi2(j))then
      do i = is, ie
       if(pr>=xi1(i-1))then;if(pr<xi1(i))then
        if(e(i,j,ks)+grvphi(i,j,ks)*d(i,j,ks)>0d0)then
         ptci(2,n) = 1
        else
         ptci(2,n) = 0
        end if
        extflag=.true.
        exit
       end if;end if
      end do
      if(extflag)exit
     end if;end if
    end do
   end do
  end if

  open(unit=70,file = ptcfile, status='replace')

  write(70,'(a,i7,a,1PE12.4e2,2(a5,i9))')&
       '#tn =',tn,'  time= ',time,'np= ',np,'npl=',npl
  write(70,'(a7,2a3,3a15)')'label','ej','ub','mass','x1','x3'

  do i = 1, np
   write(70,'(i7,2i3,3(1PE15.7e2))')ptci(0:2,i),ptcx(0:2,i)
  end do

  close(70)

!bptfile----------------------------------------------------------------
  if(outstyle==1)then
   write(bptfile,'(a8,i11.11,a5)')'data/bpt',int(time),'s.dat'
  elseif(outstyle==2)then
   write(bptfile,'(a8,i8.8,a4)')'data/bpt',tn,'.dat'
  end if

  open(unit=80,file = bptfile, status='replace',form='unformatted')

  write(80)np,npl
  write(80)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)

  close(80)

  end if

!angmomfile---------------------------------------------------------------
  if(time==inifile)then
   inimass = 0d0
   do k = ks, ke ; do j = js, je ; do i = is, ie
    if(e(i,j,k)+d(i,j,k)*(grvphi(i,j,k))<=0d0)then
     inimass = inimass + dvol(i,j,k)*d(i,j,k)
    end if
   end do ; end do ; end do
   write(60,'(a,1PE16.8e2)')'initial_mass = ',inimass*2d0/msun
   write(60,'(24a16)')'time(1)','totmass(2)','remmass(3)','bermass(4)','totangmom(5)','totangmoma(6)','totangmomb(7)','totenergy(8)','totenergya(9)','totenergyb(10)','tot_T/W(11)','T/W(12)','T/Wb(13)','Eint(14)','Einta(15)','Eintb(16)','Egrav(17)','Egrava(18)','Egravb(19)','Ekin(20)','Ekina(21)','Ekinb(22)','domega_dt(23)'
  end if

  remmass = 0d0 ; belmass = 0d0 ; totmass = 0d0
  totangmom = 0d0 ; totangmoma = 0d0 ; totangmomb = 0d0
  totenergy = 0d0 ; totenergya = 0d0 ; totenergyb = 0d0
  ToW = 0d0 ; ToWb = 0d0 ; ToWt = 0d0
  totintene=0d0; totinteneb=0d0; totintenet=0d0
  totgrvene=0d0; totgrveneb=0d0; totgrvenet=0d0
  totkinene=0d0; totkineneb=0d0; totkinenet=0d0
! for cylindrical coordinates
  if(crdnt==1)then
   do k = ks, ke
    do j = js, je
     do i = is, ie
      totmass = totmass + d(i,j,k) * dvol(i,j,k)
      totangmom = totangmom + d(i,j,k)*v2(i,j,k)*dvol(i,j,k)*x1(i)
      totenergy = totenergy + (e(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)
      vsq = v1(i,j,k)*v1(i,j,k)+v2(i,j,k)*v2(i,j,k)+v3(i,j,k)*v3(i,j,k)
      totintene = totintene + (e(i,j,k)-0.5d0*d(i,j,k)*vsq)*dvol(i,j,k)
      totgrvene = totgrvene + 0.5d0*d(i,j,k)*grvphi(i,j,k)*dvol(i,j,k)
      totkinene = totkinene + 0.5d0*d(i,j,k)*vsq*dvol(i,j,k)
      ToWt = ToWt + 0.5d0*d(i,j,k)*dvol(i,j,k)*v2(i,j,k)*v2(i,j,k)
      if(e(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
       remmass = remmass + d(i,j,k)*dvol(i,j,k)
       totangmoma = totangmoma + d(i,j,k)*v2(i,j,k)*dvol(i,j,k)*x1(i)
       totenergya = totenergya + (e(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)
       ToW = ToW + 0.5d0*d(i,j,k)*dvol(i,j,k)*v2(i,j,k)*v2(i,j,k)
      end if
      if(e(i,j,k)+p(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
       belmass = belmass + d(i,j,k)*dvol(i,j,k)
       totangmomb= totangmomb+ d(i,j,k)*v2(i,j,k)*dvol(i,j,k)*x1(i)
       totenergyb= totenergyb+ (e(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)
       ToWb= ToWb+ 0.5d0*d(i,j,k)*dvol(i,j,k)*v2(i,j,k)*v2(i,j,k)
      end if
     end do
    end do
   end do
! for spherical coordinates
  elseif(crdnt==2)then
   do k = ks, ke
    do j = js, je
     do i = is, ie
      totmass = totmass + d(i,j,k) * dvol(i,j,k)
      totangmom = totangmom + d(i,j,k)*v3(i,j,k)*dvol(i,j,k)*spin_coeffr(i)*spin_coefft(j)
      totenergy = totenergy + (e(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)
      vsq = v1(i,j,k)*v1(i,j,k)+v2(i,j,k)*v2(i,j,k)+v3(i,j,k)*v3(i,j,k)
      totintenet = totintenet + (e(i,j,k)-0.5d0*d(i,j,k)*vsq)*dvol(i,j,k)
      totgrvenet = totgrvenet + 0.5d0*d(i,j,k)*grvphi(i,j,k)*dvol(i,j,k)
      totkinenet = totkinenet + 0.5d0*d(i,j,k)*vsq*dvol(i,j,k)
      ToWt = ToWt + 0.5d0*d(i,j,k)*dvol(i,j,k)*v3(i,j,k)*v3(i,j,k)
      if(e(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
       remmass = remmass + d(i,j,k)*dvol(i,j,k)
       totangmoma = totangmoma + d(i,j,k)*v3(i,j,k)*dvol(i,j,k)*spin_coeffr(i)*spin_coefft(j)
       totenergya = totenergya + (e(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)
       totintene = totintene + (e(i,j,k)-0.5d0*d(i,j,k)*vsq)*dvol(i,j,k)
       totgrvene = totgrvene + 0.5d0*d(i,j,k)*grvphi(i,j,k)*dvol(i,j,k)
       totkinene = totkinene + 0.5d0*d(i,j,k)*vsq*dvol(i,j,k)
       ToW = ToW + 0.5d0*d(i,j,k)*dvol(i,j,k)*v3(i,j,k)*v3(i,j,k)
      end if
      if(e(i,j,k)+p(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
       belmass = belmass + d(i,j,k)*dvol(i,j,k)
       totangmomb= totangmomb+ d(i,j,k)*v3(i,j,k)*dvol(i,j,k)*spin_coeffr(i)*spin_coefft(j)
       totenergyb= totenergyb+ (e(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)
       totinteneb = totinteneb + (e(i,j,k)-0.5d0*d(i,j,k)*vsq)*dvol(i,j,k)
       totgrveneb = totgrveneb + 0.5d0*d(i,j,k)*grvphi(i,j,k)*dvol(i,j,k)
       totkineneb = totkineneb + 0.5d0*d(i,j,k)*vsq*dvol(i,j,k)
       ToWb= ToWb+ 0.5d0*d(i,j,k)*dvol(i,j,k)*v3(i,j,k)*v3(i,j,k)
      end if
     end do
    end do
   end do
  end if

  write(60,'(24(1PE16.8e2))') time-inifile, totmass*2d0/msun, remmass/msun*2d0, belmass/msun*2d0, totangmom*2d0, totangmoma*2d0, totangmomb*2d0, totenergy*2d0, totenergya*2d0, totenergyb*2d0, ToWt/totenergy, ToW/totenergya, ToWb/totenergyb, totintenet*2d0, totintene*2d0, totinteneb*2d0, totgrvenet*2d0, totgrvene*2d0, totgrveneb*2d0, totkinenet*2d0, totkinene*2d0, totkineneb*2d0, domega_dt

  call flush(60)
  
return
end subroutine output

