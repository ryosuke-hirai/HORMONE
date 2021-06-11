!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE AGNDISK
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for explosion in an AGN disk

subroutine agndisk

 use grid
 use settings,only:dt_out
 use physval
 use constants
 use pressure_mod
 use composition_mod,only:get_imu
 use ionization_mod

 implicit none

 real*8:: d_disk, d_amb, Eexp, Mexp, Rexp, Vexp, e_bcg, v0, t0, vedge, v0t0
 real*8:: M_SMBH, r_disk, rgrav, aspratio, hdisk, qkin, mass

!-----------------------------------------------------------------------------

! disk parameters
 d_disk = 1d-11
 d_amb  = 1d-25
 e_bcg  = 1d-30
 M_SMBH = 1d9*msun
 r_disk = 1d3 ! in gravitational radii
 aspratio = 1d-3
! disk composition
 spc(2,:,:,:) = 2.6797147940625970d-1!he4
 spc(3,:,:,:) = 2.4091251279941919d-3!c12
 spc(4,:,:,:) = 7.0571096437237744d-4!n14
 spc(5,:,:,:) = 6.5523144958419283d-3!o16
 spc(6,:,:,:) = 1.4696223348338588d-3!ne20
 spc(7,:,:,:) = 1.0078218370915769d-3!fe56
 spc(8,:,:,:) = 3.353548168440335d-3!other metals
 spc(1,:,:,:) = 1d0-spc(2,:,:,:)-spc(3,:,:,:)-spc(4,:,:,:)-spc(5,:,:,:) &
                   -spc(6,:,:,:)-spc(7,:,:,:)-spc(8,:,:,:)
 
! explosion parameters
 Eexp = 1d51
 Mexp = 1.3d0*msun
 Rexp = 6d12 ! energy injection radius
 qkin = 0.5d0 ! Fraction of energy injected as kinetic energy
 
 do i = is, ie
  if(xi1(i)>Rexp)then
   Rexp = xi1(i)
   exit
  end if
 end do
 
 rgrav = 2d0*G*M_SMBH/clight**2
 r_disk = r_disk*rgrav
 hdisk = r_disk*aspratio
 Vexp = 4d0*pi/3d0*Rexp**3
 v0 = sqrt(qkin*Eexp/(6d0*Mexp))
 vedge = v0*10d0 ! Ejecta velocity at the edge of injection region
 t0 = Rexp/vedge
 v0t0 = Rexp/10d0

! initialize variables
 v1=0d0;v2=0d0;v3=0d0;b1=0d0;b2=0d0;b3=0d0
 e = e_bcg

 imu = 1d0/muconst
 mass=0d0
 do k = ks, ke
  do j = js, je
   do i = is, ie
    imu(i,j,k) = get_imu(spc(1:2,i,j,k))
    p(i,j,k) = eos_p(d(i,j,k),e(i,j,k),T(i,j,k),imu(i,j,k))
    d(i,j,k) = d_amb + d_disk*exp(-0.5d0*(x1(i)*cosc(j)/hdisk)**2)
    if(x1(i)<Rexp)mass=mass+d(i,j,k)*dvol(i,j,k)
! supernova ejecta
    if(x1(i)<Rexp)then
     d   (i,j,k) = Mexp/(8d0*pi)/v0t0**3*exp(-x1(i)/v0t0)
     v1  (i,j,k) = x1(i)*v0/v0t0
     eint(i,j,k) = e(i,j,k) + Eexp*(1d0-qkin)/Mexp*d(i,j,k)
     spc(3,i,j,k) = 0.5d0
     spc(4,i,j,k) = 0.5d0
     spc(1:2,i,j,k) = 0d0
     spc(5:8,i,j,k) = 0d0
     imu(i,j,k) = get_imu(spc(1:2,i,j,k))
     p   (i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k))
    end if

   end do
  end do
 end do
! print*,mass/msun*2d0
! stop

 species(1) = 'h1'
 species(2) = 'he4'
 species(3) = 'c12'
 species(4) = 'n14'
 species(5) = 'o16'
 species(6) = 'ne20'
 species(7) = 'fe56'
 species(8) = 'others'
 
return
end subroutine agndisk

