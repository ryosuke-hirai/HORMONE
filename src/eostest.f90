module eostest_mod
 implicit none

contains
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EOSTEST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To test the EoS module forwards and backwards consistency

subroutine eostest

 use settings,only:compswitch,eoserr
 use grid
 use pressure_mod
 use composition_mod,only:get_imu
 use ionization_mod

 real*8:: d, p, e, T, imu
 real*8:: X, Y, Z, Xi, Xf, p2, T2, Ti, Tf
 real*8:: Qi, Qf, Q, erec, imu2, imu3, T0, cs
 real*8:: T3, p3, p4, d2, S, T4
 real*8:: rerrp1,rerrp2,rerrp3,rerrT1,rerrT2,rerrd1
 integer:: iie, jje, kke, ierr, ui
 character(len=100):: form1

!-----------------------------------------------------------------------------

! Set grid resolution for EoS unit test
 iie = 0
 jje = 100
 kke = 100
 Z = 0.02d0
 
! Set range of grid for the EoS test
 Xi = 0.6d0
 Xf = 0.8d0
 Ti = 3d0
 Tf = 6d0
 Qi = -8d0
 Qf = 0d0

 Xi=0.7d0

! Test gas + radiation EoS
 print'(a,1PE9.1e2)','Tolerance on temperature is set to eoserr = ',eoserr
 print*,''
 print*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 print*,''
 print*,'Testing ideal gas EoS forwards and backwards consistency...'
 eostype = 0
 compswitch = 0
 imu = 1d0/0.62d0
 rerrp1=0d0;rerrp2=0d0;rerrT1=0d0;rerrT2=0d0;rerrd1=0d0
 open(newunit=ui,file='data/eostest_gas.dat',status='replace')
 write(ui,'(2a5,15a23)')&
    'j','k','d','Q','e',&
    'T_from_eos_e','T_from_eos_p','T_from_eos_p_cs',&
    'p_for_eos_e','p_from_eos_p','p_from_eos_p_cs',&
    'rel_error_for_eos_p','rel_error_for_eos_p_cs',&
    'rel_error_p_from_ds','entropy_from_dp','d_from_inventropy'

 do k = 0, kke ! loop over T
  T0 = 10d0**(Ti+(Tf-Ti)/dble(kke)*dble(k))
  do j = 0, jje ! loop over Q
   Q = 10d0**(Qi+(Qf-Qi)/dble(jje)*dble(j))
   X = Xi!+(Xf-Xi)/dble(iie)*dble(i)
   Y = 1d0-X-Z
   T = T0

   d = 10d0**(log10(Q)+2d0*log10(T)-12d0)

   p = fac_pgas*imu*d*T
   T = 1d3; T2 = 1d3; T3 = 1d3; T4 = 1d3
   e = eos_e(d,p,T,imu)
   p2 = eos_p(d,e,T2,imu)
   call eos_p_cs(d,e,T3,imu,p3,cs,ierr=ierr)

   S = entropy_from_dp(d,p,T4,imu)
   d2 = get_d_from_ps(p,S,imu)
   p4 = get_p_from_ds(d,S,imu)

   rerrp1 = max(rerrp1,abs(p2/p-1d0))
   rerrp2 = max(rerrp2,abs(p3/p-1d0))
   rerrp3 = max(rerrp3,abs(p4/p-1d0))
   rerrT1 = max(rerrT1,abs(T2/T-1d0))
   rerrT2 = max(rerrT2,abs(T3/T-1d0))
   rerrd1 = max(rerrd1,abs(d2/d-1d0))
   write(ui,'(2i5,15(1PE23.15e2))')&
    j,k,d,Q,e,&
    T,T2,T3,p,p2,p3,&
    p2/p-1d0,p3/p-1d0,&
    p4/p-1d0,S,d2
   
  end do
  write(ui,'()')
 end do

 close(ui)

 print*,'Relative errors computed from each module is:'
 form1='(a,1PE10.3e2,2X,a,1PE10.3e2)'
 print form1,'eos_p    : temperature =',rerrT1,'pressure =',rerrp1
 print form1,'eos_p_cs : temperature =',rerrT2,'pressure =',rerrp2
 print form1,'get_d_from_ps: density =',rerrd1
 print form1,'get_p_from_ds: pressure =',rerrp3
 print*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 print*,''


 print*,'Testing gas+rad EoS forwards and backwards consistency...'
 eostype = 1
 compswitch = 1
 imu = 1d0/0.62d0
 rerrp1=0d0;rerrp2=0d0;rerrT1=0d0;rerrT2=0d0;rerrd1=0d0
 open(newunit=ui,file='data/eostest_gasrad.dat',status='replace')
 write(ui,'(2a5,14a23)')&
    'j','k','d','Q','e',&
    'T_from_eos_e','T_from_eos_p','T_from_eos_p_cs',&
    'p_for_eos_e','p_from_eos_p','p_from_eos_p_cs',&
    'rel_error_for_eos_p','rel_error_for_eos_p_cs',&
    'rel_error_p_from_ds','entropy_from_dp','d_from_inventropy'

 do k = 0, kke ! loop over T
  T0 = 10d0**(Ti+(Tf-Ti)/dble(kke)*dble(k))
  do j = 0, jje ! loop over Q
   Q = 10d0**(Qi+(Qf-Qi)/dble(jje)*dble(j))
   X = Xi!+(Xf-Xi)/dble(iie)*dble(i)
   Y = 1d0-X-Z
   T = T0

   d = 10d0**(log10(Q)+2d0*log10(T)-12d0)

   p = fac_pgas*imu*d*T+arad*T**4/3d0
   T = 1d3; T2 = 1d3; T3 = 1d3; T4 = 1d3
   e = eos_e(d,p,T,imu)
   p2 = eos_p(d,e,T2,imu)
   call eos_p_cs(d,e,T3,imu,p3,cs,ierr=ierr)

   S = entropy_from_dp(d,p,T4,imu)
   d2 = get_d_from_ps(p,S,imu)
   p4 = get_p_from_ds(d,S,imu)

   rerrp1 = max(rerrp1,abs(p2/p-1d0))
   rerrp2 = max(rerrp2,abs(p3/p-1d0))
   rerrp3 = max(rerrp3,abs(p4/p-1d0))
   rerrT1 = max(rerrT1,abs(T2/T-1d0))
   rerrT2 = max(rerrT2,abs(T3/T-1d0))
   rerrd1 = max(rerrd1,abs(d2/d-1d0))
   write(ui,'(2i5,15(1PE23.15e2))')&
    j,k,d,Q,e,&
    T,T2,T3,p,p2,p3,&
    p2/p-1d0,p3/p-1d0,&
    p4/p-1d0,S,d2
   
  end do
  write(ui,'()')
 end do

 close(ui)

 print*,'Relative errors computed from each module is:'
 form1='(a,1PE10.3e2,2X,a,1PE10.3e2)'
 print form1,'eos_p    : temperature =',rerrT1,'pressure =',rerrp1
 print form1,'eos_p_cs : temperature =',rerrT2,'pressure =',rerrp2
 print form1,'get_d_from_ps: density =',rerrd1
 print form1,'get_p_from_ds: pressure =',rerrp3
 print*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 print*,''
 
! Test gas + radiation + recombination EoS
 print*,'Testing gas+rad+rec EoS forwards and backwards consistency...'
 eostype = 2
 compswitch = 2
 call ionization_setup
 rerrp1=0d0;rerrp2=0d0;rerrT1=0d0;rerrT2=0d0;rerrd1=0d0
 open(newunit=ui,file='data/eostest_gasradrec.dat',status='replace')
 write(ui,'(2a5,18a23)')&
    'j','k','d','Q','e','X',&
    'mu_from_eos_e','mu_from_eos_p','mu_from_eos_p_cs',&
    'T_from_eos_e','T_from_eos_p','T_from_eos_p_cs',&
    'p_for_eos_e','p_from_eos_p','p_from_eos_p_cs',&
    'rel_error_for_eos_p','rel_error_for_eos_p_cs',&
    'rel_error_p_from_ds','entropy_from_dp','d_from_inventropy'

 do k = 0, kke ! loop over T
  T0 = 10d0**(Ti+(Tf-Ti)/dble(kke)*dble(k))
  do j = 0, jje ! loop over Q
   Q = 10d0**(Qi+(Qf-Qi)/dble(jje)*dble(j))
   X = Xi!+(Xf-Xi)/dble(iie)*dble(i)
   Y = 1d0-X-Z
   T = T0

   d = 10d0**(log10(Q)+2d0*log10(T)-12d0)
   call get_erec_imurec(log10(d),T,X,Y,erec,imu)

   p = fac_pgas*imu*d*T+arad*T**4/3d0
   T = 1d3; T2 = 1d3; T3 = 1d3; T4 = 1d3
   e = eos_e(d,p,T,imu,X,Y)
   p2 = eos_p(d,e,T2,imu2,X,Y)
   call eos_p_cs(d,e,T3,imu3,p3,cs,X,Y,ierr)

! Currently not ready for entropy calculations
!   S = entropy_from_dp(d,p,T4,imu4,X,Y)
!   d2 = get_d_from_ps(p,S,imu4,X,Y)
!   p4 = get_p_from_ds(d,S,imu)


   rerrp1 = max(rerrp1,abs(p2/p-1d0))
   rerrp2 = max(rerrp2,abs(p3/p-1d0))
   rerrp3 = max(rerrp3,abs(p4/p-1d0))
   rerrT1 = max(rerrT1,abs(T2/T-1d0))
   rerrT2 = max(rerrT2,abs(T3/T-1d0))
   rerrd1 = max(rerrd1,abs(d2/d-1d0))
   write(ui,'(2i5,18(1PE23.15e2))')&
    j,k,d,Q,e,X,&
    1d0/imu,1d0/imu2,1d0/imu3,&
    T,T2,T3,p,p2,p3,&
    p2/p-1d0,p3/p-1d0,&
    p4/p-1d0,S,d2
   
  end do
  write(ui,'()')
 end do

 close(ui)

 print*,'Relative errors computed from each module is:'
 print form1,'eos_p    : temperature =',rerrT1,'pressure =',rerrp1
 print form1,'eos_p_cs : temperature =',rerrT2,'pressure =',rerrp2
! print form1,'get_d_from_ps: density =',rerrd1
! print form1,'get_p_from_ds: pressure =',rerrp3
 print*,'Entropy calculations for eostype=2 currently under construction'

 stop
 
return
end subroutine eostest

end module eostest_mod
