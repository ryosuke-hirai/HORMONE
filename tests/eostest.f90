!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE EOSTEST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To test the EoS module forwards and backwards consistency

subroutine eostest

 use settings,only:compswitch
 use grid
 use pressure_mod
 use composition_mod,only:get_imu
 use ionization_mod

 implicit none

 real*8:: d, di, df, p, pi, pf, e, T, imu
 real*8:: X, Y, Z, Xi, Xf, p2, T2, Ti, Tf
 real*8:: Qi, Qf, Q, erec, imu2, imu3, T0, cf
 real*8:: T3, p3
 integer:: iie, jje, kke, ierr, ui

!-----------------------------------------------------------------------------
 
 iie =0
 jje = 100
 kke = 100

 di = -10d0
 df = 0d0
 pi = 2d0
 pf = 12d0
 Xi = 0.6d0
 Xf = 0.8d0
 Z = 0.02d0

 Ti = 3d0
 Tf = 6d0
 Qi = -8d0
 Qf = 0d0

 Xi=0.7d0

 eostype = 2
 compswitch = 2
 call ionization_setup
 open(newunit=ui,file='data/eos.dat',status='replace')
 write(ui,'(2a5,16a23)')&
    'j','k','d','Q','e','X',&
    'mu_from_eos_e','mu_from_eos_p','mu_from_eos_p_cf',&
    'T_from_eos_e','T_from_eos_p','T_from_eos_p_cf',&
    'p_for_eos_e','p_from_eos_p','p_from_eos_p_cf',&
    'rel_error_for_eos_p','rel_error_for_eos_p_cf'
 
 do k = 0, kke ! loop over T
  T0 = 10d0**(Ti+(Tf-Ti)/dble(kke)*dble(k))
  do j = 0, jje ! loop over Q
   Q = 10d0**(Qi+(Qf-Qi)/dble(jje)*dble(j))
   X = Xi!+(Xf-Xi)/dble(iie)*dble(i)
   Y = 1d0-X-Z
   T = T0

   d = 10d0**(log10(Q)+2d0*log10(T)-12d0)
   call get_erec_imurec(log10(d),T,X,Y,erec,imu)
!   imu = 1d0/0.62d0

   p = fac_pgas*imu*d*T+arad*T**4/3d0
   T = 1d3; T2 = 1d3; T3 = 1d3
   e = eos_e(d,p,T,imu,X,Y)
   p2 = eos_p(d,e,T2,imu2,X,Y)
   call eos_p_cf(d,0d0,0d0,0d0,e,T3,imu3,p3,cf,X,Y,ierr)
   write(ui,'(2i5,16(1PE23.15e2))')&
    j,k,d,Q,e,X,&
    1d0/imu,1d0/imu2,1d0/imu3,&
    T,T2,T3,p,p2,p3,&
    p2/p-1d0,p3/p-1d0
  end do
  write(ui,'()')
 end do

 close(ui)
 
return
end subroutine eostest

