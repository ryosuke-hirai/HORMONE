!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE EOSTEST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To test the EoS module forwards and backwards consistency

subroutine eostest

 use grid
 use physval
 use constants
 use pressure_mod
 use composition_mod,only:get_imu
 use ionization_mod

 implicit none

 real*8:: dtest, dtesti, dtestf, ptest, ptesti, ptestf, etest, Ttest, imutest
 real*8:: Xtest, Ytest, Ztest, Xtesti, Xtestf, ptest2, Ttest2, Ttesti, Ttestf
 real*8:: Qtesti, Qtestf, Qtest, ere, imutest2, Ttest0, cftest
 integer:: iie, jje, kke

!-----------------------------------------------------------------------------

 iie =0
 jje = 100
 kke = 100

 dtesti = -10d0
 dtestf = 0d0
 ptesti = 2d0
 ptestf = 12d0
 Xtesti = 0.6d0
 Xtestf = 0.8d0
 Ztest = 0.02d0

 Ttesti = 3d0
 Ttestf = 6d0
 Qtesti = -8d0
 Qtestf = 0d0

 Xtesti=0.7d0

 open(unit=123,file='eostest.dat',status='replace')
 
 do k = 0, kke
  Ttest0 = 10d0**(Ttesti+(Ttestf-Ttesti)/dble(kke)*dble(k))
  do j = 0, jje
   Qtest = 10d0**(Qtesti+(Qtestf-Qtesti)/dble(jje)*dble(j))
   Xtest = Xtesti!+(Xtestf-Xtesti)/dble(iie)*dble(i)
   Ytest = 1d0-Xtest-Ztest
   Ttest = Ttest0

   dtest = 10d0**(log10(Qtest)+2d0*log10(Ttest)-12d0)
   call get_erec_imurec(log10(dtest),Ttest,Xtest,Ytest,ere,imutest)
!   imutest = 1d0/0.62d0

   ptest = fac_pgas*imutest*dtest*Ttest+arad*Ttest**4/3d0
   Ttest = 1d3;Ttest2=1d3
   etest = eos_e(dtest,ptest,Ttest,imutest,Xtest,Ytest)
   ptest2 = eos_p(dtest,etest,Ttest2,imutest2,Xtest,Ytest)
   call eos_p_cf(dtest,0d0,0d0,0d0,0d0,0d0,0d0,etest,Ttest3,imutest3,ptest3,cftest,Xtest,Ytest)
   write(123,'(2i5,16(1PE23.15e2))')&
    j,k,dtest,ptest,Qtest,etest,Xtest,&
    1d0/imutest,1d0/imutest2,1d0/imutest3,&
    Ttest,Ttest2,Ttest3,ptest2,ptest3,&
    ptest2/ptest-1d0,ptest3/ptest-1d0
  end do
  write(123,'()')
 end do

 close(123)

 stop
 
return
end subroutine eostest

