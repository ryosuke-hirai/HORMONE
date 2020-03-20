module amr_evolve_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE AMR_EVOLVE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: For one step of time integration in AMR mode.

subroutine amr_evolve

 use settings,only:rktype
 use grid,only:rungen,dt,time,tn
 use amr_module
 use amr_boundary_mod
 use amr_source_mod
 use amr_rungekutta_mod

 implicit none

 integer dum(1),n

!-----------------------------------------------------------------------------

 allocate( schedule(0:curlvl) )

 schedule = 0.d0

 do ! Loop for one step of lowest level

  dum = minloc(schedule) - 1 ; lvl = dum(1)

  dt_local = dt / dble(cuts**lvl)
  schedule(lvl) = schedule(lvl) + 1.d0/dble(cuts**lvl)
  if(lvl==0)then
   sync = .false.
  elseif(schedule(lvl)==schedule(lvl-1))then
   sync = .true.  ! if child time has synchronized with parent
  else
   sync = .false. ! if child has not synchronized with parent yet
  end if

  do rungen = 1, rktype
   do bkn = 1, totbloks
    if(bk(bkn)%level<=lvl)then
     call amr_boundary
    end if
   end do
   do bkn = 1, totbloks ! second call to update corner boundaries
    if(bk(bkn)%level<=lvl)then
     call amr_boundary
    end if
   end do
   do bkn = 1, totbloks
    if(bk(bkn)%level==lvl)then
!print *,rungen,lvl,sync,bkn
!if(bkn>=1)write(*,'(a2,2(i3),3(1PE13.5e2))'),'b',rungen,bkn,lf(bkn)%flux1(1,1,1,8)*lf(bkn)%detg1(1)*lf(bkn)%idetg1(1)*dt_local,lf(bkn)%e(1:2,1,1)!*lf(bkn)%v1(1:2,1,1)**2.*5.d-1
!print *,bkn,'before',lf(bkn)%v1(:,1,1)
     call amr_numflux(lf(bkn))
     call amr_source(lf(bkn))
     call amr_rungekutta(lf(bkn))
!call amr_output
!if(tn==11)stop
!print *,bkn,'after',lf(bkn)%v1(:,1,1)
    end if
   end do
  end do
!!$print *,'u',lf(16)%u(7:9,1,1,1)
!!$print *,'m',lf(16)%midu1(-2:0,1,1,1)
!!$print *,'o',lf(16)%uorg(7:9,1,1,1)

  if(sync.and.lvl>0)then
   do n = 1, lvl
    dt_local = dt / dble(cuts**(lvl-n))
    if(schedule(lvl)==schedule(lvl-n)) call amr_restriction(lvl-n+1)
   end do
  end if

  if(schedule(curlvl)==1.d0)exit ! End of timestep for lowest level

 end do

 deallocate(schedule)

return
end subroutine amr_evolve

end module amr_evolve_mod
