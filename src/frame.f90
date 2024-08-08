module frame_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SET_FRAME_ACC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get the frame acceleration in case of noninertial frame

 subroutine set_frame_acc

  use settings,only:frame

!-----------------------------------------------------------------------------

  select case(frame)
  case(0) ! inertial frame
   return
  case(1:) ! centre on sink
   call centre_on_sink(frame)
  case default
   print*,'Error in value of frame; frame=',frame
   stop
  end select

  return
 end subroutine set_frame_acc


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE CENTRE_ON_SINK
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set the frame to centre on a sink particle

subroutine centre_on_sink(n)

 use grid,only:frame_acc
 use sink_mod,only:sink,get_sink_loc,get_sinkgas_acc,get_sinksink_acc

 integer,intent(in):: n

!-----------------------------------------------------------------------------

 call get_sink_loc(sink(n))
 call get_sinkgas_acc(sink(n))
 call get_sinksink_acc(sink)

 frame_acc = -sink(n)%a

return
end subroutine centre_on_sink

end module frame_mod
