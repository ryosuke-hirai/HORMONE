module dirichlet_mod

 implicit none

contains
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE DIRICHLETBOUND
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary condition for Dirichlet type.

 subroutine dirichletbound

  use settings,only:mag_on
  use grid
  use ejectamod
  use physval,only:d0,p0,v10,v20,v30,b10,b20,b30,spc0

  real(8):: v0, mej
  integer:: nn, mm

!-----------------------------------------------------------------------------

  if(.not.mag_on)then
   b10 = 0d0
   b20 = 0d0
   b30 = 0d0
  end if

  do j = js,je
   do i = is,ie
    do k = ke+1, ke+2
     if(time*nsdfr(i,j,k)>=tstart.and.time*nsdfr(i,j,k)<=t_ej(count))then
      do nn = 1, count-1
       if(time*nsdfr(i,j,k)>=t_ej(nn).and.time*nsdfr(i,j,k)<t_ej(nn+1))then
        d0 (i,j,k) = (d_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                    - d_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                   / (t_ej(nn+1)-t_ej(nn)) &
                    * nsdfr(i,j,k)**3.d0
!!$        p0 (i,j,k) = (p_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
!!$                    - p_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
!!$                   / (t_ej(nn+1)-t_ej(nn)) &
!!$                   * nsdfr(i,j,k)**(3.d0)
        v0         = (v_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))   &
                    - v_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k)))  &
                   / (t_ej(nn+1)-t_ej(nn))
        v10(i,j,k) = v0 *   nssin(i,j,k)
        v20(i,j,k) = 0d0
        v30(i,j,k) = v0 * (-nscos(i,j,k))
        p0 (i,j,k) = 0.5d0*d0(i,j,k)*(v10(i,j,k)**2+v30(i,j,k)**2)*1d-2
        mej = (m_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
             - m_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
             / (t_ej(nn+1)-t_ej(nn))
        do mm = 1, compsize
         if(mej>comp_ej(0,mm))then
          spc0(2:8,i,j,k) = comp_ej(2:8,mm)
          spc0(1,i,j,k) = 1d0-sum(comp_ej(2:8,mm))
          exit
         end if
        end do
        
       end if
      end do
     elseif(time*nsdfr(i,j,k)<tstart)then !free until the ejecta comes
      call freeboundary(i,k)
     elseif(time*nsdfr(i,j,k)>t_ej(count))then ! extrapolation of data
      d0 (i,j,k) = d_ej(count) * (t_ej(count)/(time*nsdfr(i,j,k)))**2d0 &
                 * nsdfr(i,j,k)**3d0
!!$      p0 (i,j,k) = p_ej(count) * (t_ej(count)/(time*nsdfr(i,j,k)))**4d0 &
!!$                 * nsdfr(i,j,k)**3d0
      v0         = v_ej(count) * (t_ej(count)/(time*nsdfr(i,j,k)))
      v10(i,j,k) = v0 *   nssin(i,j,k)
      v20(i,j,k) = 0d0
      v30(i,j,k) = v0 * (-nscos(i,j,k))
      p0 (i,j,k) = 0.5d0*d0(i,j,k)*(v10(i,j,k)**2+v30(i,j,k)**2)*1d-2
      mej = m_ej(count)
      do mm = 1, compsize
       if(mej>comp_ej(0,mm))then
        spc0(2:8,i,j,k) = comp_ej(2:8,mm)
        spc0(1,i,j,k) = 1d0-sum(comp_ej(2:8,mm))
        exit
       end if
      end do
     else
      print *,'error from dirichletbound.f tn = ',tn
      stop
      exit
     end if

    end do
   end do
  end do

  return
  
end subroutine dirichletbound



!----------------------------------------------------------------------
subroutine outgoingboundary(i_,j_)

  use grid
  use physval

  implicit none

  integer,intent(in):: i_,j_
  integer:: kk
!----------------------------------------------------------------------
  do kk = ks,ke
   d0 (i_,j_,kk) = d(ie,j_,kk)
   p0 (i_,j_,kk) = p(ie,j_,kk)
   v10(i_,j_,kk) = (sign(0.5d0,v1(i_,j,ke))+0.5d0)*v1(ie,j_,kk)
   v20(i_,j_,kk) = (sign(0.5d0,v1(i_,j,ke))+0.5d0)*v2(ie,j_,kk)
   spc0(1:8,i_,j_,kk) = spc(1:8,ie,j_,kk)
  end do

  return
end subroutine outgoingboundary

!----------------------------------------------------------------------
subroutine freeboundary(i_,m_)

  use grid
  use physval

  implicit none

  integer,intent(in):: i_,m_
  integer:: j_
!----------------------------------------------------------------------
  do j_ = js,je
   d0 (i_,j_,m_) = d(i_,j_,ke)
   p0 (i_,j_,m_) = p(i_,j_,ke)
   v10(i_,j_,m_) = v1(i_,j_,ke)
   v30(i_,j_,m_) = (sign(0.5d0,v3(i_,j_,ke))+0.5d0)*v3(i_,j_,ke)
   spc0(1:8,i_,j_,m_) = spc(1:8,i_,j_,ke)
  end do

  return
end subroutine freeboundary

end module dirichlet_mod
