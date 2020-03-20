!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE DIRICHLETBOUND
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary condition for Dirichlet type.

subroutine dirichletbound

  use grid
  use dirichlet
  use ejectamod
  use physval,only:gamma

  implicit none

  real*8 v0
  integer nn, mm

!-----------------------------------------------------------------------------

  do j = js,je
   do i = is,ie
    do mm = 1,2
     k = ke+mm
     if(time*nsdfr(i,j,k)>=tstart.and.time*nsdfr(i,j,k)<=t_ej(count))then
      do nn = 1, count-1
       if(time*nsdfr(i,j,k)>=t_ej(nn).and.time*nsdfr(i,j,k)<t_ej(nn+1))then
        d0 (i,j,k) = (d_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                    - d_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                   / (t_ej(nn+1)-t_ej(nn)) &
                    * nsdfr(i,j,k)**3.d0
        p0 (i,j,k) = (p_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                    - p_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                   / (t_ej(nn+1)-t_ej(nn)) &
                   * nsdfr(i,j,k)**(3.d0)
        v0         = (v_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))   &
                    - v_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k)))  &
                   / (t_ej(nn+1)-t_ej(nn))
        v10(i,j,k) = v0 *   nssin(i,j,k)
        v30(i,j,k) = v0 * (-nscos(i,j,k))

       end if
      end do
     elseif(time*nsdfr(i,j,k)<tstart)then !free until the ejecta comes
      call freeboundary(i,k)
!      print *,i,k,j
     elseif(time*nsdfr(i,j,k)>t_ej(count))then ! extrapolation of data
      d0 (i,j,k) = d_ej(count) * (t_ej(count)/(time*nsdfr(i,j,k)))**2d0 &
                 * nsdfr(i,j,k)**3d0
      p0 (i,j,k) = p_ej(count) * (t_ej(count)/(time*nsdfr(i,j,k)))**4d0 &
                 * nsdfr(i,j,k)**3d0
      v0         = v_ej(count) * (t_ej(count)/(time*nsdfr(i,j,k)))
      v10(i,j,k) = v0 *   nssin(i,j,k)
      v30(i,j,k) = v0 * (-nscos(i,j,k))
     else
      print *,'error from dirichletbound.f tn = ',tn
      stop
      exit
     end if

    end do
   end do
  end do


  return
contains
!----------------------------------------------------------------------
subroutine freeboundary(i_,m_)

  use grid
  use dirichlet
  use physval

  implicit none

  integer,intent(in):: i_,m_
  integer j_
!----------------------------------------------------------------------
  do j_ = js,je
   d0 (i_,j_,m_) = d(i_,j_,ke)
   p0 (i_,j_,m_) = p(i_,j_,ke)
   v10(i_,j_,m_) = v1(i_,j_,ke)
   v30(i_,j_,m_) = (sign(0.5d0,v3(i_,j_,ke))+0.5d0)*v3(i_,j_,ke)
  end do

  return
end subroutine freeboundary

end subroutine dirichletbound
