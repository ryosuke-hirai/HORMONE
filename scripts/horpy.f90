module readbin_python

 implicit none

 ! Re-define sink particle properties here because f2py cannot deal with derived types
 integer,allocatable,dimension(:):: sink_i, sink_j, sink_k
 real(8),allocatable,dimension(:,:):: sink_x, sink_v, sink_a, sink_xpol, sink_jspin, sink_jdot, sink_jet_dir
 real(8),allocatable,dimension(:):: sink_mass, sink_mdot, sink_softfac, sink_lsoft, sink_laccr, sink_locres
 real(8),allocatable,dimension(:):: sink_dt, sink_racc, sink_facc, sink_jet_ang, sink_acclum

 interface get_interior
  module procedure get_interior_1D
  module procedure get_interior_3D
  module procedure get_interior_4D
 end interface get_interior

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE SETUP_PYTHON
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set up variables for python interface

subroutine setup_python(dir)

 use settings,only:parafile,include_sinks,nsink
 use physval,only:d
 use setup_mod
 use checksetup_mod
 use mpi_domain,only:domain_decomp
 use allocation_mod

 character(len=*),intent(in):: dir
 character(len=1000):: cwd

!-----------------------------------------------------------------------------

 call getcwd(cwd)
 call chdir(dir)

! Read startfile
 call read_startfile

! Read default parameter file
 call read_default

! Reading parameters
 call read_parameters(parafile)

 call checksetup
 call domain_decomp

! call allocate_subset
 if(allocated(d))call deallocate_all
 call allocations

 if(include_sinks)then
    if(allocated(sink_x))deallocate(sink_i,sink_j,sink_k,sink_x,sink_v,sink_a,sink_xpol,&
                                    sink_jspin,sink_jdot,sink_jet_dir,sink_mass,sink_mdot,sink_softfac,&
                                    sink_lsoft,sink_laccr,sink_locres,sink_dt,sink_racc,sink_facc,&
                                    sink_jet_ang,sink_acclum)
  allocate(sink_x(1:3,1:nsink),sink_mass(1:nsink),sink_i(1:nsink))
  
  !allocate integer variables
  allocate(sink_j,sink_k,mold=sink_i)
  
  !allocate multi-d variables
  allocate(sink_v,sink_a,sink_xpol,sink_jspin,&
           sink_jdot,sink_jet_dir,mold=sink_x)

  !allocate real 1-d variables
  allocate(sink_mdot,sink_softfac,sink_lsoft,&
           sink_laccr,sink_locres,sink_dt,sink_racc,&
           sink_facc,sink_jet_ang,sink_acclum,mold=sink_mass)

 end if

 call chdir(cwd)

return
end subroutine setup_python

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READ_GRIDFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read gridfile in python

subroutine read_gridfile(gridfile)

 use grid
 use readbin_mod
 use metric_mod

 character(len=*),intent(in):: gridfile
 integer :: j

!-----------------------------------------------------------------------------

 call readgrid(gridfile)

 idxi1(gis-1:gie+2)=1d0/dxi1(gis-1:gie+2)
 idx1 (gis-1:gie+2)=1d0/dx1 (gis-1:gie+2)
 idxi2(gjs-1:gje+2)=1d0/dxi2(gjs-1:gje+2)
 idx2 (gjs-1:gje+2)=1d0/dx2 (gjs-1:gje+2)
 idxi3(gks-1:gke+2)=1d0/dxi3(gks-1:gke+2)
 idx3 (gks-1:gke+2)=1d0/dx3 (gks-1:gke+2)

 xi1e = xi1(ie)
 xi1s = xi1(is-1)

 if(crdnt==2)then
  if(allocated(sinc))deallocate(sinc,sini,cosc,cosi)
  allocate( sinc, sini, cosc, cosi, mold=x2 )
  do j = js_global-2, je_global+2
   sinc(j)=real(sin(real(x2 (j),kind=16)),kind=8)!sin0(x2 (j))
   sini(j)=real(sin(real(xi2(j),kind=16)),kind=8)!sin0(xi2(j))
   cosc(j)=real(cos(real(x2 (j),kind=16)),kind=8)!cos0(x2 (j))
   cosi(j)=real(cos(real(xi2(j),kind=16)),kind=8)!cos0(xi2(j))
  end do
 end if

 call metric

 return
end subroutine read_gridfile

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READ_BINFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read a bin file in python

subroutine read_binfile(binfile)

 use settings,only:include_sinks,nsink
 use readbin_mod
 use sink_mod,only:sink

 character(len=*),intent(in):: binfile
 integer:: n

!-----------------------------------------------------------------------------

 call readbin(binfile)

 if(include_sinks)then
   do n = 1, nsink

   !integers
   sink_i(n)   = sink(n)%i
   sink_j(n)   = sink(n)%j
   sink_k(n)   = sink(n)%k

   !reals, 3d
   sink_x(1:3,n)       = sink(n)%x(1:3) 
   sink_v(1:3,n)       = sink(n)%v(1:3)
   sink_a(1:3,n)       = sink(n)%a(1:3)
   sink_xpol(1:3,n)    = sink(n)%xpol(1:3)
   sink_jspin(1:3,n)   = sink(n)%jspin(1:3)
   sink_jdot(1:3,n)    = sink(n)%jdot(1:3)
   sink_jet_dir(1:3,n) = sink(n)%jet_dir(1:3)

   !reals, 1d
   sink_mass(n)    = sink(n)%mass
   sink_mdot(n)    = sink(n)%mdot
   sink_softfac(n) = sink(n)%softfac
   sink_lsoft(n)   = sink(n)%lsoft
   sink_laccr(n)   = sink(n)%laccr
   sink_locres(n)  = sink(n)%locres
   sink_dt(n)      = sink(n)%dt
   sink_racc(n)    = sink(n)%racc
   sink_facc(n)    = sink(n)%facc
   sink_acclum(n)  = sink(n)%acclum
   sink_jet_ang(n) = sink(n)%jet_ang
   
  end do
 end if

return
end subroutine read_binfile

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GET_INTERIOR_1D
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To slice of ghost cells from 1D variables

subroutine get_interior_1D(val,sliced_val)

 real(8),intent(in):: val(:)
 real(8),intent(out):: sliced_val(size(val,1)-4)
 integer:: is,ie

!-----------------------------------------------------------------------------

 is = lbound(val,1)+2
 ie = ubound(val,1)-2

 sliced_val = val(is:ie)

return
end subroutine get_interior_1D

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GET_INTERIOR_3D
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To slice of ghost cells from 3D variables

subroutine get_interior_3D(val,sliced_val)

 real(8),intent(in):: val(:,:,:)
 real(8),intent(out):: sliced_val(size(val,1)-4,size(val,2)-4,size(val,3)-4)
 integer:: is,ie,js,je,ks,ke

!-----------------------------------------------------------------------------

 is = lbound(val,1)+2
 ie = ubound(val,1)-2
 js = lbound(val,2)+2
 je = ubound(val,2)-2
 ks = lbound(val,3)+2
 ke = ubound(val,3)-2

 sliced_val = val(is:ie,js:je,ks:ke)

return
end subroutine get_interior_3D

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GET_INTERIOR_4D
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To slice of ghost cells from 3D variables

subroutine get_interior_4D(val,sliced_val)

 real(8),intent(in):: val(:,:,:,:)
 real(8),intent(out):: sliced_val(size(val,1),size(val,2)-4,size(val,3)-4,size(val,4)-4)
 integer:: is,ie,js,je,ks,ke,ls,le

!-----------------------------------------------------------------------------

 is = lbound(val,1)
 ie = ubound(val,1)
 js = lbound(val,2)+2
 je = ubound(val,2)-2
 ks = lbound(val,3)+2
 ke = ubound(val,3)-2
 ls = lbound(val,4)+2
 le = ubound(val,4)-2

 sliced_val = val(is:ie,js:je,ks:ke,ls:le)

return
end subroutine get_interior_4D

end module readbin_python
