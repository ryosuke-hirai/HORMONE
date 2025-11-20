module output_mod
 implicit none

 integer:: ievo,iskf
 public:: output,terminal_output,set_file_name,write_extgrv,evo_output,&
          scaling_output,write_grid,write_plt,write_grid_dat
 private:: write_bin,get_header,add_column, add_directory_to_filename,&
           write_val,write_vertical_slice

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output snapshot data

subroutine output

 use settings,only:is_test,output_ascii
 use grid,only:tn
 use profiler_mod
 use mpi_utils,only:barrier_mpi

 integer:: wtind

!----------------------------------------------------------------------------

 if(is_test) return ! do not bother outputting anything if it is a test

 if(tn==0)then
  wtind = wtou1 ! Record time of initial output separately
 else
  wtind = wtout
 end if

 call start_clock(wtind)

 if(tn==0)call write_grid

 call write_bin
 if(output_ascii)call write_plt

! Tracer particle outputs
 call write_bpt
 call write_ptc

 call evo_output
 call sink_output

 call stop_clock(wtind)

 call profiler_output

 call barrier_mpi

 return
end subroutine output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE TERMINAL_OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Output things on the terminal

subroutine terminal_output

 use mpi_utils,only:myrank
 use settings,only:dt_unit_in_sec,dt_unit
 use grid,only:tn,time,dt

!-----------------------------------------------------------------------------

 if (myrank==0) then
  print'(a,i10,2(3X,a,1PE13.5e2,1X,a))',&
    'tn =',tn,&
    'time =',time/dt_unit_in_sec,dt_unit,&
    'dt =',dt/dt_unit_in_sec,dt_unit
 end if

return
end subroutine terminal_output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OPEN_EVOFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Start or open evofile

subroutine open_evofile

 use settings,only:crdnt,sigfig,gravswitch,mag_on,crdnt,write_evo,start
 use grid,only:dim
 use mpi_utils, only:myrank

 integer::ierr
 character(len=70):: forma, filename

!-----------------------------------------------------------------------------

 if (myrank/=0) return
 if(.not.write_evo)return

 write(forma,'("(a",i2,")")')sigfig+8 ! for strings
 call add_directory_to_filename('evo.dat',filename)

 if(start/=0)then
  open(newunit=ievo,file=filename,status='old',&
       position='append',iostat=ierr)
  if(ierr==0)return ! Return if evofile already exists.
 end if

! Write headers if it needs to be freshly made.
 open(newunit=ievo,file=filename,status='replace')

 write(ievo,'(a10)',advance='no')'tn'
 write(ievo,forma,advance="no")'time'
 write(ievo,forma,advance="no")'tot_mass'
 write(ievo,forma,advance="no")'tot_e'
 write(ievo,forma,advance="no")'tot_eint'
 write(ievo,forma,advance="no")'tot_ekin'
 if(gravswitch>=1)write(ievo,forma,advance="no")'tot_egrv'
 if(gravswitch>=1)write(ievo,forma,advance="no")'bound_mass'
 if(gravswitch>=1)write(ievo,forma,advance="no")'bound_e'
 if(gravswitch>=1)write(ievo,forma,advance="no")'bound_angmom'
 if(mag_on)write(ievo,forma,advance="no")'tot_emag'
 if(dim>=2.and.crdnt>=1)write(ievo,forma,advance="no")'tot_angmom'
 write(ievo,'()')

 return
end subroutine open_evofile



!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE EVO_OUTPUT
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output quantities as a function of time

subroutine evo_output

 use settings
 use grid
 use physval
 use gravmod
 use mpi_utils, only:myrank, allreduce_mpi

 implicit none

 integer:: i,j,k
 character(len=50):: forme
 real(8):: Mtot, Etot, Eitot, Ektot, Egtot, Ebtot, Jtot, Mbound, Ebound, Jbound

!-----------------------------------------------------------------------------

 if(.not.write_evo)return

 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers

 Mtot  = sum(d   (is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))
 Etot  = sum(e   (is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))
 Eitot = sum(eint(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))
 Ektot = 0.5d0*sum(d(is:ie,js:je,ks:ke)&
                   * ( v1(is:ie,js:je,ks:ke)**2 &
                     + v2(is:ie,js:je,ks:ke)**2 &
                     + v3(is:ie,js:je,ks:ke)**2 ) &
                   * dvol(is:ie,js:je,ks:ke) )
 if(eq_sym)then
  Mtot  = 2d0*Mtot
  Etot  = 2d0*Etot
  Eitot = 2d0*Eitot
  Ektot = 2d0*Ektot
 end if

 if(mag_on)then
  Ebtot = 0.5d0*sum( ( b1(is:ie,js:je,ks:ke)**2 &
                     + b2(is:ie,js:je,ks:ke)**2 &
                     + b3(is:ie,js:je,ks:ke)**2 ) &
                   * dvol(is:ie,js:je,ks:ke) )
  if(eq_sym) Ebtot = 2d0*Ebtot
 end if
 if(gravswitch>=1)then
! grvphi is self-gravity so is double counted, whereas extgrv is not
  Egtot = sum(d(is:ie,js:je,ks:ke) &
              *(totphi(is:ie,js:je,ks:ke)-0.5d0*grvphi(is:ie,js:je,ks:ke)) &
              *dvol(is:ie,js:je,ks:ke))
  if (is==is_global) Mtot = Mtot + mc(is-1)

  Mbound=0d0;Ebound=0d0;Jbound=0d0
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(totphi(i,j,k)+(e(i,j,k)-eint(i,j,k))/d(i,j,k)<=0d0)then
      Mbound = Mbound + d(i,j,k)*dvol(i,j,k)
      Ebound = Ebound + (e(i,j,k)+d(i,j,k)*(totphi(i,j,k)-0.5d0*grvphi(i,j,k)))&
                        *dvol(i,j,k)
      select case(crdnt)
      case(1)
       Jbound = Jbound + d(i,j,k)*x1(i)*v2(i,j,k)*dvol(i,j,k)
      case(2)
       Jbound = Jbound + d(i,j,k)*x1(i)*sinc(j)*v3(i,j,k)*dvol(i,j,k)
      end select
     end if
    end do
   end do
  end do

  if(eq_sym)then
   Egtot  = 2d0*Egtot
   Mbound = 2d0*Mbound
   Ebound = 2d0*Ebound
   Jbound = 2d0*Jbound
  end if
  if (include_extgrv .and. is==is_global) Mbound = Mbound+mc(is-1)
 end if

 if(dim>=2.and.crdnt>=1)then
  Jtot = 0d0
  select case(crdnt)
  case(1)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      Jtot = Jtot + d(i,j,k)*x1(i)*v2(i,j,k)*dvol(i,j,k)
     end do
    end do
   end do
  case(2)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      Jtot = Jtot + d(i,j,k)*x1(i)*sinc(j)*v3(i,j,k)*dvol(i,j,k)
     end do
    end do
   end do
  end select
  if(eq_sym) Jtot = 2d0*Jtot
 end if

 call allreduce_mpi('sum', Mtot)
 call allreduce_mpi('sum', Etot)
 call allreduce_mpi('sum', Eitot)
 call allreduce_mpi('sum', Ektot)
 call allreduce_mpi('sum', Egtot)
 call allreduce_mpi('sum', Mbound)
 call allreduce_mpi('sum', Ebound)
 call allreduce_mpi('sum', Jbound)
 call allreduce_mpi('sum', Ebtot)
 call allreduce_mpi('sum', Jtot)

 if (myrank==0) then
   write(ievo,'(i10)',advance='no')tn
   call write_anyval(ievo,forme,time,1)
   call write_anyval(ievo,forme,Mtot,1)
   call write_anyval(ievo,forme,Etot,1)
   call write_anyval(ievo,forme,Eitot,1)
   call write_anyval(ievo,forme,Ektot,1)
   if(gravswitch>=1)call write_anyval(ievo,forme,Egtot ,1)
   if(gravswitch>=1)call write_anyval(ievo,forme,Mbound,1)
   if(gravswitch>=1)call write_anyval(ievo,forme,Ebound,1)
   if(gravswitch>=1)call write_anyval(ievo,forme,Jbound,1)
   if(mag_on)call write_anyval(ievo,forme,Ebtot,1)
   if(dim>=2.and.crdnt>=1)call write_anyval(ievo,forme,Jtot,1)
   write(ievo,'()')
   flush(ievo)
 end if

return
end subroutine evo_output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OPEN_SINKFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Start or open sinkfile

subroutine open_sinkfile

 use settings,only:include_sinks,sigfig,start,include_accretion,frame
 use constants,only:msun
 use sink_mod,only:nsink,sink
 use mpi_utils, only:myrank

 integer:: ierr,n
 character(len=70):: forma, forme, header, filename

!-----------------------------------------------------------------------------

 if (myrank/=0) return
 if(.not.include_sinks)return

 write(forma,'("(a",i2,")")')sigfig+8 ! for strings
 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers
 call add_directory_to_filename('sinks.dat',filename)

 if(start/=0)then
  open(newunit=iskf,file=filename,status='old',&
       position='append',iostat=ierr)
  if(ierr==0)return ! Return if sinks file already exists.
 end if

! Write headers if it needs to be freshly made.
 open(newunit=iskf,file=filename,status='replace')

 do n = 1, nsink
  write(iskf,'(2x,a,i0,a)',advance="no")"Msink_",n,"/Msun="
  write(iskf,forme,advance="no")sink(n)%mass/msun
 end do
 write(iskf,'()')

 write(iskf,'(a10)',advance='no')'tn'
 write(iskf,forma,advance="no")'time'
 if(frame>0)then
  write(iskf,forma,advance="no")'frame_acc_x1'
  write(iskf,forma,advance="no")'frame_acc_x2'
  write(iskf,forma,advance="no")'frame_acc_x3'
  write(iskf,forma,advance="no")'frame_vel_x1'
  write(iskf,forma,advance="no")'frame_vel_x2'
  write(iskf,forma,advance="no")'frame_vel_x3'
  write(iskf,forma,advance="no")'frame_pos_x1'
  write(iskf,forma,advance="no")'frame_pos_x2'
  write(iskf,forma,advance="no")'frame_pos_x3'
 end if
 do n = 1, nsink
  write(header,'("sink_",i0)')n
  write(iskf,forma,advance="no")trim(header)//'_x1'
  write(iskf,forma,advance="no")trim(header)//'_x2'
  write(iskf,forma,advance="no")trim(header)//'_x3'
  write(iskf,forma,advance="no")trim(header)//'_v1'
  write(iskf,forma,advance="no")trim(header)//'_v2'
  write(iskf,forma,advance="no")trim(header)//'_v3'
  if(include_accretion.and.sink(n)%laccr>0d0)then
   write(iskf,forma,advance="no")trim(header)//'_mass'
   write(iskf,forma,advance="no")trim(header)//'_mdot'
   write(iskf,forma,advance="no")trim(header)//'_facc'
   write(iskf,forma,advance="no")trim(header)//'_Lacc'
   write(iskf,forma,advance="no")trim(header)//'_Jspin'
   write(iskf,forma,advance="no")trim(header)//'_jdot'
  end if
 end do
 write(iskf,'()')
 flush(iskf)

 return
end subroutine open_sinkfile

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE SINK_OUTPUT
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output sink properties as a function of time

subroutine sink_output

 use constants, only: msun,year
 use settings,  only: sigfig,include_sinks,include_accretion,frame
 use grid,      only: tn,time,frame_acc,frame_vel,frame_pos
 use sink_mod,  only: nsink,sink
 use mpi_utils, only: myrank

 integer:: n
 character(len=50):: forme

!-----------------------------------------------------------------------------

 if (myrank/=0) return
 if(.not.include_sinks)return

 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers

 write(iskf,'(i10)',advance='no')tn
 call write_anyval(iskf,forme,time,1)
 if(frame>0)then
  call write_anyval(iskf,forme,frame_acc(1),1)
  call write_anyval(iskf,forme,frame_acc(2),1)
  call write_anyval(iskf,forme,frame_acc(3),1)
  call write_anyval(iskf,forme,frame_vel(1),1)
  call write_anyval(iskf,forme,frame_vel(2),1)
  call write_anyval(iskf,forme,frame_vel(3),1)
  call write_anyval(iskf,forme,frame_pos(1),1)
  call write_anyval(iskf,forme,frame_pos(2),1)
  call write_anyval(iskf,forme,frame_pos(3),1)
 end if
 do n = 1, nsink
  call write_anyval(iskf,forme,sink(n)%x(1),1)
  call write_anyval(iskf,forme,sink(n)%x(2),1)
  call write_anyval(iskf,forme,sink(n)%x(3),1)
  call write_anyval(iskf,forme,sink(n)%v(1),1)
  call write_anyval(iskf,forme,sink(n)%v(2),1)
  call write_anyval(iskf,forme,sink(n)%v(3),1)
  if(include_accretion.and.sink(n)%laccr>0d0)then
   call write_anyval(iskf,forme,sink(n)%mass/msun,1)
   call write_anyval(iskf,forme,sink(n)%mdot/msun*year,1)
   call write_anyval(iskf,forme,sink(n)%facc,1)
   call write_anyval(iskf,forme,sink(n)%acclum,1)
   call write_anyval(iskf,forme,sink(n)%Jspin(3),1)
   call write_anyval(iskf,forme,sink(n)%jdot(3),1)
  end if
 end do
 write(iskf,'()')
 flush(iskf)

return
end subroutine sink_output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE WRITE_GRID
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output gridfile.bin and gridfile.dat
subroutine write_grid
 use mpi_utils, only:myrank,barrier_mpi
 use settings,only:output_ascii

 if (myrank==0) &
  call write_grid_bin
 call barrier_mpi
 if(output_ascii)call write_grid_dat
 call barrier_mpi

 return
end subroutine write_grid

subroutine write_grid_bin
 use grid
 integer :: ui
 character(len=70):: filename

 call add_directory_to_filename('gridfile.bin',filename)

 open(newunit=ui,file=filename,status='replace',form='unformatted')

 write(ui)x1 (gis_global-2:gie_global+2), xi1 (gis_global-2:gie_global+2), &
          dx1(gis_global-2:gie_global+2), dxi1(gis_global-2:gie_global+2), &
          x2 (gjs_global-2:gje_global+2), xi2 (gjs_global-2:gje_global+2), &
          dx2(gjs_global-2:gje_global+2), dxi2(gjs_global-2:gje_global+2), &
          x3 (gks_global-2:gke_global+2), xi3 (gks_global-2:gke_global+2), &
          dx3(gks_global-2:gke_global+2), dxi3(gks_global-2:gke_global+2)

 close(ui)

 print*,"Outputted: ",'gridfile.bin'

end subroutine write_grid_bin

subroutine write_grid_dat
 use settings
 use grid
 use mpi_utils, only: myrank
 use mpi_domain, only: is_my_domain
 use io, only: write_string_master, open_file_write_ascii, close_file

 character(len=50):: formhead,formval,formnum
 character(len=200) :: str,filename
 integer:: i,j,k,ui

 call add_directory_to_filename('gridfile.dat',filename)
 call open_file_write_ascii(filename,ui)
 call write_string_master(ui, '')
 write(formhead,'("(",i1,"a5,",i1,"a",i2,")")')dim,dim+1,sigfig+8
 write(formval ,'("(",i1,"i5,",i1,"(1x,1PE",i2,".",i2,"e2))")')&
     dim,dim+1,sigfig+7,sigfig-1

 select case (dim)
! 1D outputs
 case(1)
  write(formnum,'("(",a4,"i4,2i",i2,")")')'"#",',sigfig+8
  write(str,formnum)1,2,3
  call write_string_master(ui,str)

  if(ie_global>is_global)then
   write(str,formhead)'  i','x1','dvol'
   call write_string_master(ui,str)
   call write_string_master(ui, '')
   j=js_global;k=ks_global
   do i = is_global, ie_global
    call write_my_grid_1d(ui,formval,i,x1(i),dvol,i,j,k)
   end do

  elseif(je_global>js_global)then
   write(str,formhead)'  j','x2','dvol'
   call write_string_master(ui,str)
   call write_string_master(ui, '')
   i=is_global;k=ks_global
   do j = js_global, je_global
    call write_my_grid_1d(ui,formval,j,x2(j),dvol,i,j,k)
   end do

  elseif(ke_global>ks_global)then
   write(str,formhead)'  k','x3','dvol'
   call write_string_master(ui,str)
   call write_string_master(ui, '')
   i=is_global;j=js_global
   do k = ks_global, ke_global
    call write_my_grid_1d(ui,formval,k,x3(k),dvol,i,j,k)
   end do

  end if

! 2D outputs
 case(2)
  write(formnum,'("(",a4,"i4,i5,3i",i2,")")')'"#",',sigfig+8
  write(str,formnum)1,2,3,4,5
  call write_string_master(ui,str)

  if(ke_global==ks_global)then! For 2D Cartesian, polar or axisymmetrical spherical
   write(str,formhead)'  i','j','x1','x2','dvol'
   call write_string_master(ui,str)
   call write_string_master(ui, '')
   k=ks_global

   ! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=js_global-1
    do i = is_global, ie_global
     call write_my_grid_2d(ui,formval,i,j,x1(i),xi2s,dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end if

   do j = js_global, je_global
    do i = is_global, ie_global
     call write_my_grid_2d(ui,formval,i,j,x1(i),x2(j),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do

   ! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=je_global+1
    do i = is_global, ie_global
     call write_my_grid_2d(ui,formval,i,j,x1(i),xi2e,dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end if

  elseif(je_global==js_global)then! mainly for 2D Cartesian or axisymmetrical cylindrical
   write(str,formhead)'  i','k','x1','x3','dvol'
   call write_string_master(ui,str)
   call write_string_master(ui, '')
   j=js_global
   do k = ks_global, ke_global, outres

    ! writing inner boundary for polar coordinates
    if(crdnt==1.or.crdnt==2) then
     call write_my_grid_2d(ui,formval,is_global-1,k,xi1(is_global-1),x3(k),dvol,is_global,j,k)
    end if

    do i = is_global, ie_global, outres
     call write_my_grid_2d(ui,formval,i,k,x1(i),x3(k),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do

  elseif(ie_global==is_global)then! For 2D Cartesian
  !CAUTION: Not designed for cylindrical or spherical yet
   write(str,formhead)'  j','k','x2','x3','dvol'
   call write_string_master(ui,str)
   call write_string_master(ui, '')
   i=is_global
   do k = ks_global, ke_global
    do j = js_global, je_global
     call write_my_grid_2d(ui,formval,j,k,x2(j),x3(k),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do
  end if

 case(3)
  write(formnum,'("(",a4,"i4,2i5,4i",i2,")")')'"#",',sigfig+8
  write(str,formnum)1,2,3,4,5,6,7
  call write_string_master(ui,str)

  write(str,formhead)'  i','j','k','x1','x2','x3','dvol'
  call write_string_master(ui,str)
  call write_string_master(ui, '')
  if(crdnt==2)then
   k=ks_global-1
   do j = je_global, je_global
    do i = is_global, ie_global
     call write_my_grid_3d(ui,formval,i,j,k,x1(i),x2(j),xi3(k),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do
  end if

  do k = ks_global, ke_global
   do j = je_global, je_global
    do i = is_global, ie_global
     call write_my_grid_3d(ui,formval,i,j,k,x1(i),x2(j),x3(k),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do
  end do

  if(crdnt==2)then
   k=ke_global+1
   do j = je_global, je_global
    do i = is_global, ie_global
     call write_my_grid_3d(ui,formval,i,j,k,x1(i),x2(j),xi3(k-1),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do
  end if

 case default
  stop 'Something wrong with dimension'
 end select

 call close_file(ui)

 if (myrank==0) print*,"Outputted: ",'gridfile.dat'

!othergrid--------------------------------------------------------------------

 if(write_other_slice)then

  call add_directory_to_filename('othergrid.dat',filename)
  call open_file_write_ascii(filename,ui)
  call write_string_master(ui, '')
  write(formhead,'("(",i1,"a5,",i1,"a",i2,")")')2,2+1,sigfig+8
  write(formval ,'("(",i1,"i5,",i1,"(1x,1PE",i2,".",i2,"e2))")')&
        2,2+1,sigfig+7,sigfig-1

  write(formnum,'("(",a4,"i4,i5,4i",i2,")")')'"#",',sigfig+8
  write(str,formnum)1,2,3,4,5
  call write_string_master(ui,str)

  write(str,formhead)'  i','j','x1','x2','dvol'
  call write_string_master(ui,str)
  call write_string_master(ui, '')

  if(crdnt==2)then
   k = ks_global

   do i = is_global, ie_global
    call write_my_grid_2d(ui,formval,i,-je_global-1,x1(i),-xi2(je_global),dvol,i,je_global,k)
   end do
   call write_string_master(ui, '')

   do j = je_global, js_global, -1
    do i = is_global, ie_global
     call write_my_grid_2d(ui,formval,i,-j,x1(i),-x2(j),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do

   do i = is_global, ie_global
    call write_my_grid_2d(ui,formval,i,0,x1(i),-xi2(js_global-1),dvol,i,js_global,k)
   end do
   call write_string_master(ui, '')

   k = (ks_global+ke_global-1)/2
   do i = is_global, ie_global
    call write_my_grid_2d(ui,formval,i,0,x1(i),xi2(js_global-1),dvol,i,js_global,k)
   end do
   call write_string_master(ui, '')

   do j = js_global, je_global
    do i = is_global, ie_global
     call write_my_grid_2d(ui,formval,i,j,x1(i),x2(j),dvol,i,j,k)
    end do
    call write_string_master(ui, '')
   end do

   do i = is_global, ie_global
    call write_my_grid_2d(ui,formval,i,je_global+1,x1(i),xi2(je_global),dvol,i,je_global,k)
   end do

  end if

  call close_file(ui)

  if (myrank==0) print*,"Outputted: ",'othergrid.dat'

 end if

 return
end subroutine write_grid_dat

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE WRITE_BIN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output binary full dump file

subroutine write_bin

 use settings
 use mpi_utils, only:myrank
 use io, only:open_file_write,close_file, write_var, write_dummy_recordmarker
 use grid,only:is,ie,js,je,ks,ke,gis,gie,gjs,gje,gks,gke,time,tn,frame_acc,frame_vel,frame_pos
 use physval
 use gravmod,only:grvphi,grvpsi,cgrav_old
 use sink_mod,only:nsink,sink

 implicit none

 character(len=70):: binfile
 integer :: un

!-----------------------------------------------------------------------------

 call set_file_name('bin',tn,time,binfile)
 call open_file_write(binfile, un)

 call write_dummy_recordmarker(un)
 call write_var(un, tn)
 call write_var(un, time)
 call write_dummy_recordmarker(un)

 call write_dummy_recordmarker(un)
 call write_var(un,  d, is, ie, js, je, ks, ke)
 call write_var(un, v1, is, ie, js, je, ks, ke)
 call write_var(un, v2, is, ie, js, je, ks, ke)
 call write_var(un, v3, is, ie, js, je, ks, ke)
 call write_var(un,  e, is, ie, js, je, ks, ke)
 call write_dummy_recordmarker(un)

 if(eostype>=1)then
  call write_dummy_recordmarker(un)
  call write_var(un, T, is, ie, js, je, ks, ke)
  call write_dummy_recordmarker(un)
 end if

 if(gravswitch>=2) then
   call write_dummy_recordmarker(un)
   call write_var(un, grvphi, gis, gie, gjs, gje, gks, gke, grav=.true.)
   call write_dummy_recordmarker(un)
 end if

 if(gravswitch==3) then
   call write_dummy_recordmarker(un)
   call write_var(un, grvpsi, gis, gie, gjs, gje, gks, gke, grav=.true.)
   call write_var(un, cgrav_old)
   call write_dummy_recordmarker(un)
 endif

 if(compswitch==1) then
   call write_dummy_recordmarker(un)
   call write_var(un, ubound(mudata,1))
   call write_var(un, mudata, 0, ubound(mudata,1), 0, 1)
   call write_dummy_recordmarker(un)
 elseif(compswitch>=2) then
   call write_dummy_recordmarker(un)
   call write_var(un, spn)
   call write_var(un, spc, 1, spn, is, ie, js, je, ks, ke)
   call write_var(un, species, 1, spn)
   call write_dummy_recordmarker(un)
 endif

 if(mag_on) then
  call write_dummy_recordmarker(un)
  call write_var(un, b1, is, ie, js, je, ks, ke)
  call write_var(un, b2, is, ie, js, je, ks, ke)
  call write_var(un, b3, is, ie, js, je, ks, ke)
  call write_var(un, phi, is, ie, js, je, ks, ke)
  call write_dummy_recordmarker(un)
 end if

 if(radswitch>0) then
  call write_dummy_recordmarker(un)
  call write_var(un, erad, is, ie, js, je, ks, ke)
  call write_dummy_recordmarker(un)
 end if

 if(include_sinks) then
  call write_dummy_recordmarker(un)
  call write_var(un, sink, 1, nsink)
  call write_dummy_recordmarker(un)
  if(frame>0) then
   call write_dummy_recordmarker(un)
   call write_var(un, frame_acc, 1, 3)
   call write_var(un, frame_vel, 1, 3)
   call write_var(un, frame_pos, 1, 3)
   call write_dummy_recordmarker(un)
  endif
 endif

 call close_file(un)

 if (myrank==0) print*,"Outputted: ",trim(binfile)

return
end subroutine write_bin

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE WRITE_PLT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output ascii file for plotting

subroutine write_plt
 use settings
 use grid
 use utils,only:gravpot1d
 use shockfind_mod,only:shockfind
 use mpi_utils, only:myrank
 use io, only:open_file_write_ascii, write_string_master, close_file
 character(len=50):: pltfile
 character(len=20):: header(50)='aaa',forma,forme,formi
 character(len=200):: str
 integer:: i,j,k,n,ui,columns

! Calculate gravitational potential here if spherical
 if(gravswitch==1) call gravpot1d
! Calculate shock position if required
 if(write_shock)call shockfind

! Set format
 write(forma,'("(a",i2,")")')sigfig+8 ! for strings
 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers
 write(formi,'("(i",i2,")")')sigfig+8 ! for integers

! Open file
 call set_file_name('plt',tn,time,pltfile)
 call open_file_write_ascii(pltfile,ui)

! Write time and time step
 write(str,'(a,i10,a,1PE12.4e2,a)')&
  '#tn =',tn,'  time= ',time/dt_unit_in_sec,dt_unit
 call write_string_master(ui,str)

! Decide what quantities to write
 call get_header(header,columns)
 do n = 1, columns
  write(str,formi) n+2*dim+1
  call write_string_master(ui,str,advance=.false.)
 end do
 call write_string_master(ui, '')
 do n = 1, columns
  write(str,forma) trim(adjustl(header(n)))
  call write_string_master(ui,str,advance=.false.)
 end do
 call write_string_master(ui, '')
 call write_string_master(ui, '')

! Write file
 select case (dim)
! 1D outputs
 case(1)

  if(ie_global>is_global)then
   j=js_global;k=ks_global
   do i = is_global, ie_global
    call write_val(ui,i,j,k,forme,header)
   end do
  elseif(je_global>js_global)then
   i=is_global;k=ks_global
   do j = js_global, je_global
    call write_val(ui,i,j,k,forme,header)
   end do
  elseif(ke_global>ks_global)then
   i=is_global;j=js_global
   do k = ks_global, ke_global
    call write_val(ui,i,j,k,forme,header)
   end do
  end if

! 2D outputs
 case(2)

  if(ke_global==ks_global)then! For 2D Cartesian, polar or axisymmetrical spherical
   k=ks_global
! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=js_global
    do i = is_global, ie_global
     call write_val(ui,i,j,k,forme,header)
    end do
    call write_string_master(ui, '')
   end if

   do j = js_global, je_global
    do i = is_global, ie_global
     call write_val(ui,i,j,k,forme,header)
    end do
    call write_string_master(ui, '')
   end do

! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=je_global
    do i = is_global, ie_global
     call write_val(ui,i,j,k,forme,header)
    end do
    call write_string_master(ui, '')
   end if

  elseif(je_global==js_global)then! mainly for 2D Cartesian or axisymmetrical cylindrical
   j=js_global
   do k = ks_global, ke_global, outres
! writing inner boundary for polar coordinates
    if(crdnt==1.or.crdnt==2)call write_val(ui,is_global,j,k,forme,header)
    do i = is_global, ie_global, outres
     call write_val(ui,i,j,k,forme,header)
    end do
    call write_string_master(ui, '')
   end do

  elseif(ie_global==is_global)then! For 2D Cartesian
!CAUTION: Not designed for cylindrical or spherical yet
   i=is_global
   do k = ks_global, ke_global
    do j = js_global, je_global
     call write_val(ui,i,j,k,forme,header)
    end do
    call write_string_master(ui, '')
   end do
  end if

 case(3)

  if(crdnt==2)then
   do j = je_global, je_global
    do i = is_global, ie_global
     call write_val(ui,i,j,ks_global,forme,header)
    end do
    call write_string_master(ui, '')
   end do
  end if

  do k = ks_global, ke_global
   do j = je_global, je_global
    do i = is_global, ie_global
     call write_val(ui,i,j,k,forme,header)
    end do
    call write_string_master(ui, '')
   end do
  end do

  if(crdnt==2)then
   do j = je_global, je_global
    do i = is_global, ie_global
     call write_val(ui,i,j,ke_global,forme,header)
    end do
    call write_string_master(ui, '')
   end do
  end if

 case default
  stop 'Something wrong with dimension'
 end select

 call close_file(ui)

 if (myrank==0) print*,"Outputted: ",trim(pltfile)

!othfile----------------------------------------------------------------------

 if(write_other_slice)then
  call write_vertical_slice(ks_global                          ,'oth')
  call write_vertical_slice(ks_global+(ke_global-ks_global+1)/4,'ver')
 end if

 return
end subroutine write_plt

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE WRITE_VERTICAL_SLICE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output ascii file for vertical slices in 3D

subroutine write_vertical_slice(slice,prefix)

 use physval
 use settings
 use grid,only:is_global,ie_global,js_global,je_global,ks_global,ke_global,&
               time,tn,dim
 use shockfind_mod,only:shockfind
 use io, only:open_file_write_ascii, write_string_master, close_file
 use mpi_utils, only:myrank

 integer,intent(in):: slice
 character(len=*),intent(in):: prefix
 character(len=50):: verfile
 character(len=20):: header(50)='aaa',forma,forme,formi
 character(len=200):: str
 integer:: i,j,k,n,ui,columns,kn

!-----------------------------------------------------------------------------

! Set format
 write(forma,'("(a",i2,")")')sigfig+8 ! for strings
 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers
 write(formi,'("(i",i2,")")')sigfig+8 ! for integers

! Open file
 call set_file_name(prefix,tn,time,verfile)
 call open_file_write_ascii(verfile,ui)

! Write time and time step
 write(str,'(a,i7,a,1PE12.4e2,a)')&
      '#tn =',tn,'  time= ',time/dt_unit_in_sec,dt_unit
 call write_string_master(ui,str)

! Decide what quantities to write
 call get_header(header,columns)
 do n = 1, columns
  write(str,formi) n+2*dim+1
  call write_string_master(ui,str,advance=.false.)
 end do
 call write_string_master(ui, '')
 do n = 1, columns
  write(str,forma) trim(adjustl(header(n)))
  call write_string_master(ui,str,advance=.false.)
 end do
 call write_string_master(ui, '')
 call write_string_master(ui, '')

 if(crdnt==2)then

  kn = ke_global-ks_global+1
  k = slice
  do i = is_global, ie_global
   call write_val(ui,i,je_global,k,forme,header)
  end do
  call write_string_master(ui, '')
  do j = je_global, js_global, -1
   do i = is_global, ie_global
    call write_val(ui,i,j,k,forme,header)
   end do
   call write_string_master(ui, '')
  end do
  do i = is_global, ie_global
   call write_val(ui,i,js_global,k,forme,header)
  end do
  call write_string_master(ui, '')

  k = mod(slice+kn/2,kn)
  do i = is_global, ie_global
   call write_val(ui,i,js_global,k,forme,header)
  end do
  call write_string_master(ui, '')
  do j = js_global, je_global
   do i = is_global, ie_global
    call write_val(ui,i,j,k,forme,header)
   end do
   call write_string_master(ui, '')
  end do
  do i = is_global, ie_global
   call write_val(ui,i,je_global,k,forme,header)
  end do

 end if

 call close_file(ui)

 if (myrank==0) print*,"Outputted: ",trim(verfile)

return
end subroutine write_vertical_slice


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE WRITE_EXTGRV
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To write extgrv file

subroutine write_extgrv
 use grid
 use gravmod, only:extgrv, mc
 use mpi_utils, only:myrank, allreduce_mpi
 use io, only:open_file_write,close_file, write_var, write_extgrv_array, write_dummy_recordmarker
 integer :: ui
 real(8) :: mcore
 character(len=70):: filename

 ! get mcore and ensure each task has the same value
 if (is==is_global) then
   mcore = mc(is-1)
 else
   mcore = -1.
 endif
 call allreduce_mpi('max', mcore)

 call add_directory_to_filename('extgrv.bin',filename)
 call open_file_write(filename, ui)

 call write_dummy_recordmarker(ui)
 call write_var(ui, mcore)
 call write_dummy_recordmarker(ui)

 call write_dummy_recordmarker(ui)
 call write_extgrv_array(ui, extgrv)
 call write_dummy_recordmarker(ui)

 call close_file(ui)

 if (myrank==0) print*,"Outputted: ",'extgrv.bin'

return
end subroutine write_extgrv

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                             SUBROUTINE WRITE_PTC
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Output ascii file for tracer particles
! WARNING: This routine is not tested yet

subroutine write_ptc

 use settings,only:crdnt,include_particles
 use grid,only:tn,time,xi1,xi2,xi3,is,ie,js,je,ks,ke
 use physval,only:d,e
 use particle_mod
 use gravmod,only:grvphi
 use mpi_utils,only:myrank

 character(len=50):: ptcfile
 integer::ui,i,j,k,n
 logical:: extflag
 real(8):: pr,pt

!-----------------------------------------------------------------------------

 if(.not.include_particles) return
 if(myrank/=0) return

 call set_file_name('ptc',tn,time,ptcfile)

 if(crdnt==1.and.je==js)then
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

 elseif(crdnt==2.and.ke==ks)then
  do n = 1, np
   extflag = .false.
   pr = sqrt(ptcx(1,n)**2+ptcx(2,n)**2)
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

 open(newunit=ui,file = ptcfile, status='replace')

 write(ui,'(a,i7,a,1PE12.4e2,2(a5,i9))')&
  '#tn =',tn,'  time= ',time,'np= ',np,'npl=',npl
 write(ui,'(a7,2a3,3a15)')'label','ej','ub','mass','x1','x3'

 do i = 1, np
  write(ui,'(i7,2i3,3(1PE15.7e2))')ptci(0:2,i),ptcx(0:2,i)
 end do

 close(ui)

 return
end subroutine write_ptc

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                             SUBROUTINE WRITE_BPT
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Output binary file for tracer particles

subroutine write_bpt

 use settings,only:include_particles
 use grid,only:tn,time
 use particle_mod
 use mpi_utils,only:myrank

 character(len=50):: bptfile
 integer::ui

!-----------------------------------------------------------------------------

 if(.not.include_particles) return
 if(myrank/=0) return

 call set_file_name('bpt',tn,time,bptfile)

 open(newunit=ui,file = bptfile, status='replace',form='unformatted')

 write(ui)np,npl
 write(ui)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)

 close(ui)

 return
end subroutine write_bpt


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE PROFILER_OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output wall time for each subroutine

subroutine profiler_output

 use settings
 use profiler_mod
 use omp_lib
 use mpi_utils, only:myrank

 integer:: ui,i,j
 character(len=30)::form1,forml

!-----------------------------------------------------------------------------

! Stop all active clocks
 do i = 0, n_wt
  if(clock_on(i))wtime(i) = wtime(i) + omp_get_wtime()
 end do

! Reduce clocks across MPI tasks
 call reduce_clocks_mpi

 if (myrank==0) then
  call get_maxlbl
  write(form1,'("(",i2,"X,4a12)")')maxlbl+1
  write(forml,'("(",i2,"a)")')maxlbl+1 + 4*12

  open(newunit=ui,file='walltime.dat',status='replace')
  write(ui,form1)'wall_time','frac_loop','frac_tot','imbalance'

  do i = 1, n_wt
    if(get_layer(i)==wttot)write(ui,forml)('-',j=1,maxlbl+1+4*12)
    call profiler_output1(ui,i)
  end do
  write(ui,forml)('-',j=1,maxlbl+1+4*12)
  call profiler_output1(ui,wttot)

  close(ui)
 endif

! Reactivate clocks
 do i = 0, n_wt
  if(clock_on(i))wtime(i) = wtime(i) - omp_get_wtime()
 end do

 return
end subroutine profiler_output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SCALING_OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output scaling.dat file for scaling tests

subroutine scaling_output

 use settings
 use profiler_mod
 use omp_lib

 integer:: i,un,n
 character(len=4):: is_scaling_test
 character(len=30)::form1

!-----------------------------------------------------------------------------

 call get_environment_variable("HORMONE_SCALING_TEST",is_scaling_test)

 if(is_scaling_test/='true')return

  call stop_clock(wttot)

  open(newunit=un,file='scaling.dat',status='old',position='append',iostat=i)

  if(i/=0)then
   open(newunit=un,file='scaling.dat',status='new')
   write(form1,'("(a8,",i2,"a14)")')n_wt+1
   write(un,form1)'threads',(trim(routine_name(n)),n=0,n_wt)
  end if

!$omp parallel
  n=omp_get_num_threads()
!$omp end parallel

  write(form1,'("(i8,",i2,"F14.6)")')n_wt+1
  write(un,form1)n,wtime(0:n_wt)

  close(un)

  stop
return
end subroutine scaling_output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE SET_FILE_NAME
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To create a file name with a given prefix

subroutine set_file_name(prefix,tn,time,filename)

 use settings,only:dt_unit_in_sec,dt_unit,outstyle

 implicit none

 character(len=*),intent(in):: prefix
 character(len=*),intent(out)::filename
 character(len=50)::filename1
 integer,intent(in):: tn
 real(8),intent(in):: time

!-----------------------------------------------------------------------------
 select case(outstyle)
 case(1)
  write(filename1,'(a,i11.11,a,a4)')&
   trim(prefix),nint(time/dt_unit_in_sec),trim(dt_unit),'.dat'
  if(time/dt_unit_in_sec>2147483647d0)then ! if integer is overflowed
   write(filename1,'(a,i9.9,"00",a,a4)')&
    trim(prefix),nint(time/dt_unit_in_sec*0.01d0),trim(dt_unit),'.dat'
  end if
 case(2)
  write(filename1,'(a,i8.8,a4)')trim(prefix),tn,'.dat'
 end select

 call add_directory_to_filename(filename1,filename)

 return
end subroutine set_file_name

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                     SUBROUTINE ADD_DIRECTORY_TO_FILENAME
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To add the directory to the filename string

subroutine add_directory_to_filename(filename,fullname)

 use settings,only:outdir

 character(len=*),intent(in):: filename
 character(len=*),intent(out):: fullname

!-----------------------------------------------------------------------------

 write(fullname,'(a,"/",a)') trim(outdir),trim(filename)

return
end subroutine add_directory_to_filename

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE GET_HEADER
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To determine what to output

subroutine get_header(header,columns)

 use settings
 use physval,only:species
 use grid,only:dim,is,ie,js,je,ks,ke

 implicit none

 character(len=*),intent(inout):: header(50)
 integer,intent(out):: columns
 integer::n

!-----------------------------------------------------------------------------
! Output density/internal energy/pressure by default
 columns = 0
 call add_column('d',columns,header)
 call add_column('e',columns,header)
 call add_column('p',columns,header)

! Decide which velocity components to output
 select case(dim)
 case(1)
  if(ie>is)then
   call add_column('v1',columns,header)
  elseif(je>js)then
   call add_column('v2',columns,header)
  elseif(ke>ks)then
   call add_column('v3',columns,header)
  end if
 case(2)
  if(ke==ks)then
   call add_column('v1',columns,header)
   call add_column('v2',columns,header)
   if(write_other_vel)call add_column('v3',columns,header)
  elseif(je==js)then
   call add_column('v1',columns,header)
   if(write_other_vel)call add_column('v2',columns,header)
   call add_column('v3',columns,header)
  elseif(ie==is)then
   if(write_other_vel)call add_column('v1',columns,header)
   call add_column('v2',columns,header)
   call add_column('v3',columns,header)
  end if
 case(3)
  call add_column('v1',columns,header)
  call add_column('v2',columns,header)
  call add_column('v3',columns,header)
 case default
  stop 'Something wrong with dimensions'
 end select

! Output magnetic field if mag_on
 if(mag_on)then
  call add_column('b1',columns,header)
  call add_column('b2',columns,header)
  call add_column('b3',columns,header)
 end if

! Output temperature
 if(write_temp)then
  if(radswitch>0)then
   call add_column('Tgas',columns,header)
   call add_column('Trad',columns,header)
  else
   call add_column('T',columns,header)
  end if
 end if

! Output gravitational potential if gravswitch>=1
 if(gravswitch>=1)then
  call add_column('phi',columns,header)
  if(include_extgrv)then
   call add_column('extphi',columns,header)
  end if
  if(include_sinks)then
   call add_column('totphi',columns,header)
  end if
 end if

! Output radiation energy density if radswitch>=1
 if(radswitch>0) &
  call add_column('erad',columns,header)

! Output mean molecular weight if compswitch>=1
 if(compswitch>=1)then
  call add_column('mu',columns,header)
  if(compswitch>=2)then
   do n = 1, spn
    call add_column(species(n),columns,header)
   end do
  end if
 end if

! Output shock position if write_shock=.true.
 if(write_shock)call add_column('shock',columns,header)

! Output mass coordinate if write_mc=.true.
 if(write_mc)call add_column('mc',columns,header)

return
end subroutine get_header

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE ADD_COLUMN
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To add a column of output variable into the header array

pure subroutine add_column(string,columns,header)

 implicit none

 character(len=*),intent(inout):: header(50)
 character(len=*),intent(in):: string
 integer,intent(inout):: columns

!-----------------------------------------------------------------------------

 columns = columns + 1
 header(columns) = string

return
end subroutine add_column

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                             SUBROUTINE WRITE_VAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To write a single row of values

subroutine write_val(ui,i,j,k,forme,header)
 use settings,only:spn
 use physval
 use gravmod,only:grvphi,extgrv,totphi,mc
 use pressure_mod,only:Trad
 use mpi_domain,only:is_my_domain
 use mpi_utils,only:barrier_mpi
 use io,only:write_string

 integer,intent(in):: ui,i,j,k
 character(len=*),intent(in):: forme, header(:)
 integer:: n, nn

!-----------------------------------------------------------------------------

 call barrier_mpi
 if (.not.is_my_domain(i,j,k)) return

 do n = 1, 50
  select case(header(n))
  case('d')!density
   call write_anyval(ui,forme,d(i,j,k))
  case('e')!internal energy
   call write_anyval(ui,forme,eint(i,j,k))
  case('p')!pressure
   call write_anyval(ui,forme,p(i,j,k))
  case('v1')!velocity 1
   call write_anyval(ui,forme,v1(i,j,k))
  case('v2')!velocity 2
   call write_anyval(ui,forme,v2(i,j,k))
  case('v3')!velocity 3
   call write_anyval(ui,forme,v3(i,j,k))
  case('b1')!magnetic field 1
   call write_anyval(ui,forme,b1(i,j,k))
  case('b2')!magnetic field 2
   call write_anyval(ui,forme,b2(i,j,k))
  case('b3')!magnetic field 3
   call write_anyval(ui,forme,b3(i,j,k))
  case('T')!temperature
   call write_anyval(ui,forme,T(i,j,k))
  case('Tgas')!gas temperature
   call write_anyval(ui,forme,T(i,j,k))
  case('Trad')!radiation temperature
   call write_anyval(ui,forme,Trad(erad(i,j,k)))
  case('phi')!gravitational potential
   call write_anyval(ui,forme,grvphi(i,j,k))
  case('extphi')!external gravitational potential
   call write_anyval(ui,forme,extgrv(i,j,k))
  case('totphi')!total gravitational potential
   call write_anyval(ui,forme,totphi(i,j,k))
  case('erad')!radiation energy
   call write_anyval(ui,forme,erad(i,j,k))
  case('mu')!mean molecular weight
   call write_anyval(ui,forme,1d0/imu(i,j,k))
  case('shock')!shock position
   call write_anyval(ui,forme,dble(shock(i,j,k)))
  case('mc')!mass coordinate
   call write_anyval(ui,forme,mc(i))
  case('aaa')!end of line
   call write_string(ui,'')
   exit
  case default!chemical composition
   do nn = 1, spn
    if(header(n)==species(nn))then
     call write_anyval(ui,forme,spc(nn,i,j,k))
    end if
   end do
  end select
 end do

return
end subroutine write_val

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE WRITE_ANYVAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To write one value

subroutine write_anyval(ui,forme,val,direct)
 use io, only:write_string
 integer,intent(in):: ui
 integer,intent(in),optional:: direct
 character(len=*),intent(in):: forme
 real(8),intent(in):: val
 character(len=30):: str
!-----------------------------------------------------------------------------

 if(abs(val)>=1d-99)then
  write(str,forme)val
 else
  write(str,forme)0d0
 end if
 if(present(direct))then
  write(ui,'(a)',advance='no') trim(str)
 else
  call write_string(ui,str,advance=.false.)
 end if

return
end subroutine write_anyval

subroutine write_my_grid_1d(ui,form,ii,x1,dvol,i,j,k)
  use mpi_domain, only:is_my_domain
  use mpi_utils, only:barrier_mpi
  use io, only:write_string
  integer, intent(in) :: ui,i,j,k,ii
  character(len=*), intent(in) :: form
  real(8), intent(in) :: x1
  real(8), allocatable, intent(in) :: dvol(:,:,:)
  character(len=100) :: str

  if (is_my_domain(i,j,k)) then
    write(str,form) ii,x1,dvol(i,j,k)
    call write_string(ui,str)
  end if

  call barrier_mpi

end subroutine write_my_grid_1d

subroutine write_my_grid_2d(ui,form,ii,jj,x1,x2,dvol,i,j,k)
  use mpi_domain, only:is_my_domain
  use mpi_utils, only:barrier_mpi
  use io, only:write_string
  integer, intent(in) :: ui,i,j,k,ii,jj
  character(len=*), intent(in) :: form
  real(8), intent(in) :: x1,x2
  real(8), allocatable, intent(in) :: dvol(:,:,:)
  character(len=200) :: str

  if (is_my_domain(i,j,k)) then
    write(str,form) ii,jj,x1,x2,dvol(i,j,k)
    call write_string(ui,str)
  end if

  call barrier_mpi

end subroutine write_my_grid_2d

subroutine write_my_grid_3d(ui,form,ii,jj,kk,x1,x2,x3,dvol,i,j,k)
  use mpi_domain, only:is_my_domain
  use mpi_utils, only:barrier_mpi
  use io, only:write_string
  integer, intent(in) :: ui,i,j,k,ii,jj,kk
  character(len=*), intent(in) :: form
  real(8), intent(in) :: x1,x2,x3
  real(8), allocatable, intent(in) :: dvol(:,:,:)
  character(len=300) :: str

  if (is_my_domain(i,j,k)) then
    write(str,form) ii,jj,kk,x1,x2,x3,dvol(i,j,k)
    call write_string(ui,str)
  end if

  call barrier_mpi

end subroutine write_my_grid_3d

end module output_mod
