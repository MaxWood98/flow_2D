!Flow solver general subroutines module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.6.4
!Updated: 14/08/23

!Module
module flow_general_mod
use flow_data_mod
contains 


!Allocate mesh subroutine =========================
subroutine allocate_mesh(mesh,ncell,nedge,nvtx)
implicit none 

!Variables - Import
integer(in) :: ncell,nedge,nvtx
type(mesh_data) :: mesh

!Allocate
mesh%ncell = ncell
mesh%nedge = nedge 
mesh%nvtx = nvtx 
allocate(mesh%edge_1(nedge))
allocate(mesh%edge_2(nedge))
allocate(mesh%cell_left(nedge))
allocate(mesh%cell_right(nedge))
allocate(mesh%zone_left(nedge))
allocate(mesh%zone_right(nedge))
allocate(mesh%cell_nedge(ncell))
allocate(mesh%cell_nedgei(ncell))
allocate(mesh%cell_nedgeiw(ncell))
allocate(mesh%cell_elentot(ncell))
allocate(mesh%cell_elenint(ncell))
allocate(mesh%cell_link(ncell))
allocate(mesh%cell_mz(ncell))
allocate(mesh%vtx_x(nvtx))
allocate(mesh%vtx_y(nvtx))
allocate(mesh%edge_dx(nedge))
allocate(mesh%edge_dy(nedge))
allocate(mesh%edgelen(nedge))
allocate(mesh%cell_vol(ncell))
allocate(mesh%deltaT(ncell))
allocate(mesh%deltaTD(ncell))
allocate(mesh%edge_nx(nedge))
allocate(mesh%edge_ny(nedge))
allocate(mesh%specrad(ncell))
allocate(mesh%psensor(ncell))
allocate(mesh%psensor_num(ncell))
allocate(mesh%psensor_dnum(ncell))
allocate(mesh%psensor_temp(ncell))
allocate(mesh%cellcenx(ncell))
allocate(mesh%cellceny(ncell))
allocate(mesh%edge_mx(nedge))
allocate(mesh%edge_my(nedge))
allocate(mesh%outflow_zone(nedge))
allocate(mesh%nedge_of_zone(nedge))
allocate(mesh%bc_state(nedge))
return
end subroutine allocate_mesh




!Deallocate mesh subroutine =========================
subroutine deallocate_mesh(mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh

!Deallocate 
deallocate(mesh%edge_1)
deallocate(mesh%edge_2)
deallocate(mesh%cell_left)
deallocate(mesh%cell_right)
deallocate(mesh%zone_left)
deallocate(mesh%zone_right)
deallocate(mesh%cell_nedge)
deallocate(mesh%cell_nedgei)
deallocate(mesh%cell_nedgeiw)
deallocate(mesh%cell_elentot)
deallocate(mesh%cell_elenint)
deallocate(mesh%cell_link)
deallocate(mesh%cell_mz)
deallocate(mesh%vtx_x)
deallocate(mesh%vtx_y)
deallocate(mesh%edge_dx)
deallocate(mesh%edge_dy)
deallocate(mesh%cell_vol)
deallocate(mesh%deltaT)
deallocate(mesh%deltaTD)
deallocate(mesh%edge_nx)
deallocate(mesh%edge_ny)
deallocate(mesh%specrad)
deallocate(mesh%psensor)
deallocate(mesh%psensor_num)
deallocate(mesh%psensor_dnum)
deallocate(mesh%psensor_temp)
deallocate(mesh%edgeonwcell)
deallocate(mesh%cellonwall)
deallocate(mesh%cellcenx)
deallocate(mesh%cellceny)
deallocate(mesh%edge_mx)
deallocate(mesh%edge_my)
deallocate(mesh%outflow_zone)
deallocate(mesh%nedge_of_zone)
deallocate(mesh%bc_state)
return
end subroutine deallocate_mesh




!Allocate flow variables subroutine =========================
subroutine allocate_flow_variables(flowvars,mesh,fcoefs,options)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars
type(coeficient_data) :: fcoefs
type(options_data) :: options

!Allocate 
allocate(flowvars%rho(mesh%ncell))
allocate(flowvars%u(mesh%ncell))
allocate(flowvars%v(mesh%ncell))
allocate(flowvars%p(mesh%ncell))
allocate(flowvars%e(mesh%ncell))
allocate(flowvars%mach(mesh%ncell))
allocate(flowvars%cp(mesh%ncell))
allocate(flowvars%resid(mesh%ncell))
allocate(flowvars%paverage(mesh%ncell))
allocate(flowvars%pgrad(mesh%ncell,2))
allocate(flowvars%u_outflow(mesh%nzone_outflow))
allocate(flowvars%v_outflow(mesh%nzone_outflow))
allocate(flowvars%rho_outflow(mesh%nzone_outflow))
allocate(flowvars%u_outflow_t(mesh%nzone_outflow))
allocate(flowvars%v_outflow_t(mesh%nzone_outflow))
allocate(flowvars%rho_outflow_t(mesh%nzone_outflow))
if (.NOT.allocated(fcoefs%cl_stencil)) then 
    allocate(fcoefs%cl_stencil(options%avstencil_size))
    fcoefs%cl_stencil(:) = 0.0d0 
end if 
if (.NOT.allocated(fcoefs%cd_stencil)) then 
    allocate(fcoefs%cd_stencil(options%avstencil_size))
    fcoefs%cd_stencil(:) = 0.0d0 
end if 
if (.NOT.allocated(fcoefs%cm_stencil)) then 
    allocate(fcoefs%cm_stencil(options%avstencil_size))
    fcoefs%cm_stencil(:) = 0.0d0 
end if 
if (.NOT.allocated(fcoefs%mflux_in_stencil)) then 
    allocate(fcoefs%mflux_in_stencil(options%avstencil_size))
    fcoefs%mflux_in_stencil(:) = 0.0d0 
end if 
return
end subroutine allocate_flow_variables




!Deallocate flow variables subroutine =========================
subroutine deallocate_flow_variables(flowvars,fcoefs)
implicit none 

!Variables - Import
type(flow_var_data) :: flowvars
type(coeficient_data) :: fcoefs

!Dellocate rho,u,v,p,e
deallocate(flowvars%rho)
deallocate(flowvars%u)
deallocate(flowvars%v)
deallocate(flowvars%p)
deallocate(flowvars%e)
deallocate(flowvars%mach)
deallocate(flowvars%cp)
deallocate(flowvars%resid)
deallocate(flowvars%pgrad)
deallocate(flowvars%u_outflow)
deallocate(flowvars%v_outflow)
deallocate(flowvars%rho_outflow)
deallocate(flowvars%u_outflow_t)
deallocate(flowvars%v_outflow_t)
deallocate(flowvars%rho_outflow_t)
deallocate(fcoefs%cl_stencil)
deallocate(fcoefs%cd_stencil)
deallocate(fcoefs%cm_stencil)
return
end subroutine deallocate_flow_variables




!Set RK coeficients subroutine ========================= 
subroutine set_RK_coeficients(fcoefs,options)
implicit none 

!Variables - Import
type(options_data) :: options
type(coeficient_data) :: fcoefs

!Set coeficients
fcoefs%RK_alfa(:) = 0
if (options%nRKstage == 3) then
    fcoefs%RK_alfa(1) = 2.0d0/3.0d0 
    fcoefs%RK_alfa(2) = 2.0d0/3.0d0 
    fcoefs%RK_alfa(3) = 1.0d0 
elseif (options%nRKstage == 4) then
    fcoefs%RK_alfa(1) = 1.0d0/4.0d0
    fcoefs%RK_alfa(2) = 1.0d0/3.0d0 
    fcoefs%RK_alfa(3) = 1.0d0/2.0d0
    fcoefs%RK_alfa(4) = 1.0d0
else 
    write(*,'(A)') '    invalid RK scheme selected'
    stop
end if 
return 
end subroutine set_RK_coeficients




!Parallel cell data transfer subroutines ========================= ** execute in OMP single only!!!!!!!! **
subroutine push_to_par_transfer_1D(par_transfer,dat2add,cell_map) !Push 1D
implicit none 

!Variables - Import
integer(in), dimension(:,:) :: cell_map
real(dp), dimension(:) :: par_transfer,dat2add

!Accumulate from thread to transfer array
par_transfer(cell_map(:,2)) = par_transfer(cell_map(:,2)) + dat2add(cell_map(:,1))
return 
end subroutine push_to_par_transfer_1D

subroutine pull_from_par_transfer_1D(par_transfer,dat2take,cell_map) !Pull 1D
implicit none 

!Variables - Import
integer(in), dimension(:,:) :: cell_map
real(dp), dimension(:) :: par_transfer,dat2take

!Accumulate from thread to transfer array
dat2take(cell_map(:,1)) = par_transfer(cell_map(:,2))
return 
end subroutine pull_from_par_transfer_1D

subroutine push_to_par_transfer_2D(par_transfer,dat2add,cell_map) !Push 2D
implicit none 

!Variables - Import
integer(in), dimension(:,:) :: cell_map
real(dp), dimension(:,:) :: par_transfer,dat2add

!Accumulate from thread to transfer array
par_transfer(cell_map(:,2),:) = par_transfer(cell_map(:,2),:) + dat2add(cell_map(:,1),:)
return 
end subroutine push_to_par_transfer_2D

subroutine pull_from_par_transfer_2D(par_transfer,dat2take,cell_map) !Pull 2D
implicit none 

!Variables - Import
integer(in), dimension(:,:) :: cell_map
real(dp), dimension(:,:) :: par_transfer,dat2take

!Accumulate from thread to transfer array
dat2take(cell_map(:,1),:) = par_transfer(cell_map(:,2),:)
return 
end subroutine pull_from_par_transfer_2D


end module flow_general_mod