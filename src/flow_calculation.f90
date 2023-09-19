!Flow solver calculation module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering
!Version: 0.7.1
!Updated: 19/09/23

!Boundary condition flags -------
! -> Wall = -1
! -> Far field = -2
! -> Mass flux inflow = -3 
! -> Mass flux outflow = -4
! -> Stagnation state inflow = -5
! -> Freestram inflow = -6
! -> Back pressure outflow = -7
!Boundary condition flags -------

!Module
module flow_calculation_mod
use flow_io_mod
! use flow_data_mod
use flow_general_mod
use flow_boundary_cond_mod
contains 

!Initialise flow variables subroutine =========================
subroutine initialise_flow(options,flowvars,flowvarsAV,mesh,fcoefs,dispt)
implicit none 

!Variables - Import
integer(in) :: dispt
type(options_data) :: options 
type(flow_var_data) :: flowvars,flowvarsAV
type(coeficient_data) :: fcoefs
type(mesh_data) :: mesh

!Variables - Local 
real(dp) :: pinf_AC,sosinf_AC

!Display
if ((options%csdisp) .AND. (dispt == 1)) then
    write(*,'(A)') '--> initialising flow'
end if

!Set angle of atack in radians 
options%aoa_rad = options%aoa*(pi/180.0d0)

!Set gas constants (with temperature scaling)
flowvars%gam = options%gamma
flowvars%R = options%R
! flowvarsAV%R = flowvars%R
! flowvarsAV%gam = flowvars%gam 

!Set scaled density and pressure from the input values 
flowvars%rhoinf = options%rhoinf/1.225d0
flowvars%tinf = (options%tinf/273.15)*(1.0d0/(options%gamma*options%R))

!Set initial flow variables at normalised scale
flowvars%cinf = sqrt(options%gamma*options%R*flowvars%tinf)
flowvars%pinf = (flowvars%cinf*flowvars%cinf*flowvars%rhoinf)/options%gamma
flowvars%machinf = options%machinf 
flowvars%velinf = options%machinf*flowvars%cinf

!Set mass flow properties
flowvars%mflux_pset = 0.0d0 
flowvars%mflux_pset0 = 0.0d0 

!Set backpressure static pressure
flowvars%backp_pset = options%backp_pratio*flowvars%pinf

!Set freestram velocity
flowvars%uinf = flowvars%velinf*cos(options%aoa_rad)
flowvars%vinf = flowvars%velinf*sin(options%aoa_rad)

!Set freestram total pressure, density and temperature
flowvars%p0inf = flowvars%pinf*(1.0d0 + 0.5d0*(flowvars%gam - 1.0d0)*flowvars%machinf*flowvars%machinf)**&
                               (flowvars%gam/(flowvars%gam - 1.0d0))
flowvars%rho0inf = flowvars%rhoinf*(1.0d0 + 0.5d0*(flowvars%gam - 1.0d0)*flowvars%machinf*flowvars%machinf)**&
                               (flowvars%gam/(flowvars%gam - 1.0d0))
flowvars%t0inf = flowvars%tinf*(1.0d0 + 0.5d0*(flowvars%gam - 1.0d0)*flowvars%machinf*flowvars%machinf)

!Set initial flow primative variables 
if (options%init_state == 'f') then !Freestream
    flowvars%rho(:) = flowvars%rhoinf
    flowvars%u(:) = flowvars%uinf
    flowvars%v(:) = flowvars%vinf
    flowvars%p(:) = flowvars%pinf
    ! flowvars%e(:) = (flowvars%pinf/((flowvars%gam - 1.0d0)*flowvars%rhoinf)) + 0.5d0*flowvars%velinf*flowvars%velinf
    flowvars%e(:) = (flowvars%pinf/(flowvars%gam - 1.0d0)) + 0.5d0*flowvars%rhoinf*flowvars%velinf*flowvars%velinf
    flowvars%mach(:) = flowvars%machinf
    flowvars%cp(:) = ((flowvars%p(:)/flowvars%pinf) - 1.0d0)/(0.5d0*flowvars%gam*flowvars%machinf*flowvars%machinf)
elseif (options%init_state == 's') then !Static 
    flowvars%rho(:) = flowvars%rhoinf
    flowvars%u(:) = 0.0d0 
    flowvars%v(:) = 0.0d0 
    flowvars%p(:) = flowvars%pinf
    ! flowvars%e(:) = (flowvars%pinf/((flowvars%gam - 1.0d0)*flowvars%rhoinf))
    flowvars%e(:) = (flowvars%pinf/(flowvars%gam - 1.0d0))
    flowvars%mach(:) = 0.0d0 
    flowvars%cp(:) = ((flowvars%p(:)/flowvars%pinf) - 1.0d0)/(0.5d0*flowvars%gam*flowvars%machinf*flowvars%machinf)
elseif (options%init_state == 'r') then !Restart 
    call read_restart_file('flow2d_restart.dat',flowvars,mesh,options)
else
    options%status = 3
    flowvars%rho_res = 0.0d0 
    flowvarsAV%rho_res = 0.0d0 
    call export_status(options,flowvars,flowvarsAV)
    write(*,'(A)') '--> ** invalid initialisation state selected **'
    stop 
end if 
flowvars%pgrad(:,:) = 0.0d0 

!Set freestream conservative variables
flowvars%Winf(1) = flowvars%rhoinf
flowvars%Winf(2) = flowvars%rhoinf*flowvars%uinf 
flowvars%Winf(3) = flowvars%rhoinf*flowvars%vinf 
flowvars%Winf(4) = (flowvars%pinf/(flowvars%gam - 1.0d0)*flowvars%rhoinf) + 0.5d0*flowvars%velinf*flowvars%velinf

!Initialise average force coeficient stencils
fcoefs%cl_stencil(:) = 0.0d0 
fcoefs%cd_stencil(:) = 0.0d0 
fcoefs%cm_stencil(:) = 0.0d0 

!Calculate flow property scale factors 
pinf_AC = options%rhoinf*options%R*options%tinf
sosinf_AC = sqrt(options%gamma*options%R*options%tinf)
! flowvars%pscale = pinf_AC/flowvars%pinf
! flowvars%velscale = sosinf_AC/flowvars%cinf
! flowvars%rhoscale = options%rhoinf/flowvars%rhoinf

!Display actual flow properties 
if ((options%csdisp) .AND. (dispt == 1)) then
    write(*,'(A)') '    Actual flow properties: '
    write(*,'(A,F12.5)') '     pressure (Pa): ',pinf_AC
    write(*,'(A,F12.5)') '     density (Kg/m^3): ',options%rhoinf
    write(*,'(A,F12.5)') '     speed of sound (m/s): ',sosinf_AC
    write(*,'(A,F12.5)') '     velocity (m/s): ',options%machinf*sosinf_AC
    write(*,'(A)') '    Scaled flow properties: '
    write(*,'(A,F12.5)') '     pressure : ',flowvars%pinf
    write(*,'(A,F12.5)') '     density : ',flowvars%rhoinf
    write(*,'(A,F12.5)') '     speed of sound : ',flowvars%cinf
    write(*,'(A,F12.5)') '     velocity : ',flowvars%velinf
end if 
return
end subroutine initialise_flow




!Set flow variable scale subroutine =========================
subroutine set_flowvars_scale(flowvars,options)
implicit none 

!Variables - Import
type(options_data) :: options 
type(flow_var_data) :: flowvars

!Variables - Local
real(dp) :: c_a,p_a,vel_scale,rho_scale,p_scale

!Actual pressure and speed of sound 
c_a = sqrt(options%gamma*options%R*options%tinf)
p_a = (options%rhoinf*c_a*c_a)/options%gamma

!Scale velocity 
vel_scale = c_a/flowvars%cinf
flowvars%u(:) = flowvars%u(:)*vel_scale
flowvars%v(:) = flowvars%v(:)*vel_scale

!Scale pressure 
p_scale = p_a/flowvars%pinf
flowvars%p(:) = flowvars%p(:)*p_scale

!Scale density 
rho_scale = options%rhoinf/flowvars%rhoinf
flowvars%rho(:) = flowvars%rho(:)*rho_scale
return 
end subroutine set_flowvars_scale




!Set boundary condition states subroutine =========================
subroutine set_fixedbc_states(mesh_par,flowvars_par,tr_idx,bc_type)
implicit none 

!Variables - Import
integer(in) :: tr_idx,bc_type
type(mesh_data), dimension(:), allocatable :: mesh_par
type(flow_var_data), dimension(:), allocatable :: flowvars_par

!Variables - Local 
integer(in) :: ee 
real(dp) :: Vnorm

!Check each edge 
do ee=1,mesh_par(tr_idx)%nedge
    if ((mesh_par(tr_idx)%cell_left(ee) == bc_type) .OR. (mesh_par(tr_idx)%cell_right(ee) == bc_type)) then 
        Vnorm = (flowvars_par(tr_idx)%uinf*mesh_par(tr_idx)%edge_nx(ee) + &
        flowvars_par(tr_idx)%vinf*mesh_par(tr_idx)%edge_ny(ee))
        if (Vnorm .GT. 0.0d0) then !set inflow
            mesh_par(tr_idx)%bc_state(ee) = 1
        elseif (Vnorm .LT. 0.0d0) then !set outflow
            mesh_par(tr_idx)%bc_state(ee) = -1
        elseif (Vnorm == 0.0d0) then !set symmetry 
            mesh_par(tr_idx)%bc_state(ee) = 0 
        end if 
    end if 
end do 
return 
end subroutine set_fixedbc_states




!Calculate conservative variables subroutines =========================
subroutine get_conservative_vars(W,flowvars) 
implicit none 

!Variables - Import
real(dp), dimension(:,:) :: W
type(flow_var_data) :: flowvars

!Find conservative variables in each cell 
W(:,1) = flowvars%rho(:)
W(:,2) = flowvars%rho(:)*flowvars%u(:)
W(:,3) = flowvars%rho(:)*flowvars%v(:)
! W(:,4) = flowvars%rho(:)*flowvars%e(:)
W(:,4) = flowvars%e(:)
return 
end subroutine get_conservative_vars




!Calculate primative variables subroutine =========================
subroutine get_primative_vars(W,flowvars,ncell,nanflag) 
implicit none 

!Variables - Import
integer(in) :: ncell,nanflag
real(dp), dimension(:,:) :: W
type(flow_var_data) :: flowvars

!Variables - Local 
integer(in) :: ii
real(dp) :: Vabs2,SoS

!Find primative variables in each cell 
nanflag = 0 
do ii=1,ncell
    flowvars%rho(ii) = W(ii,1)
    if (isnan(W(ii,1))) then
        nanflag = 1
        write(*,'(A)') '    nan rho'
        return
    end if
    flowvars%u(ii) = W(ii,2)/W(ii,1)
    flowvars%v(ii) = W(ii,3)/W(ii,1)
    ! flowvars%e(ii) = W(ii,4)/W(ii,1)
    flowvars%e(ii) = W(ii,4)
    Vabs2 = flowvars%u(ii)**2 + flowvars%v(ii)**2
    ! flowvars%p(ii) = (flowvars%e(ii) - 0.5d0*Vabs2)*(flowvars%gam - 1.0d0)*W(ii,1)
    flowvars%p(ii) = (flowvars%e(ii) - 0.5d0*W(ii,1)*Vabs2)*(flowvars%gam - 1.0d0)
    SoS = sqrt(flowvars%gam*(flowvars%p(ii)/flowvars%rho(ii)))
    flowvars%mach(ii) = sqrt(Vabs2)/SoS
    flowvars%cp(ii) = ((flowvars%p(ii)/flowvars%pinf) - 1.0d0)/(0.5d0*flowvars%gam*flowvars%machinf*flowvars%machinf)
    ! flowvars%cp(ii) = (flowvars%p(ii) - flowvars%pinf)/(0.5d0*flowvars%rhoinf*flowvars%velinf*flowvars%velinf)
end do 
return
end subroutine get_primative_vars




!Edge flux calculation subroutine =========================
subroutine edge_flux(Flxe,cpe,cl,cr,zl,zr,edge_idx,mesh,flowvars,options)
implicit none 

!Variables - Import
integer(in) :: cl,cr,zl,zr,edge_idx
real(dp) :: cpe
real(dp) :: Flxe(4)
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:)  :: mesh
type(options_data) :: options

!Variables - Local 
real(dp) :: Vnorm,pedge,ue,ve,mache,SoS
real(dp) :: Fleft(4),Fright(4),Gleft(4),Gright(4)
real(dp) :: Fedge(4),Gedge(4)

!Evaluate flux
if (cl .GT. 0) then    !Other cell
    call cell_flux(Fright,Gright,flowvars,cr,zr)
    call cell_flux(Fleft,Gleft,flowvars,cl,zl)
    pedge = 0.5d0*(flowvars(zl)%p(cl) + flowvars(zr)%p(cr))
elseif (cl == -1) then !Object surface wall --> enforce flux wall boundary condition
    ! call cell_flux(Fright,Gright,flowvars,cr,zr)
    if (options%wall_bc_type == 'l') then 
        ! call wall_flux_bc_Lextrp(pedge,Fleft,Gleft,flowvars,mesh,edge_idx)
    elseif (options%wall_bc_type == 'c') then 
        call wall_flux_bc_Cextrp(pedge,Fleft,Gleft,flowvars,cr,zr)
    else
        options%status = 1
        write(*,*) '** warning -> invalid wall boundary condition specified (must be either constant (c) or linear (l))'
        return 
    end if 
    Fright(:) = Fleft(:)
    Gright(:) = Gleft(:)
elseif (cl == -2) then !Far field --> enforce far field boundary conditions

    !Cell flux on internal cell 
    call cell_flux(Fright,Gright,flowvars,cr,zr)

    !Edge normal flow velocity (if Vnorm > 0 flow from cl -> cr thus inflow to cr)
    ue = flowvars(zr)%u(cr) 
    ve = flowvars(zr)%v(cr)
    Vnorm = (ue*mesh(zr)%edge_nx(edge_idx) + ve*mesh(zr)%edge_ny(edge_idx))

    !Edge normal mach number
    SoS = sqrt(flowvars(zr)%gam*(flowvars(zr)%p(cr)/flowvars(zr)%rho(cr)))
    mache = abs(Vnorm)/SoS

    !Force boundary condition inflow/outflow selection if enabled 
    if (options%force_fixed_ff_ifof_state) then 
        if (mesh(zr)%bc_state(edge_idx) == 1) then !force inflow 
            Vnorm = abs(Vnorm)
        elseif (mesh(zr)%bc_state(edge_idx) == -1) then !force outflow 
            Vnorm = -1.0d0*abs(Vnorm) 
        elseif (mesh(zr)%bc_state(edge_idx) == 0) then !force symmetry 
            Vnorm = 0.0d0 
        end if 
    end if 

    !Determine inflow or outflow on cell edge to cell right 
    if (Vnorm .GT. 0.0d0) then !Inflow to cell right
        if (mache .GE. 1.0d0) then !Supersonic
            call supersonic_inflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh)
        else !Subsonic
            if (options%ff_bc_type == 'r') then 
                call subsonic_inflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh)
            elseif (options%ff_bc_type == 'c') then 
                call subsonic_bc_char(Fleft,Gleft,cr,zr,flowvars,1_in)
            endif
        end if
    elseif (Vnorm .LT. 0.0d0) then !Outflow from cell right
        if (mache .GE. 1.0d0) then !Supersonic
            call supersonic_outflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh)
        else !Subsonic
            if (options%ff_bc_type == 'r') then 
                call subsonic_outflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh)
            elseif (options%ff_bc_type == 'c') then 
                call subsonic_bc_char(Fleft,Gleft,cr,zr,flowvars,-1_in)
            end if 
        end if
    else !flow is parallel to the boundary (use outflow conditions)
        if (mache .GE. 1.0d0) then !Supersonic
            call supersonic_outflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh)
        else !Subsonic
            if (options%ff_bc_type == 'r') then 
                call subsonic_outflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh)
            elseif (options%ff_bc_type == 'c') then 
                call subsonic_bc_char(Fleft,Gleft,cr,zr,flowvars,-1_in)
            end if 
        end if
    end if
    pedge = 0.5d0*(flowvars(zr)%pinf + flowvars(zr)%p(cr))
elseif (cl == -3) then !Inflow --> mass flow inflow boundary condition enforced elsewhere
    Fright(:) = 0.0d0 
    Gright(:) = 0.0d0 
    Fleft(:) = 0.0d0 
    Gleft(:) = 0.0d0 
    pedge = flowvars(zr)%p(cr)
elseif (cl == -4) then !Outflow --> mass flow outflow boundary condition enforced elsewhere
    Fright(:) = 0.0d0 
    Gright(:) = 0.0d0 
    Fleft(:) = 0.0d0 
    Gleft(:) = 0.0d0 
    pedge = flowvars(zr)%p(cr)
elseif (cl == -5) then !Inflow --> enforce stagnation state inflow boundary condition
    call cell_flux(Fright,Gright,flowvars,cr,zr)
    call stagnation_inflow_bc(Fleft,Gleft,flowvars,mesh,edge_idx,cr,zr)
    pedge = flowvars(zr)%p(cr)
elseif (cl == -6) then !Inflow --> enforce 'freestream' state inflow boundary condition
    call cell_flux(Fright,Gright,flowvars,cr,zr)
    call freestream_inflow_c(Fleft,Gleft,cr,zr,flowvars,mesh,options,edge_idx)
    pedge = flowvars(zr)%p(cr)
elseif (cl == -7) then !Outflow --> enforce back pressure outflow boundary condition
    call cell_flux(Fright,Gright,flowvars,cr,zr)
    if (options%ff_bc_type == 'r') then 
        call backpressure_outflow_bc_ri(Fleft,Gleft,flowvars,mesh,edge_idx,cr,zr)
    elseif (options%ff_bc_type == 'c') then 
        call backpressure_outflow_bc_c(Fleft,Gleft,flowvars,mesh,options,edge_idx,cr,zr)
    end if 
    pedge = flowvars(zr)%p(cr)
end if 

!Edge flux value ([dy,-dx] normal)
Fedge(:) = 0.5d0*(Fleft(:) + Fright(:))*mesh(zr)%edge_dy(edge_idx)
Gedge(:) = 0.5d0*(Gleft(:) + Gright(:))*mesh(zr)%edge_dx(edge_idx)
Flxe(:) = Fedge(:) - Gedge(:)

!Pressure coefficient on edge
cpe = ((pedge/flowvars(zr)%pinf) - 1.0d0)/(0.5d0*flowvars(zr)%gam*flowvars(zr)%machinf*flowvars(zr)%machinf)
! cpe = (pedge - flowvars%pinf)/(0.5d0*flowvars%rhoinf*flowvars%velinf*flowvars%velinf)
return 
end subroutine edge_flux




!Cell spectral radius calculation subroutine =========================
subroutine cell_spectral_radius(mesh,flowvars,zone)
implicit none 

!Variables - Import
integer(in) :: zone
type(mesh_data), dimension(:) :: mesh
type(flow_var_data), dimension(:) :: flowvars

!Variables - Local 
integer(in) :: ii 
real(dp) :: lamedge

!Calculate spectral radius on each edge and accumulate to cells
do ii=1,mesh(zone)%nedge
    
    !Calculate edge flow spectral radius
    call edge_spectral_radius(lamedge,mesh(zone),flowvars,ii)

    !Add edge contribution to cells 
    if (mesh(zone)%cell_mz(mesh(zone)%cell_right(ii)) == 0) then 
        mesh(zone)%specrad(mesh(zone)%cell_right(ii)) = mesh(zone)%specrad(mesh(zone)%cell_right(ii)) + lamedge
    else    
        !$OMP atomic update 
        mesh(zone)%specrad(mesh(zone)%cell_right(ii)) = mesh(zone)%specrad(mesh(zone)%cell_right(ii)) + lamedge
    end if 
    if (mesh(zone)%cell_left(ii) .GT. 0) then
        if (mesh(mesh(zone)%zone_left(ii))%cell_mz(mesh(zone)%cell_left(ii)) == 0) then 
            mesh(mesh(zone)%zone_left(ii))%specrad(mesh(zone)%cell_left(ii)) = &
            mesh(mesh(zone)%zone_left(ii))%specrad(mesh(zone)%cell_left(ii)) + lamedge
        else
            !$OMP atomic update
            mesh(mesh(zone)%zone_left(ii))%specrad(mesh(zone)%cell_left(ii)) = &
            mesh(mesh(zone)%zone_left(ii))%specrad(mesh(zone)%cell_left(ii)) + lamedge
        end if 
    end if 
end do
return
end subroutine cell_spectral_radius




!Edge spectral radius calculation subroutine =========================
subroutine edge_spectral_radius(lamedge,mesh,flowvars,edg_idx)
implicit none 

!Variables - Import 
integer(in) :: edg_idx
real(dp) :: lamedge
type(mesh_data) :: mesh
type(flow_var_data), dimension(:) :: flowvars

!Variables - Local 
integer(in) :: cl,cr,zl,zr
real(dp) :: pl,pr,rhor,rhol
real(dp) :: sosl,sosr,ul,ur,vl,vr,ue,ve,sose

!Left and right cells and zones 
cl = mesh%cell_left(edg_idx)
cr = mesh%cell_right(edg_idx)
zl = mesh%zone_left(edg_idx)
zr = mesh%zone_right(edg_idx)

!Find right cell primatives 
pr = flowvars(zr)%p(cr)
rhor = flowvars(zr)%rho(cr)
ur = flowvars(zr)%u(cr)
vr = flowvars(zr)%v(cr)
sosr = sqrt(flowvars(zr)%gam*(pr/rhor))

!Consturct edge values ue,ve,ce
if (cl .GT. 0) then !Normal cell

    !Take values from left cell
    pl = flowvars(zl)%p(cl)
    rhol = flowvars(zl)%rho(cl)
    ul = flowvars(zl)%u(cl)
    vl = flowvars(zl)%v(cl)

    !Calculate speeds of sound in each cell
    sosl = sqrt(flowvars(zl)%gam*(pl/rhol))

    !Edge flow parameters
    ue = 0.5d0*(ul + ur)
    ve = 0.5d0*(vl + vr)
    sose = 0.5d0*(sosl + sosr)
else !Other cell
    ue = ur
    ve = vr
    sose = sosr
end if

!Calculate edge flow spectral radius
lamedge = (abs(mesh%edge_nx(edg_idx)*ue + mesh%edge_ny(edg_idx)*ve) + sose)*mesh%edgelen(edg_idx)
return 
end subroutine edge_spectral_radius




!Edge pressure sensor calculation subroutine =========================
subroutine edge_pressure_sensor(psensor_num,psensor_dnum,cl,cr,zl,zr,eidx,flowvars,mesh)
implicit none 

!Variables - Import
integer(in) :: cl,cr,zl,zr,eidx
real(dp) :: psensor_num,psensor_dnum
type(mesh_data), dimension(:) :: mesh
type(flow_var_data), dimension(:) :: flowvars

!Variables - Local
real(dp) :: pl,pr

!Calculate pressure sensor value
if (cl .GT. 0) then !Internal 
    pl = abs(flowvars(zl)%p(cl))
    pr = abs(flowvars(zr)%p(cr))
    psensor_num = (pr - pl)*mesh(zr)%edgelen(eidx)
    psensor_dnum = (pr + pl)*0.5d0*(mesh(zl)%cell_elenint(cl) +  mesh(zr)%cell_elenint(cr))
else !Other cells
    psensor_num = 0.0d0 
    psensor_dnum = 0.0d0 
end if
return
end subroutine edge_pressure_sensor


    

!Pressure averageing subroutine =========================
subroutine accumulate_average_pressure(flowvars,mesh,zone)
implicit none 

!Variables - Import 
integer(in) :: zone
type(mesh_data), dimension(:) :: mesh
type(flow_var_data), dimension(:) :: flowvars

!Variables - Local 
integer(in) :: ii,cl,cr,zl,zr
real(dp) :: pl,pr

!Smooth 
do ii=1,mesh(zone)%nedge

    !Left and right cells and zones for this edge
    cl = mesh(zone)%cell_left(ii)
    cr = mesh(zone)%cell_right(ii)
    zl = mesh(zone)%zone_left(ii)
    zr = mesh(zone)%zone_right(ii)

    !Cases
    if (cl .GT. 0) then !Internal cell

        !Edge pressures
        pl = mesh(zl)%psensor_temp(cl) 
        pr = mesh(zone)%psensor_temp(cr) 

        !Accumulate
        if (mesh(zr)%cell_mz(cr) == 0) then 
            flowvars(zr)%paverage(cr) = flowvars(zr)%paverage(cr) + pl 
        else    
            !$OMP atomic  
            flowvars(zr)%paverage(cr) = flowvars(zr)%paverage(cr) + pl 
        end if 
        if (cl .GT. 0) then
            if (mesh(zl)%cell_mz(cl) == 0) then 
                flowvars(zl)%paverage(cl) = flowvars(zl)%paverage(cl) + pr 
            else
                !$OMP atomic 
                flowvars(zl)%paverage(cl) = flowvars(zl)%paverage(cl) + pr 
            end if 
        end if
    end if
end do 
return 
end subroutine accumulate_average_pressure




!Pressure sensor smoothing subroutine =========================
subroutine smooth_psensor(mesh,zone) 
implicit none 

!Variables - Import 
integer(in) :: zone
type(mesh_data), dimension(:) :: mesh

!Variables - Local 
integer(in) :: ii,cl,cr,zl,zr
real(dp) :: pl,pr

!Smooth 
do ii=1,mesh(zone)%nedge

    !Left and right cells and zones for this edge
    cl = mesh(zone)%cell_left(ii)
    cr = mesh(zone)%cell_right(ii)
    zl = mesh(zone)%zone_left(ii)
    zr = mesh(zone)%zone_right(ii)

    !Cases
    if (cl .GT. 0) then !Internal cell
        pl = mesh(zl)%psensor(cl) 
        pr = mesh(zone)%psensor(cr) 
    else !Other cell
        pl = mesh(zone)%psensor(cr) 
        pr = mesh(zone)%psensor(cr)
    end if

    !Accumulate
    if (mesh(zr)%cell_mz(cr) == 0) then 
        mesh(zr)%psensor_temp(cr) = mesh(zr)%psensor_temp(cr) + pl 
    else    
        !$OMP atomic  
        mesh(zr)%psensor_temp(cr) = mesh(zr)%psensor_temp(cr) + pl 
    end if 
    if (cl .GT. 0) then
        if (mesh(zl)%cell_mz(cl) == 0) then 
            mesh(zl)%psensor_temp(cl) = mesh(zl)%psensor_temp(cl) + pr 
        else
            !$OMP atomic 
            mesh(zl)%psensor_temp(cl) = mesh(zl)%psensor_temp(cl) + pr 
        end if 
    end if
end do 
return 
end subroutine smooth_psensor




!Evaluate wall adjacent pressure gradient subroutine =========================
subroutine evaluate_wadj_pgrad(flowvars,mesh)
implicit none 

!Variables - Import 
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Variables - Local 
integer(in) :: ii,edg_idx,cl,cr 
real(dp) :: pl,pr,pe,gfx,gfy 

!Reset 
flowvars%pgrad(mesh%cellonwall,:) = 0.0d0 

!Evaluate pressure gradients 
do ii=1,mesh%npwalledge

    !Target edge 
    edg_idx = mesh%edgeonwcell(ii)
    cl = mesh%cell_left(edg_idx)
    cr = mesh%cell_right(edg_idx)

    !Get edge pressures
    pr = flowvars%p(cr)
    if (cl .GT. 0) then !other cell 
        pl = flowvars%p(cl)  
    else !boundary condition
        pl = pr 
    end if
    pe = 0.5d0*(pl + pr)

    !Evaluate pressure flux 
    gfx = mesh%edge_dy(edg_idx)*pe
    gfy = -mesh%edge_dx(edg_idx)*pe
    flowvars%pgrad(cr,1) = flowvars%pgrad(cr,1) - gfx
    flowvars%pgrad(cr,2) = flowvars%pgrad(cr,2) - gfy
    if (cl .GT. 0) then 
        flowvars%pgrad(cl,1) = flowvars%pgrad(cl,1) + gfx
        flowvars%pgrad(cl,2) = flowvars%pgrad(cl,2) + gfy
    end if 
end do 
flowvars%pgrad(mesh%cellonwall,1) = flowvars%pgrad(mesh%cellonwall,1)/mesh%cell_vol(mesh%cellonwall)
flowvars%pgrad(mesh%cellonwall,2) = flowvars%pgrad(mesh%cellonwall,2)/mesh%cell_vol(mesh%cellonwall)
return 
end subroutine evaluate_wadj_pgrad




!Edge flow laplacian calculation subroutine =========================
subroutine edge_laplacian(Lp_edge,cl,cr,zl,zr,fres_par)
implicit none 

!Variables - Import 
integer(in) :: cl,cr,zl,zr
real(dp) :: Lp_edge(4)
type(parallel_residual), dimension(:) :: fres_par

!Variables - Local 
real(dp) :: Wl(4),Wr(4)

!Cases
if (cl .GT. 0) then !Internal cell

    !Conservative variables in each adjacent cell
    Wl(:) = fres_par(zl)%W(cl,:)
    Wr(:) = fres_par(zr)%W(cr,:)

    !Edge laplacian 
    Lp_edge(:) = Wr(:) - Wl(:)
else !Other cell
    Lp_edge(:) = 0.0d0 
end if
return 
end subroutine edge_laplacian

    


!Edge artificial dissipation coeficient calculation subroutine =========================    
subroutine edge_ADcoeff(dk2,dk4,dpsilam,cl,cr,zl,zr,mesh,options,flowvars,edg_idx,zone)
implicit none 

!Variables - Import 
integer(in) :: cl,cr,zl,zr,edg_idx,zone
real(dp) :: dk2,dk4,dpsilam
type(mesh_data), dimension(:) :: mesh
type(options_data) :: options 
type(flow_var_data), dimension(:) :: flowvars

!Variables - Local 
real(dp) :: lamedge,psi0,psi1,psi01,S2,S4,si

!Cell type cases
if (cl .GT. 0) then !Internal cell

    !Calculate edge flow spectral radius
    call edge_spectral_radius(lamedge,mesh(zone),flowvars,edg_idx)

    !4th order psi coefficient 
    psi0 = sqrt(mesh(zr)%specrad(cr)/(4.0d0*lamedge))!**(2.0d0/3.0d0)
    psi1 = sqrt(mesh(zl)%specrad(cl)/(4.0d0*lamedge))!**(2.0d0/3.0d0)
    psi01 = (4.0d0*(psi0*psi1))/(psi0 + psi1)

    !Cell edge quantity scaling values 
    S2 = 2.0d0*(real(mesh(zl)%cell_nedgei(cl),dp) + real(mesh(zr)%cell_nedgei(cr),dp))/&
               (real(mesh(zl)%cell_nedgei(cl),dp)*real(mesh(zr)%cell_nedgei(cr),dp))
    S4 = 0.25d0*(S2*S2)
    ! S2 = mesh(zr)%edgelen(edg_idx)/(0.5d0*(mesh(zr)%cell_elenint(cr) + mesh(zl)%cell_elenint(cl)))
    ! S4 = 0.25d0*(S2*S2)

    !Pressure sensor and base second order dissipation
    si = 0.5d0*abs(mesh(zl)%psensor(cl) + mesh(zr)%psensor(cr)) !+ options%sp

    !Construct dissipation coeficients for this edge
    dk2 = options%k2*si*S2  
    dk4 = max(0.0d0,(options%k4 - options%c2*dK2))*S4
    dk2 = dk2 + options%sp !dk4*options%sp
    dpsilam = lamedge*psi01
else !Other cell
    dk2 = 0.0d0 
    dk4 = 0.0d0 
    dpsilam = 0.0d0 
end if
return 
end subroutine edge_ADcoeff




!Edge artificial dissipation calculation subroutine =========================
subroutine edge_dissipation(Dedge,dK2,dK4,dpsilam,cl,cr,zl,zr,fres_par)
implicit none 

!Variables - Import 
integer(in) :: cl,cr,zl,zr
real(dp) :: dK2,dK4,dpsilam,Dedge(4)
type(parallel_residual), dimension(:) :: fres_par

!Cases
if (cl .GT. 0) then !Internal cell
    Dedge(:) = (dK2*(fres_par(zr)%W(cr,:) - fres_par(zl)%W(cl,:)) - dK4*&
                (fres_par(zr)%cell_Lp(cr,:) - fres_par(zl)%cell_Lp(cl,:)))*dpsilam
else
    Dedge(:) = 0.0d0 
end if
return 
end subroutine edge_dissipation




!Force coeficient averaging subroutine =========================
subroutine average_fcoefficients(fcoefs,options,iter)
implicit none 

!Variables - Import
integer(in) :: iter
type(options_data) :: options
type(coeficient_data) :: fcoefs 

!Variables - Local 
integer(in) :: stins
real(dp) :: cl_sd,cd_sd,cm_sd

!Stencil index
stins = mod(iter-1,options%avstencil_size) + 1

!Populate stencils
fcoefs%cl_stencil(stins) = fcoefs%cl
fcoefs%cd_stencil(stins) = fcoefs%cd
fcoefs%cm_stencil(stins) = fcoefs%cm

!Average components
if (iter .LT. options%avstencil_size) then
    fcoefs%cl_av = sum(fcoefs%cl_stencil(1:iter))/real(iter,dp)
    fcoefs%cd_av = sum(fcoefs%cd_stencil(1:iter))/real(iter,dp)
    fcoefs%cm_av = sum(fcoefs%cm_stencil(1:iter))/real(iter,dp)
else
    fcoefs%cl_av = sum(fcoefs%cl_stencil(:))/real(options%avstencil_size,dp)
    fcoefs%cd_av = sum(fcoefs%cd_stencil(:))/real(options%avstencil_size,dp)
    fcoefs%cm_av = sum(fcoefs%cm_stencil(:))/real(options%avstencil_size,dp)
end if

!Standard deviation of each coefficient and the average of all  
if (iter .LT. options%avstencil_size) then
    cl_sd = sqrt((1.0d0/real(iter,dp))*sum((fcoefs%cl_stencil(1:iter) - fcoefs%cl_av)**2))
    cd_sd = sqrt((1.0d0/real(iter,dp))*sum((fcoefs%cd_stencil(1:iter) - fcoefs%cd_av)**2))
    cm_sd = sqrt((1.0d0/real(iter,dp))*sum((fcoefs%cm_stencil(1:iter) - fcoefs%cm_av)**2))
else
    cl_sd = sqrt((1.0d0/real(options%avstencil_size,dp))*sum((fcoefs%cl_stencil(:) - fcoefs%cl_av)**2))
    cd_sd = sqrt((1.0d0/real(options%avstencil_size,dp))*sum((fcoefs%cd_stencil(:) - fcoefs%cd_av)**2))
    cm_sd = sqrt((1.0d0/real(options%avstencil_size,dp))*sum((fcoefs%cm_stencil(:) - fcoefs%cm_av)**2))
end if
fcoefs%force_sd = (cl_sd + cd_sd + cm_sd)/3.0d0 
return 
end subroutine average_fcoefficients




!Long term residual trend calculation subroutine =========================
subroutine long_residual_sd(res_sd,rhores_hist,rhores_hist_av,longres_sum,iter)
implicit none 

!Variables - Import
integer(in) :: iter
real(dp) :: longres_sum,res_sd
real(dp) , dimension(:) :: rhores_hist,rhores_hist_av

!Variables - Local 
integer(in) :: itlow,ithigh,avlow,avhigh,avdelta,avrange
real(dp) :: avres

!Averaging range
avrange = 1000

!Set bound iterations for the average 
itlow = iter - avrange
ithigh = iter 
avlow = itlow
if (avlow .LT. 1) then 
    avlow = 1
end if
avhigh = ithigh
avdelta = avhigh - avlow
if (avdelta == 0) then 
    avdelta = 1
end if

!Accumulate to long sum
longres_sum = longres_sum + rhores_hist(iter)
if (itlow .GE. 1) then 
    longres_sum = longres_sum - rhores_hist(itlow)
end if 

!Set current average 
avres = longres_sum/real(avdelta,dp)
rhores_hist_av(iter) = avres

!Find standard deviation of the average 
res_sd = sqrt((1.0d0/real(avdelta,dp))*sum((rhores_hist_av(avlow:avhigh) - avres)**2))
return 
end subroutine long_residual_sd




!Evaluate flowfield properties subroutine =========================
subroutine evaluate_flowfield_properties(flowvars,mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Evaluate vorticity
call get_flowfield_vorticity(flowvars,mesh)

!Evaluate velocity 
call get_flowfield_velocity(flowvars,mesh)

!Evaluate temperature
call get_flowfield_temp(flowvars,mesh)
return 
end subroutine evaluate_flowfield_properties




!Evaluate flowfield vorticity subroutine =========================
subroutine get_flowfield_vorticity(flowvars,mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Variables - Local
integer(in) :: ii,cl,cr
real(dp) :: uedge,vedge,gfx,gfy
real(dp) :: ugrad_y(mesh%ncell),vgrad_x(mesh%ncell)

!Allocate vorticity array
allocate(flowvars%vorticity(mesh%ncell))

!Initialise
ugrad_y(:) = 0.0d0 
vgrad_x(:) = 0.0d0 
flowvars%vorticity(:) = 0.0d0 

!Evaluate velocity gradients 
do ii=1,mesh%nedge 

    !Left and right cells for this edge
    cl = mesh%cell_left(ii)
    cr = mesh%cell_right(ii)

    !Edge velocity 
    if (cl .GT. 0) then 
        uedge = 0.5d0*(flowvars%u(cr) + flowvars%u(cl))
        vedge = 0.5d0*(flowvars%v(cr) + flowvars%v(cl))
    else
        uedge = flowvars%u(cr)
        vedge = flowvars%v(cr)
    end if

    !Gradients 
    gfx = vedge*mesh%edge_dy(ii)
    gfy = -uedge*mesh%edge_dx(ii)

    !Accumulate gradients 
    ugrad_y(cr) = ugrad_y(cr) - gfy
    vgrad_x(cr) = vgrad_x(cr) - gfx
    if (cl .GT. 0) then 
        ugrad_y(cl) = ugrad_y(cl) + gfy
        vgrad_x(cl) = vgrad_x(cl) + gfx
    end if
end do 

!Find actual values (comment out for area normalised vorticity)
do ii=1,mesh%ncell
    ugrad_y(ii) = ugrad_y(ii)/mesh%cell_vol(ii)
    vgrad_x(ii) = vgrad_x(ii)/mesh%cell_vol(ii)
end do 

!Evaluate vorticity 
flowvars%vorticity(:) = -(ugrad_y(:) - vgrad_x(:))
return 
end subroutine get_flowfield_vorticity




!Evaluate flowfield velocity subroutine =========================
subroutine get_flowfield_velocity(flowvars,mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Variables - Local
integer(in) :: ii

!Allocate velocity
allocate(flowvars%vel(mesh%ncell))

!Evaluate velocity
do ii=1,mesh%ncell
    flowvars%vel(ii) = sqrt(flowvars%u(ii)**2 + flowvars%v(ii)**2)
end do 
return 
end subroutine get_flowfield_velocity




!Evaluate flowfield temperature subroutine =========================
subroutine get_flowfield_temp(flowvars,mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Variables - Local
integer(in) :: ii

!Allocate temperature
allocate(flowvars%temp(mesh%ncell))

!Evaluate temperature
do ii=1,mesh%ncell
    flowvars%temp(ii) = flowvars%p(ii)/(flowvars%R*flowvars%rho(ii))
end do
return 
end subroutine get_flowfield_temp


end module flow_calculation_mod 