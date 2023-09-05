!Flow solver solution itteration module (parallel capable)
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering
!Version: 0.2.1
!Updated: 21/08/23

!Module
module flow_iterator_mod
use omp_lib
use flow_io_mod
! use flow_data_mod
use flow_mesh_mod
use flow_general_mod
use flow_calculation_mod
contains 


!Flow timestepping subroutine =========================
subroutine flow_solve(mesh,options,flowvars,flowvarsAV,fcoefs)
use ieee_arithmetic, only: ieee_value,IEEE_QUIET_NAN
implicit none 

!Variables - Import
type(options_data) :: options
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars,flowvarsAV
type(coeficient_data) :: fcoefs

!Variables - Local 
integer(in) :: cc,tt,fitt,nanflag,tr_idx,resconv
real(dp) :: rho_res,res_sd,longres_sum,mfluxin
real(dp) :: rhores_hist(options%ittlim),rhores_hist_av(options%ittlim),massfluxsten(options%massflux_niterav)
type(mesh_data), dimension(:), allocatable :: mesh_par
type(parallel_residual), dimension(:), allocatable :: fres_par
type(coeficient_data), dimension(:), allocatable :: fcoefs_par
type(flow_var_data), dimension(:), allocatable :: flowvars_par,flowvarsAV_par

!Initialise error state
nanflag = 0 

!Set number of threads
call omp_set_num_threads(options%num_threads)

!Segment mesh into parallel regions 
call flow_par_segment_mesh(mesh_par,mesh,options,nanflag)

!Initialise the global flow domain 
call allocate_flow_variables(flowvars,mesh,fcoefs,options)
call allocate_flow_variables(flowvarsAV,mesh,fcoefs,options)
call initialise_flow(options,flowvars,flowvarsAV,mesh,fcoefs,1_in)

!Initialise the segmented flow domain 
allocate(fres_par(options%num_threads))
allocate(fcoefs_par(options%num_threads)) 
allocate(flowvars_par(options%num_threads))
allocate(flowvarsAV_par(options%num_threads))
do tt=1,options%num_threads
    call allocate_flow_variables(flowvars_par(tt),mesh_par(tt),fcoefs_par(tt),options)
    call initialise_flow(options,flowvars_par(tt),flowvarsAV_par(tt),mesh_par(tt),fcoefs_par(tt),tt)
    allocate(fres_par(tt)%W(mesh_par(tt)%ncell,4))
    allocate(fres_par(tt)%W0(mesh_par(tt)%ncell,4))
    allocate(fres_par(tt)%R(mesh_par(tt)%ncell,4))
    allocate(fres_par(tt)%D(mesh_par(tt)%ncell,4))
    allocate(fres_par(tt)%cell_Lp(mesh_par(tt)%ncell,4))
    fres_par(tt)%W(:,:) = 0.0d0 
    fres_par(tt)%W0(:,:) = 0.0d0 
    fres_par(tt)%R(:,:) = 0.0d0 
    fres_par(tt)%D(:,:) = 0.0d0 
    fres_par(tt)%cell_Lp(:,:) = 0.0d0 
    call get_conservative_vars(fres_par(tt)%W,flowvars_par(tt)) 
    call get_primative_vars(fres_par(tt)%W,flowvars_par(tt),mesh_par(tt)%ncell,nanflag) 
    ! print *, '-----------'
    ! print *, 'mesh_par(tt)%ncell = ',mesh_par(tt)%ncell
    ! print *, 'mesh_par(tt)%nedge = ',mesh_par(tt)%nedge
end do 

!Set RK coeficients
do tt=1,options%num_threads
    call set_RK_coeficients(fcoefs_par(tt),options)
end do 

!Preprocess required boundary conditions 
if ((mesh%bc_active(3) == 1) .OR. (mesh%bc_active(4) == 1)) then !Evaluate intake mass flow rate
    call set_initial_massflux(flowvars,mesh,options)
end if 

!Timestepping
if (options%csdisp) then
    write(*,'(A)') '--> timestepping ------ parallel mode:'
end if 
flowvars%res_conv = 0
if (options%csdisp) then
    write(*,'(A,A,A,A,A,A)') '    ittn ','     res ','            cl ','             cd ','          Fsd ','          Rsd '
    write(*,'(A)') '    ----------------------------------------------------------------------------'
end if

!Initialise forces 
fcoefs%cl = 0.0d0  
fcoefs%cd = 0.0d0 
fcoefs%cm = 0.0d0 

!Initialise residual trend 
longres_sum = 0.0d0 
rhores_hist(:) = 0.0d0 
rhores_hist_av(:) = 0.0d0 

!Initialise additional arrays 
if (options%evalMassflux) then !Mass flux through inflow 
    massfluxsten(:) = 0.0d0 
end if 

!Initialise parallel zones 
resconv = 0 
res_sd = 0.0d0
rho_res = 0.0d0  
!$OMP parallel private(tr_idx)

!Get thread index for each zone 
!$OMP critical 
tr_idx = omp_get_thread_num()
tr_idx = tr_idx + 1
tr_idx = mesh_par(tr_idx)%zone_right(1)
!$OMP end critical 

!Push boundary condition properties to each thread
!$OMP barrier
!$OMP single 
if ((mesh%bc_active(3) == 1) .OR. (mesh%bc_active(4) == 1)) then 
    do tt=1,options%num_threads
        flowvars_par(tt)%A_in = flowvars%A_in
        flowvars_par(tt)%mflow_in = flowvars%mflow_in
        flowvars_par(tt)%vnset_of = flowvars%vnset_of 
    end do 
end if 
!$OMP end single 

!Timestep each zone
!$OMP barrier
do fitt=1,options%ittlim

    !Convergence cycle condition
    if (resconv == 1) then
        cycle
    end if

    !NAN value cycle condition 
    if (nanflag == 1) then 
        cycle 
    end if

    !Evaluate cell timesteps
    mesh_par(tr_idx)%specrad(:) = 0.0d0 
    !$OMP barrier
    call cell_spectral_radius(mesh_par,flowvars_par,tr_idx)
    !$OMP barrier
    call timestep_p(mesh_par,options,nanflag,tr_idx)

    !Rk iterate 
    call rk_iterate_p(fres_par,fcoefs_par,mesh_par,flowvars_par,options,rho_res,tr_idx,fitt,nanflag)
    !$OMP barrier

    !Evaluate any additional target properties 
    if (options%evalMassflux) then !Mass flux through inflow 
        call evaluate_inflow_mass_flux(flowvars_par,mesh_par,options,tr_idx)
    end if

    !Process step ---
    !$OMP single 
        
        !Gather force coeficients
        fcoefs%cl = 0.0d0 
        fcoefs%cd = 0.0d0 
        do tt=1,options%num_threads
            fcoefs%cl = fcoefs%cl + fcoefs_par(tt)%cl 
            fcoefs%cd = fcoefs%cd + fcoefs_par(tt)%cd
        end do 

        !Accumulate mass flux through inflow to stencil
        if (options%evalMassflux) then !Mass flux through inflow 
            mfluxin = sum(flowvars_par(:)%mflux_in)/real(options%num_threads,dp)
            massfluxsten(mod(fitt,options%massflux_niterav)+1) = mfluxin
        end if 
        
        !Process step 
        rhores_hist(fitt) = rho_res
        call long_residual_sd(res_sd,rhores_hist,rhores_hist_av,longres_sum,fitt)
        if (res_sd .NE. 0.0d0) then 
            res_sd = log10(abs(res_sd))
        else
            res_sd = res_sd
        end if 
        if (rho_res .NE. 0.0d0) then 
            rho_res = log10(abs(rho_res))
        else
            rho_res = rho_res
        end if 
        call average_fcoefficients(fcoefs,options,fitt)
        flowvars%rho_res = rho_res
        if (fcoefs%force_sd .NE. 0.0d0) then 
            flowvars%force_res = log10(abs(fcoefs%force_sd))
        else
            flowvars%force_res = fcoefs%force_sd
        end if 

        !NAN error handling
        if (nanflag == 1) then
            write(*,'(A)') '    nan value exit'
            fcoefs%cl_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
            fcoefs%cd_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
            fcoefs%cm_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
            options%status = nanflag
        elseif (options%status == 1) then 
            write(*,'(A)') '    invalid options exit'
            fcoefs%cl_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
            fcoefs%cd_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
            fcoefs%cm_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
            options%status = 1
            nanflag = 1
        end if

        !History display
        if (mod(fitt,100) == 0) then
            if (options%csdisp) then
                write(*,'(A,A,A,A,A,A)') '    ittn ','     res ','            cl ','             cd ','          Fsd ',&
                                         '          Rsd '
                write(*,'(A)') '    ----------------------------------------------------------------------------'
            end if
        end if
        if (options%csdisp) then
            write(*,'(I8,A,F10.6,A,F14.8,A,F14.8,A,F10.6,A,F10.6)') fitt,'  ',rho_res,'  ',fcoefs%cl,'  ',fcoefs%cd,'   ',&
                flowvars%force_res,'   ',res_sd
        end if

        !Residual convergence check
        if (fitt .GE. options%ittlim_min) then 
            if (rho_res .LE. options%term_res) then
                flowvars%res_conv = 1
                if (options%csdisp) then
                    write(*,'(A)') '    rho_res convergence'
                end if
                resconv = 1 
            elseif (flowvars%force_res .LE. options%term_res)then 
                flowvars%res_conv = 1
                if (options%csdisp) then
                    write(*,'(A)') '    force_sd convergence'
                end if
                resconv = 1 
            elseif (res_sd .LE. options%term_res) then 
                flowvars%res_conv = 1
                if (options%csdisp) then
                    write(*,'(A)') '    residual stagnation convergence'
                end if
                resconv = 1  
            end if
        end if 
    !$OMP end single 
end do 
!$OMP end parallel 

!Construct final flow variables fields 
do tt=1,options%num_threads
    do cc=1,mesh_par(tt)%ncell
        flowvars%rho(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%rho(cc)
        flowvars%u(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%u(cc)
        flowvars%v(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%v(cc)
        flowvars%e(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%e(cc)
        ! flowvars%p(mesh_par(tt)%cell_link(cc)) = mesh_par(tt)%psensor(cc)
        flowvars%p(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%p(cc)
        flowvars%cp(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%cp(cc)
        flowvars%mach(mesh_par(tt)%cell_link(cc)) = flowvars_par(tt)%mach(cc)
        flowvars%resid(mesh_par(tt)%cell_link(cc)) = fres_par(tt)%R(cc,1)
        mesh%psensor(mesh_par(tt)%cell_link(cc)) = mesh_par(tt)%psensor(cc)
    end do 
end do 
if (options%evalMassflux) then
    flowvars%mflux_in = sum(massfluxsten(:))/real(options%massflux_niterav,dp)
    fcoefs%mflux_in = flowvars%mflux_in
else
    flowvars%mflux_in = 0.0d0 
    fcoefs%mflux_in = 0.0d0 
end if 

!nan error handling
if (nanflag == 1) then
    write(*,'(A)') '    nan value exit'
    fcoefs%cl_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
    fcoefs%cd_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
    fcoefs%cm_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
    fcoefs%mflux_in = ieee_value(1.0d0,IEEE_QUIET_NAN)
    return 
end if

!Display mass flux error
if ((mesh%bc_active(3) == 1) .OR. (mesh%bc_active(4) == 1)) then 
    call mass_flux_error(mesh,flowvars)
    if (options%csdisp) then
        write(*,'(A,F8.4,A)') '    mass flux error = ',flowvars%mflux_error*100.0d0,'%'
    end if 
end if

!Display inflow mass flux 
if (options%evalMassflux) then
    if (options%csdisp) then
        write(*,'(A,F12.8)') '    inflow mass flux = ',flowvars%mflux_in
    end if 
end if 
return 
end subroutine flow_solve




!Runge-Kutta timestep subroutine (parallel) =========================
subroutine rk_iterate_p(fres_par,fcoefs_par,mesh_par,flowvars_par,options,rho_res,tr_idx,iter,nanflag)
implicit none 

!Variables - Import
integer(in) :: nanflag,tr_idx,upd_allow,iter
real(dp) :: rho_res
type(options_data) :: options 
type(mesh_data), dimension(:) :: mesh_par
type(flow_var_data), dimension(:) :: flowvars_par
type(parallel_residual), dimension(:) :: fres_par
type(coeficient_data), dimension(:), allocatable :: fcoefs_par

!Variables - Local 
integer(in) :: tt,ss,rk_stage,vari

!Set current step initial conservative variables
upd_allow = 0 
fres_par(tr_idx)%W0(:,:) = fres_par(tr_idx)%W(:,:)

!Perform itteration
!$OMP barrier
do rk_stage=1,options%nRKstage

    !Initialise this thread
    fres_par(tr_idx)%R(:,:) = 0.0d0 
    fcoefs_par(tr_idx)%cx = 0.0d0
    fcoefs_par(tr_idx)%cy = 0.0d0 
    if (options%RKstagediss(rk_stage) == 1) then 
        fres_par(tr_idx)%D(:,:) = 0.0d0 
        fres_par(tr_idx)%cell_Lp(:,:) = 0.0d0 
        mesh_par(tr_idx)%psensor_num(:) = 0.0d0
        mesh_par(tr_idx)%psensor_dnum(:) = 0.0d0
        mesh_par(tr_idx)%psensor(:) = 0.0d0
    end if
    !$OMP barrier

    !Evaluate dissipation 
    if (options%RKstagediss(rk_stage) == 1) then 

        !Evaluate pressure sensor
        call accumulate_pressure_sensor_p(mesh_par,flowvars_par,tr_idx)
        !$OMP barrier
        mesh_par(tr_idx)%psensor(:) = abs(mesh_par(tr_idx)%psensor_num(:))/mesh_par(tr_idx)%psensor_dnum(:)
        !$OMP barrier
    
        !Smooth pressure sensor
        if (options%psNsmooth .GT. 0) then 
            do ss=1,options%psNsmooth
                mesh_par(tr_idx)%psensor_temp(:) = 0.0d0 
                !$OMP barrier
                call smooth_psensor(mesh_par,tr_idx) 
                !$OMP barrier
                mesh_par(tr_idx)%psensor(:) = abs(mesh_par(tr_idx)%psensor_temp(:) + mesh_par(tr_idx)%psensor(:))/&
                (real(mesh_par(tr_idx)%cell_nedgei(:),dp) + 1.0d0)
            end do 
            !$OMP barrier
        end if 

        !Evaluate cell laplacian
        call accumulate_laplacian_p(fres_par,mesh_par,tr_idx)
        !$OMP barrier
        
        !Evaluate dissipation
        call accumulate_dissipation_p(fres_par,mesh_par,flowvars_par,options,tr_idx)
        !$OMP barrier
    end if

    !Evaluate residual 
    call accumulate_residual_p(fres_par,mesh_par,flowvars_par,fcoefs_par,options,iter,rk_stage,upd_allow,tr_idx)
    !$OMP barrier

    !Incorporate dissipation to residual 
    fres_par(tr_idx)%R(:,:) = fres_par(tr_idx)%R(:,:) - fres_par(tr_idx)%D(:,:) 
    !!$OMP barrier

    !Perform step
    do vari=1,4
        fres_par(tr_idx)%W(:,vari) = fres_par(tr_idx)%W0(:,vari) - fcoefs_par(tr_idx)%RK_alfa(rk_stage)*&
                                    (mesh_par(tr_idx)%deltaT(:)/mesh_par(tr_idx)%cell_vol(:))*fres_par(tr_idx)%R(:,vari)
    end do 
    !!$OMP barrier

    !Update flow primative variables
    call get_primative_vars(fres_par(tr_idx)%W,flowvars_par(tr_idx),mesh_par(tr_idx)%ncell,nanflag)  

    !Syncronise threads
    !$OMP barrier 

    !Exit on nan value 
    if (nanflag == 1) then
        exit  
    end if
end do 

!Calculate log10 density residual 2norm
!$OMP single
rho_res = 0.0d0 
do tt=1,options%num_threads
    rho_res = rho_res + sqrt(sum((fres_par(tt)%R(:,1))**2))
end do 
!$OMP end single
flowvars_par(tr_idx)%rho_res = rho_res
return
end subroutine rk_iterate_p




!Cell residual accumulation subroutine (paralell) =========================
subroutine accumulate_residual_p(fres_par,mesh_par,flowvars_par,fcoefs_par,options,iter,rk_stage,upd_allow,tr_idx)
implicit none 

!Variables - Import
integer(in) :: tr_idx,iter,rk_stage,upd_allow
type(parallel_residual), dimension(:) :: fres_par
type(flow_var_data), dimension(:) :: flowvars_par
type(coeficient_data), dimension(:) :: fcoefs_par
type(mesh_data), dimension(:) :: mesh_par
type(options_data) :: options

!Variables - Local
integer(in) :: ii,tt,cl,cr,zl,zr
real(dp) :: cpw
real(dp) :: Flxe(4)
real(dp) :: u_outflow(mesh_par(tr_idx)%nzone_outflow),v_outflow(mesh_par(tr_idx)%nzone_outflow)
real(dp) :: rho_outflow(mesh_par(tr_idx)%nzone_outflow)

!Evaluate pressure gradients in wall cells 
if (options%wall_bc_type == 'l') then 
    !call evaluate_wadj_pgrad(flowvars,mesh)
end if 

!Evaluate average outflow state in each zone 
if ((mesh_par(tr_idx)%bc_active(7) == 1) .AND. (rk_stage == 1)) then 
! if (mesh_par(tr_idx)%bc_active(7) == 1) then 
    flowvars_par(tr_idx)%u_outflow_t(:) = 0.0d0 
    flowvars_par(tr_idx)%v_outflow_t(:) = 0.0d0 
    flowvars_par(tr_idx)%rho_outflow_t(:) = 0.0d0 
    do ii=1,mesh_par(tr_idx)%nedge
        if (mesh_par(tr_idx)%cell_left(ii) == -7) then
            cr = mesh_par(tr_idx)%cell_right(ii)
            zr = mesh_par(tr_idx)%outflow_zone(ii)
            flowvars_par(tr_idx)%u_outflow_t(zr) = flowvars_par(tr_idx)%u_outflow_t(zr) + &
            (flowvars_par(tr_idx)%u(cr)/real(mesh_par(tr_idx)%nedge_of_zone(zr),dp))
            flowvars_par(tr_idx)%v_outflow_t(zr) = flowvars_par(tr_idx)%v_outflow_t(zr) + &
            (flowvars_par(tr_idx)%v(cr)/real(mesh_par(tr_idx)%nedge_of_zone(zr),dp))
            flowvars_par(tr_idx)%rho_outflow_t(zr) = flowvars_par(tr_idx)%rho_outflow_t(zr) + &
            (flowvars_par(tr_idx)%rho(cr)/real(mesh_par(tr_idx)%nedge_of_zone(zr),dp))
        end if 
    end do
    !$OMP barrier
    u_outflow(:) = 0.0d0 
    v_outflow(:) = 0.0d0 
    rho_outflow(:) = 0.0d0 
    do tt=1,options%num_threads
        do ii=1,mesh_par(tr_idx)%nzone_outflow
            u_outflow(ii) = u_outflow(ii) + flowvars_par(tt)%u_outflow_t(ii)
            v_outflow(ii) = v_outflow(ii) + flowvars_par(tt)%v_outflow_t(ii)
            rho_outflow(ii) = rho_outflow(ii) + flowvars_par(tt)%rho_outflow_t(ii)
        end do 
    end do 
    !$OMP barrier
    flowvars_par(tr_idx)%u_outflow(:) = u_outflow(:)
    flowvars_par(tr_idx)%v_outflow(:) = v_outflow(:)
    flowvars_par(tr_idx)%rho_outflow(:) = rho_outflow(:)
    ! if (iter == 1) then 
    !     flowvars_par(tr_idx)%u_outflow(:) = u_outflow(:)
    !     flowvars_par(tr_idx)%v_outflow(:) = v_outflow(:)
    !     flowvars_par(tr_idx)%rho_outflow(:) = rho_outflow(:)
    ! else
    !     ! flowvars_par(tr_idx)%u_outflow(:) = flowvars_par(tr_idx)%u_outflow(:) + &
    !     ! 0.125d0*(u_outflow(:) - flowvars_par(tr_idx)%u_outflow(:))
    !     ! flowvars_par(tr_idx)%v_outflow(:) = flowvars_par(tr_idx)%v_outflow(:) + &
    !     ! 0.125d0*(v_outflow(:) - flowvars_par(tr_idx)%v_outflow(:))
    !     ! flowvars_par(tr_idx)%rho_outflow(:) = flowvars_par(tr_idx)%rho_outflow(:) + &
    !     ! 0.125d0*(rho_outflow(:) - flowvars_par(tr_idx)%rho_outflow(:))
    ! end if 
    !$OMP barrier
end if 

!Evaluate residual 
do ii=1,mesh_par(tr_idx)%nedge

    !Left and right cells and zones for this edge
    cl = mesh_par(tr_idx)%cell_left(ii)
    cr = mesh_par(tr_idx)%cell_right(ii)
    zl = mesh_par(tr_idx)%zone_left(ii)
    zr = mesh_par(tr_idx)%zone_right(ii)

    !Find edge flux 
    call edge_flux(Flxe,cpw,cl,cr,zl,zr,ii,mesh_par,flowvars_par,options)
    if (mesh_par(zr)%cell_mz(cr) == 0) then 
        fres_par(zr)%R(cr,:) = fres_par(zr)%R(cr,:) - Flxe(:)
    else    
        !$OMP atomic  
        fres_par(zr)%R(cr,1) = fres_par(zr)%R(cr,1) - Flxe(1)
        !$OMP atomic  
        fres_par(zr)%R(cr,2) = fres_par(zr)%R(cr,2) - Flxe(2)
        !$OMP atomic  
        fres_par(zr)%R(cr,3) = fres_par(zr)%R(cr,3) - Flxe(3)
        !$OMP atomic  
        fres_par(zr)%R(cr,4) = fres_par(zr)%R(cr,4) - Flxe(4)
    end if 
    if (cl .GT. 0) then
        if (mesh_par(zl)%cell_mz(cl) == 0) then 
            fres_par(zl)%R(cl,:) = fres_par(zl)%R(cl,:) + Flxe(:)
        else
            !$OMP atomic 
            fres_par(zl)%R(cl,1) = fres_par(zl)%R(cl,1) + Flxe(1)
            !$OMP atomic 
            fres_par(zl)%R(cl,2) = fres_par(zl)%R(cl,2) + Flxe(2)
            !$OMP atomic   
            fres_par(zl)%R(cl,3) = fres_par(zl)%R(cl,3) + Flxe(3)
            !$OMP atomic   
            fres_par(zl)%R(cl,4) = fres_par(zl)%R(cl,4) + Flxe(4)
        end if 
    end if

    !Accumulate force coefficients on cell wall edges
    if (cl == -1) then !(surface inwards normal)
        fcoefs_par(tr_idx)%cx = fcoefs_par(tr_idx)%cx - cpw*mesh_par(tr_idx)%edge_dy(ii)
        fcoefs_par(tr_idx)%cy = fcoefs_par(tr_idx)%cy + cpw*mesh_par(tr_idx)%edge_dx(ii)
    end if
    if (options%status == 1) then 
        return 
    end if
end do
!$OMP barrier

!Enforce mass flux boundary conditions 
if (mesh_par(tr_idx)%bc_active(3) == 1) then     !Mass flux inflow
    ! call mass_flux_inflow_bc(R,mesh,flowvars)
elseif (mesh_par(tr_idx)%bc_active(4) == 1) then !Mass flux outflow
    if (options%mflux_bc_type == 'p') then !pressure diven flux
        ! call mass_flux_outflow_bc_P(R,mesh,flowvars,options,iter)
    elseif (options%mflux_bc_type == 'v') then !velocity driven flux 
        call mass_flux_outflow_bc_V(fres_par(tr_idx)%R,mesh_par,flowvars_par,options,iter,rk_stage,upd_allow,tr_idx)
    end if 
end if 
!$OMP barrier

!Calculate Cl Cd Cm at this iteration
fcoefs_par(tr_idx)%cl = fcoefs_par(tr_idx)%cy*cos(options%aoa_rad) - fcoefs_par(tr_idx)%cx*sin(options%aoa_rad)
fcoefs_par(tr_idx)%cd = fcoefs_par(tr_idx)%cx*cos(options%aoa_rad) + fcoefs_par(tr_idx)%cy*sin(options%aoa_rad)
! fcoefs_par(tr_idx)%cl = fcoefs_par(tr_idx)%cl/mesh%chordx
! fcoefs_par(tr_idx)%cd = fcoefs_par(tr_idx)%cd/mesh%chordx
fcoefs_par(tr_idx)%cm = 0.0d0 
return 
end subroutine accumulate_residual_p




!Cell dissipation accumulation subroutine (paralell) =========================
subroutine accumulate_dissipation_p(fres_par,mesh_par,flowvars_par,options,tr_idx)
implicit none 

!Variables - Import
integer(in) :: tr_idx
type(parallel_residual), dimension(:) :: fres_par
type(flow_var_data), dimension(:) :: flowvars_par
type(mesh_data), dimension(:) :: mesh_par
type(options_data) :: options

!Variables - Local
integer(in) :: ii,cl,cr,zl,zr
real(dp) :: dk2,dk4,dpsilam
real(dp) :: Dedge(4)

!Evaluate dissipation
do ii=1,mesh_par(tr_idx)%nedge

    !Left and right cells and zones for this edge
    cl = mesh_par(tr_idx)%cell_left(ii)
    cr = mesh_par(tr_idx)%cell_right(ii)
    zl = mesh_par(tr_idx)%zone_left(ii)
    zr = mesh_par(tr_idx)%zone_right(ii)

    !Calculate edge dissipation
    call edge_ADcoeff(dk2,dk4,dpsilam,cl,cr,zl,zr,mesh_par,options,flowvars_par,ii,tr_idx)
    call edge_dissipation(Dedge,dK2,dK4,dpsilam,cl,cr,zl,zr,fres_par)
    if (mesh_par(zr)%cell_mz(cr) == 0) then 
        fres_par(zr)%D(cr,:) = fres_par(zr)%D(cr,:) - Dedge(:)
    else    
        !$OMP atomic  
        fres_par(zr)%D(cr,1) = fres_par(zr)%D(cr,1) - Dedge(1)
        !$OMP atomic  
        fres_par(zr)%D(cr,2) = fres_par(zr)%D(cr,2) - Dedge(2)
        !$OMP atomic  
        fres_par(zr)%D(cr,3) = fres_par(zr)%D(cr,3) - Dedge(3)
        !$OMP atomic  
        fres_par(zr)%D(cr,4) = fres_par(zr)%D(cr,4) - Dedge(4)
    end if 
    if (cl .GT. 0) then
        if (mesh_par(zl)%cell_mz(cl) == 0) then 
            fres_par(zl)%D(cl,:) = fres_par(zl)%D(cl,:) + Dedge(:)
        else
            !$OMP atomic   
            fres_par(zl)%D(cl,1) = fres_par(zl)%D(cl,1) + Dedge(1)
            !$OMP atomic   
            fres_par(zl)%D(cl,2) = fres_par(zl)%D(cl,2) + Dedge(2)
            !$OMP atomic   
            fres_par(zl)%D(cl,3) = fres_par(zl)%D(cl,3) + Dedge(3)
            !$OMP atomic   
            fres_par(zl)%D(cl,4) = fres_par(zl)%D(cl,4) + Dedge(4)
        end if 
    end if
end do
return 
end subroutine accumulate_dissipation_p




!Cell pressure sensor accumulation subroutine (paralell) =========================
subroutine accumulate_pressure_sensor_p(mesh_par,flowvars_par,tr_idx)
implicit none 

!Variables - Import
integer(in) :: tr_idx
type(flow_var_data), dimension(:) :: flowvars_par
type(mesh_data), dimension(:) :: mesh_par

!Variables - Local
integer(in) :: ii,cl,cr,zl,zr
real(dp) :: psensor_num,psensor_dnum

!Evaluate dissipation
do ii=1,mesh_par(tr_idx)%nedge

    !Left and right cells and zones for this edge
    cl = mesh_par(tr_idx)%cell_left(ii)
    cr = mesh_par(tr_idx)%cell_right(ii)
    zl = mesh_par(tr_idx)%zone_left(ii)
    zr = mesh_par(tr_idx)%zone_right(ii)

    !Pressure sensor value 
    call edge_pressure_sensor(psensor_num,psensor_dnum,cl,cr,zl,zr,ii,flowvars_par,mesh_par)
    if (mesh_par(zr)%cell_mz(cr) == 0) then 
        mesh_par(zr)%psensor_num(cr) = mesh_par(zr)%psensor_num(cr) - psensor_num
        mesh_par(zr)%psensor_dnum(cr) = mesh_par(zr)%psensor_dnum(cr) + psensor_dnum
    else    
        !$OMP atomic 
        mesh_par(zr)%psensor_num(cr) = mesh_par(zr)%psensor_num(cr) - psensor_num
        !$OMP atomic  
        mesh_par(zr)%psensor_dnum(cr) = mesh_par(zr)%psensor_dnum(cr) + psensor_dnum
    end if 
    if (cl .GT. 0) then
        if (mesh_par(zl)%cell_mz(cl) == 0) then 
            mesh_par(zl)%psensor_num(cl) = mesh_par(zl)%psensor_num(cl) + psensor_num
            mesh_par(zl)%psensor_dnum(cl) = mesh_par(zl)%psensor_dnum(cl) + psensor_dnum
        else
            !$OMP atomic  
            mesh_par(zl)%psensor_num(cl) = mesh_par(zl)%psensor_num(cl) + psensor_num
            !$OMP atomic  
            mesh_par(zl)%psensor_dnum(cl) = mesh_par(zl)%psensor_dnum(cl) + psensor_dnum
        end if 
    end if 
end do
return 
end subroutine accumulate_pressure_sensor_p




!Cell laplacian accumulation subroutine (paralell) =========================
subroutine accumulate_laplacian_p(fres_par,mesh_par,tr_idx)
implicit none 

!Variables - Import
integer(in) :: tr_idx
type(parallel_residual), dimension(:) :: fres_par
type(mesh_data), dimension(:) :: mesh_par

!Variables - Local
integer(in) :: ii,cl,cr,zl,zr
real(dp) :: Lp_edge(4)

!Evaluate dissipation
do ii=1,mesh_par(tr_idx)%nedge

    !Left and right cells and zones for this edge
    cl = mesh_par(tr_idx)%cell_left(ii)
    cr = mesh_par(tr_idx)%cell_right(ii)
    zl = mesh_par(tr_idx)%zone_left(ii)
    zr = mesh_par(tr_idx)%zone_right(ii)

    !Cell laplacian value 
    call edge_laplacian(Lp_edge,cl,cr,zl,zr,fres_par)
    if (mesh_par(zr)%cell_mz(cr) == 0) then 
        fres_par(zr)%cell_Lp(cr,:) = fres_par(zr)%cell_Lp(cr,:) - Lp_edge(:)
    else    
        !$OMP atomic  
        fres_par(zr)%cell_Lp(cr,1) = fres_par(zr)%cell_Lp(cr,1) - Lp_edge(1)
        !$OMP atomic  
        fres_par(zr)%cell_Lp(cr,2) = fres_par(zr)%cell_Lp(cr,2) - Lp_edge(2)
        !$OMP atomic  
        fres_par(zr)%cell_Lp(cr,3) = fres_par(zr)%cell_Lp(cr,3) - Lp_edge(3)
        !$OMP atomic  
        fres_par(zr)%cell_Lp(cr,4) = fres_par(zr)%cell_Lp(cr,4) - Lp_edge(4)
    end if 
    if (cl .GT. 0) then
        if (mesh_par(zl)%cell_mz(cl) == 0) then 
            fres_par(zl)%cell_Lp(cl,:) = fres_par(zl)%cell_Lp(cl,:) + Lp_edge(:)
        else
            !$OMP atomic  
            fres_par(zl)%cell_Lp(cl,1) = fres_par(zl)%cell_Lp(cl,1) + Lp_edge(1)
            !$OMP atomic  
            fres_par(zl)%cell_Lp(cl,2) = fres_par(zl)%cell_Lp(cl,2) + Lp_edge(2)
            !$OMP atomic  
            fres_par(zl)%cell_Lp(cl,3) = fres_par(zl)%cell_Lp(cl,3) + Lp_edge(3)
            !$OMP atomic  
            fres_par(zl)%cell_Lp(cl,4) = fres_par(zl)%cell_Lp(cl,4) + Lp_edge(4)
        end if 
    end if
end do
return 
end subroutine accumulate_laplacian_p




!Timestep calculation subroutine =========================
subroutine timestep_p(mesh_par,options,nanflag,tr_idx)
implicit none

!Variables - Import
integer(in) :: nanflag,tr_idx
type(options_data) :: options 
type(mesh_data), dimension(:) :: mesh_par

!Variables - Local 
integer(in) :: ii

!Calculate timestep for each cell
do ii=1,mesh_par(tr_idx)%ncell
    mesh_par(tr_idx)%deltaT(ii) = options%cfl*(mesh_par(tr_idx)%cell_vol(ii)/mesh_par(tr_idx)%specrad(ii))
    if (isnan(mesh_par(tr_idx)%deltaT(ii))) then
        print *,  'nan dt -> thread: ',tr_idx
        print *, 'lam = ',mesh_par(tr_idx)%specrad(ii)
        nanflag = 1
    end if
end do

! !Set global timestep if selected
! if (options%ts_mode == 0) then 
!     mesh%deltaT(:) = minval(mesh%deltaT(:))
! end if
return 
end subroutine timestep_p


end module flow_iterator_mod