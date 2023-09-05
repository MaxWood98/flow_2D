! !Flow solver solution itteration module (single threded)
! !Max Wood (mw16116@bristol.ac.uk)
! !University of Bristol - Department of Aerospace Engineering
! !Version: 0.6.4
! !Updated: 21/04/23

! !Module
! module flow_iterator_single_mod
! use flow_io_mod
! use flow_data_mod
! use flow_general_mod
! use flow_calculation_mod
! contains 

! !Flow timestepping subroutine =========================
! subroutine flow_solve_single(mesh,options,flowvars,flowvarsAV,fcoefs)
! use ieee_arithmetic, only: ieee_value,IEEE_QUIET_NAN
! implicit none 

! !Variables - Import
! type(options_data) :: options
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars,flowvarsAV
! type(coeficient_data) :: fcoefs

! !Variables - Local 
! integer(in) :: fitt,nanflag,upd_allow
! real(dp) :: rho_res,longres_sum,res_sd
! real(dp) :: W(mesh%ncell,4),W0(mesh%ncell,4),R(mesh%ncell,4),D(mesh%ncell,4),cell_Lp(mesh%ncell,4)
! real(dp) :: rhores_hist(options%ittlim),rhores_hist_av(options%ittlim)

! !Inititalise status 
! options%status = 0 

! !Allocate flow variables 
! call allocate_flow_variables(flowvars,mesh,fcoefs,options)
! call allocate_flow_variables(flowvarsAV,mesh,fcoefs,options)

! !Initialise flow domain 
! call initialise_flow(options,flowvars,flowvarsAV,mesh,fcoefs,1_in)
! call get_conservative_vars(W,flowvars) 
! W0(:,:) = W(:,:) 
! R(:,:) = 0.0d0 
! D(:,:) = 0.0d0 
! cell_Lp(:,:) = 0.0d0 

! !Initialise residual trend 
! longres_sum = 0.0d0 
! rhores_hist(:) = 0.0d0 
! rhores_hist_av(:) = 0.0d0 

! !Update flow primative variables 
! call get_primative_vars(W,flowvars,mesh%ncell,nanflag) 

! !Set RK coeficients
! call set_RK_coeficients(fcoefs,options)

! !Preprocess required mass flux boundary conditions 
! if ((mesh%bc_active(3) == 1) .OR. (mesh%bc_active(4) == 1)) then !Evaluate intake mass flow rate
!     call set_initial_massflux(flowvars,mesh,options)
! end if 

! !Set to allow boundary condition updates
! upd_allow = 1

! !Timestepping
! if (options%csdisp) then
!     write(*,'(A)') '--> timestepping ------ single mode:'
! end if
! nanflag = 0
! flowvars%res_conv = 0
! if (options%csdisp) then
!     write(*,'(A,A,A,A,A,A)') '    ittn ','     res ','            cl ','             cd ','          Fsd ','          Rsd '
!     write(*,'(A)') '    ----------------------------------------------------------------------------'
! end if
! do fitt=1,options%ittlim

!     !Evaluate cell timesteps
!     call cell_spectral_radius(mesh,flowvars)
!     call timestep_s(mesh,flowvars,options,nanflag)

!     !Set to disallow boundary condition updates if the final iteration
!     if (fitt == options%ittlim) then 
!         upd_allow = 0
!     end if 

!     !Rk iterate 
!     call rk_iterate_s(W,W0,R,D,cell_Lp,mesh,options,flowvars,fcoefs,rho_res,nanflag,fitt,upd_allow)
    
!     !Process step 
!     rhores_hist(fitt) = rho_res
!     call long_residual_sd(res_sd,rhores_hist,rhores_hist_av,longres_sum,fitt)
!     if (res_sd .NE. 0.0d0) then 
!         res_sd = log10(abs(res_sd))
!     else
!         res_sd = res_sd
!     end if 
!     if (rho_res .NE. 0.0d0) then 
!         rho_res = log10(abs(rho_res))
!     else
!         rho_res = rho_res
!     end if 
!     call average_fcoefficients(fcoefs,options,fitt)
!     flowvars%rho_res = rho_res
!     if (fcoefs%force_sd .NE. 0.0d0) then 
!         flowvars%force_res = log10(abs(fcoefs%force_sd))
!     else
!         flowvars%force_res = fcoefs%force_sd
!     end if 

!     !NAN error handling
!     if (nanflag == 1) then
!         write(*,'(A)') '    nan value exit'
!         fcoefs%cl_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
!         fcoefs%cd_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
!         fcoefs%cm_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
!         options%status = nanflag
!         return  
!     elseif (options%status == 1) then 
!         write(*,'(A)') '    invalid options exit'
!         fcoefs%cl_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
!         fcoefs%cd_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
!         fcoefs%cm_av = ieee_value(1.0d0,IEEE_QUIET_NAN)
!         options%status = 1
!         return  
!     end if

!     !History display
!     if (mod(fitt,100) == 0) then
!         if (options%csdisp) then
!             write(*,'(A,A,A,A,A,A)') '    ittn ','     res ','            cl ','             cd ','          Fsd ','          Rsd '
!             write(*,'(A)') '    ----------------------------------------------------------------------------'
!         end if
!     end if
!     if (options%csdisp) then
!         write(*,'(I8,A,F10.6,A,F14.8,A,F14.8,A,F10.6,A,F10.6)') fitt,'  ',rho_res,'  ',fcoefs%cl,'  ',fcoefs%cd,'   ',&
!              flowvars%force_res,'   ',res_sd
!     end if 

!     !Residual convergence check
!     if (fitt .GE. options%ittlim_min) then 
!         if (rho_res .LE. options%term_res) then
!             flowvars%res_conv = 1
!             if (options%csdisp) then
!                 write(*,'(A)') '    rho_res convergence'
!             end if
!             exit 
!         elseif ((flowvars%force_res .LE. options%term_res) .AND. (fitt .GE. 10)) then 
!             flowvars%res_conv = 1
!             if (options%csdisp) then
!                 write(*,'(A)') '    force_sd convergence'
!             end if
!             exit 
!         elseif ((res_sd .LE. options%term_res) .AND. (fitt .GE. 10)) then 
!             flowvars%res_conv = 1
!             if (options%csdisp) then
!                 write(*,'(A)') '    residual stagnation convergence'
!             end if
!             exit 
!         end if
!     end if 
! end do 

! !Set to disallow boundary condition updates
! upd_allow = 0

! !Perform flow property averaging 
! if (options%avstencil_size .GT. 0) then !Perform averaging

!     !Initialise variables 
!     flowvarsAV%rho(:) = 0.0d0 
!     flowvarsAV%u(:) = 0.0d0 
!     flowvarsAV%v(:) = 0.0d0 
!     flowvarsAV%p(:) = 0.0d0 
!     flowvarsAV%e(:) = 0.0d0 
!     flowvarsAV%mach(:) = 0.0d0 
!     flowvarsAV%cp(:) = 0.0d0 
!     flowvarsAV%resid(:) = 0.0d0 
!     flowvarsAV%rho_res = 0.0d0 
!     do fitt=1,options%avstencil_size

!         !Evaluate cell timesteps
!         call cell_spectral_radius(mesh,flowvars)
!         call timestep_s(mesh,flowvars,options,nanflag)

!         !Rk iterate 
!         call rk_iterate_s(W,W0,R,D,cell_Lp,mesh,options,flowvars,fcoefs,rho_res,nanflag,fitt,upd_allow)

!         !Process step 
!         rhores_hist(fitt) = rho_res
!         call long_residual_sd(res_sd,rhores_hist,rhores_hist_av,longres_sum,fitt)
!         if (res_sd .NE. 0.0d0) then 
!             res_sd = log10(abs(res_sd))
!         else
!             res_sd = res_sd
!         end if 
!         if (rho_res .NE. 0.0d0) then 
!             rho_res = log10(abs(rho_res))
!         else
!             rho_res = rho_res
!         end if 
!         call average_fcoefficients(fcoefs,options,fitt)
!         flowvars%rho_res = rho_res
!         if (fcoefs%force_sd .NE. 0.0d0) then 
!             flowvars%force_res = log10(abs(fcoefs%force_sd))
!         else
!             flowvars%force_res = fcoefs%force_sd
!         end if 

!         !History display
!         if (mod(fitt,100) == 0) then
!             if (options%csdisp) then
!                 write(*,'(A,A,A,A,A,A)') '    ittn ','     res ','            cl ','             cd ','          Fsd ',&
!                                          '          Rsd '
!                 write(*,'(A)') '    ----------------------------------------------------------------------------'
!             end if
!         end if
!         if (options%csdisp) then
!             write(*,'(I8,A,F10.6,A,F14.8,A,F14.8,A,F10.6,A,F10.6)') fitt,'  ',rho_res,'  ',fcoefs%cl,'  ',fcoefs%cd,'   ',&
!                 flowvars%force_res,'   ',res_sd
!         end if

!         !Accumulate average variables 
!         flowvarsAV%rho(:) = flowvarsAV%rho(:) + flowvars%rho(:) 
!         flowvarsAV%u(:) = flowvarsAV%u(:) + flowvars%u(:)  
!         flowvarsAV%v(:) = flowvarsAV%v(:) + flowvars%v(:)  
!         flowvarsAV%p(:) = flowvarsAV%p(:) + flowvars%p(:)  
!         flowvarsAV%e(:) = flowvarsAV%e(:) + flowvars%e(:)  
!         flowvarsAV%mach(:) = flowvarsAV%mach(:) + flowvars%mach(:)  
!         flowvarsAV%cp(:) = flowvarsAV%cp(:) + flowvars%cp(:)  
!         flowvarsAV%resid(:) = flowvarsAV%resid(:) + R(:,1)
!         flowvarsAV%rho_res = flowvarsAV%rho_res + flowvars%rho_res
!     end do 

!     !Average
!     flowvarsAV%rho(:) = flowvarsAV%rho(:)/real(options%avstencil_size,dp)
!     flowvarsAV%u(:) = flowvarsAV%u(:)/real(options%avstencil_size,dp) 
!     flowvarsAV%v(:) = flowvarsAV%v(:)/real(options%avstencil_size,dp) 
!     flowvarsAV%p(:) = flowvarsAV%p(:)/real(options%avstencil_size,dp) 
!     flowvarsAV%e(:) = flowvarsAV%e(:)/real(options%avstencil_size,dp)
!     flowvarsAV%mach(:) = flowvarsAV%mach(:)/real(options%avstencil_size,dp)
!     flowvarsAV%cp(:) = flowvarsAV%cp(:)/real(options%avstencil_size,dp)
!     flowvarsAV%resid(:) = flowvarsAV%resid(:)/real(options%avstencil_size,dp)  
!     flowvarsAV%rho_res = flowvarsAV%rho_res/real(options%avstencil_size,dp)  
! else !Copy flow properties 
!     flowvarsAV%rho(:) = flowvars%rho(:) 
!     flowvarsAV%u(:) = flowvars%u(:)  
!     flowvarsAV%v(:) = flowvars%v(:)  
!     flowvarsAV%p(:) = flowvars%p(:)  
!     flowvarsAV%e(:) = flowvars%e(:)  
!     flowvarsAV%mach(:) = flowvars%mach(:)  
!     flowvarsAV%cp(:) = flowvars%cp(:)  
!     flowvarsAV%resid(:) = R(:,1) 
!     flowvarsAV%rho_res = flowvars%rho_res
! end if 

! !Display mass flux error
! if ((mesh%bc_active(3) == 1) .OR. (mesh%bc_active(4) == 1)) then 
!     call mass_flux_error(mesh,flowvars)
!     if (options%csdisp) then
!         write(*,'(A,F8.4,A)') '    mass flux error = ',flowvars%mflux_error*100.0d0,'%'
!     end if 
! end if

! !Store residual 
! flowvars%resid(:) = R(:,1)

! !Set status 
! options%status = nanflag
! return 
! end subroutine flow_solve_single




! !Runge-Kutta timestep subroutine (single) =========================
! subroutine rk_iterate_s(W,W0,R,D,cell_Lp,mesh,options,flowvars,fcoefs,rho_res,nanflag,iter,upd_allow)
! implicit none 

! !Variables - Import
! integer(in) :: nanflag,iter,upd_allow
! real(dp) :: rho_res
! real(dp), dimension(:,:) :: W,W0,R,D,cell_Lp
! type(options_data) :: options 
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars
! type(coeficient_data) :: fcoefs

! !Variables - Local 
! integer(in) :: rk_stage,vari

! !Set current step initial conservative variables
! W0(:,:) = W(:,:)

! !Perform itteration
! nanflag = 0
! do rk_stage=1,options%nRKstage

!     !Cell residuals at current state 
!     call cell_spectral_radius(mesh,flowvars)
!     call accumulate_residual_s(R,D,W,cell_Lp,mesh,flowvars,options,fcoefs,options%RKstagediss(rk_stage),iter,rk_stage,upd_allow)

!     !Perform RK stage step
!     do vari=1,4
!         W(:,vari) = W0(:,vari) - fcoefs%RK_alfa(rk_stage)*(mesh%deltaT(:)/mesh%cell_vol(:))*R(:,vari)
!     end do 

!     !Update flow primative variables 
!     call get_primative_vars(W,flowvars,mesh%ncell,nanflag) 

!     !Exit on nan value 
!     if (nanflag == 1) then
!         return 
!     end if
!     if (options%status == 1) then 
!         return 
!     end if
! end do 

! !Calculate log10 density residual 2norm
! ! rho_res = log10(sqrt(sum((W(:,1) - W0(:,1))**2)))
! rho_res = sqrt(sum((R(:,1))**2))
! return
! end subroutine rk_iterate_s




! !Cell residual accumulation subroutine =========================
! subroutine accumulate_residual_s(R,D,W,cell_Lp,mesh,flowvars,options,fcoefs,disstoggle,iter,rk_stage,upd_allow)
! implicit none 

! !Variables - Import
! integer(in) :: disstoggle,iter,rk_stage,upd_allow
! real(dp), dimension(:,:) :: R,D,W,cell_Lp
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars
! type(options_data) :: options
! type(coeficient_data) :: fcoefs 

! !Variables - Local
! integer(in) :: ii,ss,cl,cr
! real(dp) :: psensor_num,psensor_dnum,cpw,dk2,dk4,dpsilam
! real(dp) :: Lp_edge(4),Dedge(4),Flxe(4)

! !Evaluate pressure gradients in wall cells 
! if (options%wall_bc_type == 'l') then 
!     call evaluate_wadj_pgrad(flowvars,mesh)
! end if 

! !Calculate cell residual 
! !D(:,:) = 0.0d0 
! R(:,:) = 0.0d0 
! fcoefs%cx = 0.0d0
! fcoefs%cy = 0.0d0 
! if (disstoggle == 1) then !dissipation 

!     !Precalculate dissipation parameters
!     D(:,:) = 0.0d0 
!     cell_Lp(:,:) = 0.0d0 
!     mesh%psensor_num(:) = 0.0d0
!     mesh%psensor_dnum(:) = 0.0d0
!     do ii=1,mesh%nedge

!         !Left and right cells for this edge
!         cl = mesh%cell_left(ii)
!         cr = mesh%cell_right(ii)
        
!         !Pressure sensor value ------------------------------------------------
!         call edge_pressure_sensor(psensor_num,psensor_dnum,cl,cr,flowvars)
!         mesh%psensor_num(cr) = mesh%psensor_num(cr) - psensor_num
!         mesh%psensor_dnum(cr) = mesh%psensor_dnum(cr) + psensor_dnum
!         if (cl .GT. 0) then
!             mesh%psensor_num(cl) = mesh%psensor_num(cl) + psensor_num
!             mesh%psensor_dnum(cl) = mesh%psensor_dnum(cl) + psensor_dnum
!         end if

!         !Cell laplacian value  ------------------------------------------------
!         call edge_laplacian(Lp_edge,cl,cr,W)
!         cell_Lp(cr,:) = cell_Lp(cr,:) - Lp_edge(:)
!         if (cl .GT. 0) then
!             cell_Lp(cl,:) = cell_Lp(cl,:) + Lp_edge(:)
!         end if
!         if (options%status == 1) then 
!             return 
!         end if
!     end do
!     mesh%psensor(:) = abs(mesh%psensor_num(:))/mesh%psensor_dnum(:)
!     do ss=1,options%psNsmooth
!         call smooth_psensor_s(mesh)
!     end do 

!     !Accumulate cell residuals 
!     do ii=1,mesh%nedge

!         !Left and right cells for this edge
!         cl = mesh%cell_left(ii)
!         cr = mesh%cell_right(ii)
        
!         !Edge flux ------------------------------------------------
!         !Find edge flux 
!         call edge_flux(Flxe,cpw,cl,cr,ii,mesh,flowvars,options)
!         R(cr,:) = R(cr,:) - Flxe(:)
!         if (cl .GT. 0) then
!             R(cl,:) = R(cl,:) + Flxe(:)
!         end if
    
!         !Accumulate force coefficients on cell wall edges
!         if (cl == -1) then !(surface inwards normal)
!             fcoefs%cx = fcoefs%cx - cpw*mesh%edge_dy(ii)
!             fcoefs%cy = fcoefs%cy + cpw*mesh%edge_dx(ii)
!         end if

!         !Dissipation ----------------------------------------------
!         !Calculate edge dissipation
!         call edge_ADcoeff(dk2,dk4,dpsilam,cl,cr,mesh,options,flowvars,ii)
!         call edge_dissipation(Dedge,dK2,dK4,dpsilam,cell_Lp,cl,cr,W)
!         D(cr,:) = D(cr,:) - Dedge(:)
!         if (cl .GT. 0) then
!             D(cl,:) = D(cl,:) + Dedge(:)
!         end if
!         if (options%status == 1) then 
!             return 
!         end if
!     end do

!     !Incorporate dissipation to residual 
!     R(:,:) = R(:,:) - D(:,:)
! else !no dissipation 
!     do ii=1,mesh%nedge

!         !Left and right cells for this edge
!         cl = mesh%cell_left(ii)
!         cr = mesh%cell_right(ii)
        
!         !Edge flux ------------------------------------------------
!         !Find edge flux 
!         call edge_flux(Flxe,cpw,cl,cr,ii,mesh,flowvars,options)
!         R(cr,:) = R(cr,:) - Flxe(:)
!         if (cl .GT. 0) then
!             R(cl,:) = R(cl,:) + Flxe(:)
!         end if
    
!         !Accumulate force coefficients on cell wall edges
!         if (cl == -1) then !(surface inwards normal)
!             fcoefs%cx = fcoefs%cx - cpw*mesh%edge_dy(ii)
!             fcoefs%cy = fcoefs%cy + cpw*mesh%edge_dx(ii)
!         end if
!         if (options%status == 1) then 
!             return 
!         end if
!     end do

!     !Incorporate dissipation to residual 
!     R(:,:) = R(:,:) - D(:,:)
! end if 

! !Enforce mass flux boundary conditions 
! if (mesh%bc_active(3) == 1) then     !Mass flux inflow
!     call mass_flux_inflow_bc(R,mesh,flowvars)
! elseif (mesh%bc_active(4) == 1) then !Mass flux outflow
!     if (options%mflux_bc_type == 'p') then !pressure diven flux
!         call mass_flux_outflow_bc_P(R,mesh,flowvars,options,iter)
!     elseif (options%mflux_bc_type == 'v') then !velocity driven flux 
!         call mass_flux_outflow_bc_V(R,mesh,flowvars,options,iter,rk_stage,upd_allow)
!     end if 
! end if 

! !Calculate Cl Cd Cm at this iteration
! fcoefs%cl = fcoefs%cy*cos(options%aoa_rad) - fcoefs%cx*sin(options%aoa_rad)
! fcoefs%cd = fcoefs%cx*cos(options%aoa_rad) + fcoefs%cy*sin(options%aoa_rad)
! fcoefs%cl = fcoefs%cl/mesh%chordx
! fcoefs%cd = fcoefs%cd/mesh%chordx
! fcoefs%cm = 0.0d0 
! return 
! end subroutine accumulate_residual_s




! !Timestep calculation subroutine =========================
! subroutine timestep_s(mesh,flowvars,options,nanflag)
! implicit none

! !Variables - Import
! integer(in) :: nanflag
! type(options_data) :: options 
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars

! !Variables - Local 
! integer(in) :: ii,cl,cr

! !Calculate timestep for each cell
! call cell_spectral_radius(mesh,flowvars)
! do ii=1,mesh%ncell
!     mesh%deltaT(ii) = options%cfl*(mesh%cell_vol(ii)/mesh%specrad(ii))
!     if (isnan(mesh%deltaT(ii))) then
!         print *,  'nan dt'
!         print *, 'lam = ',mesh%specrad(ii)
!         nanflag = 1
!         return 
!     end if
! end do
! if (options%damptsteps) then 
!     do ii=1,mesh%nedge
!         cl = mesh%cell_left(ii)
!         cr = mesh%cell_right(ii)
!         if ((cl .GT. 0) .AND.(cr .GT. 0)) then 
!             mesh%deltaTD(cl) = min(mesh%deltaT(cl),mesh%deltaT(cr))
!             mesh%deltaTD(cr) = min(mesh%deltaT(cl),mesh%deltaT(cr))
!         end if 
!     end do 
!     mesh%deltaT(:) = mesh%deltaTD(:)
! end if 

! !Set global timestep if selected
! if (options%ts_mode == 0) then 
!     mesh%deltaT(:) = minval(mesh%deltaT(:))
! end if
! return 
! end subroutine timestep_s


! end module flow_iterator_single_mod