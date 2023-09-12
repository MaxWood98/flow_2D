!Flow solver data types and specifications module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.5.8
!Updated: 14/08/23

!Module
module flow_data_mod

!Integer data type 
use ISO_FORTRAN_ENV, only: in=>int64
use ISO_FORTRAN_ENV, only: in32=>int32

!Double precision data type 
use ISO_FORTRAN_ENV, only: dp=>real64 

!Define constant pi 
real(dp) :: pi = 4.0d0*atan(1.0d0)

!Mesh type
type mesh_data
    integer(in) :: ncell,nedge,nvtx,ncshare,npwallcell,npwalledge,nzone_outflow
    integer(in) :: bc_active(7)
    real(dp) :: chordx,mindt
    integer(in), dimension(:), allocatable :: edge_1,edge_2,edgeonwcell,cellonwall,cell_mz
    integer(in), dimension(:), allocatable :: cell_left,cell_right,cell_nedge,cell_nedgei,cell_nedgeiw,cell_link 
    integer(in), dimension(:), allocatable :: zone_left,zone_right,outflow_zone,nedge_of_zone
    integer(in), dimension(:,:), allocatable :: celltransfer
    real(dp), dimension(:), allocatable :: vtx_x,vtx_y,cell_vol,edge_dx,edge_dy,cellcenx,cellceny,cell_elentot,cell_elenint
    real(dp), dimension(:), allocatable :: edge_nx,edge_ny,edgelen,edge_mx,edge_my,deltaT,deltaTD
    real(dp), dimension(:), allocatable :: specrad,psensor,psensor_num,psensor_dnum,psensor_temp
end type mesh_data

!Flow options type
type options_data
    character(len=:), allocatable :: iopath,optpath
    character(len=1) :: mflux_bc_type,wall_bc_type,ff_bc_type,init_state,ffexport_scale
    logical :: csdisp,pptoggle,ffexporttoggle,damptsteps,evalMassflux,force_fixed_pratio
    integer(in) :: ittlim,ittlim_min,avstencil_size,num_threads,nRKstage,ts_mode,psNsmooth
    integer(in) :: status,Mflux_relaxiter,massflux_niterav
    integer(in) :: RKstagediss(4)
    real(dp) :: term_res,aoa,aoa_rad,cfl,k2,c2,k4,sp
    real(dp) :: rhoinf,tinf,machinf,gamma,R,Mflux_relax,Mflux_toll,backp_pratio
end type options_data

!Flow variables type
type flow_var_data
    integer(in) :: res_conv
    real(dp) :: Winf(4)
    real(dp) :: uinf,vinf,velinf,t0inf,p0inf,rho0inf,mflow_in,A_in,A_out,mflux_in!,dt_in
    real(dp) :: pinf,cinf,rhoinf,tinf,machinf
    real(dp) :: mflux_pset,mflux_pset0,backp_pset
    real(dp) :: rho_res,force_res,R,gam,vnset_of,mflux_error!,mass_intake_total,time_total
    real(dp), dimension(:), allocatable :: rho,u,v,p,e,mach,cp,resid,vorticity,vel,temp,paverage
    real(dp), dimension(:), allocatable :: u_outflow,v_outflow,rho_outflow
    real(dp), dimension(:), allocatable :: u_outflow_t,v_outflow_t,rho_outflow_t
    real(dp), dimension(:,:), allocatable :: pgrad
end type flow_var_data

!Force coeficients type
type coeficient_data
    real(dp) :: cx,cy,cl,cd,cm,mflux_in,force_sd
    real(dp) :: cl_av,cd_av,cm_av
    real(dp) :: RK_alfa(4)
    real(dp), dimension(:), allocatable :: cl_stencil,cd_stencil,cm_stencil
end type coeficient_data

!Parallel residual type
type parallel_residual
    real(dp), dimension(:,:), allocatable :: W,W0,R,D,cell_Lp
end type parallel_residual

! !Parallel accumulator type
! type paralell_acc
!     integer(in) :: nparr
!     integer(in), dimension(:), allocatable :: n_edge_pr
!     integer(in), dimension(:,:), allocatable :: edge_se
!     type(parR_cell_map), dimension(:), allocatable :: par_map
!     type(accum), dimension(:), allocatable :: par_accumulator
! end type paralell_acc

! !Parallel cell map type
! type parR_cell_map
!     integer(in), dimension(:), allocatable :: cell2acc_map
!     integer(in), dimension(:), allocatable :: acc2cell_map
! end type parR_cell_map

! !Accumulator type
! type accum
!     real(dp), dimension(:,:), allocatable :: accumr
! end type accum

! !Cell parallel structure type
! type cparstruct
!     integer(in) :: nparzone
!     integer(in), dimension(:), allocatable :: Ncell_zone,Nedge_zone
!     real(dp), dimension(:), allocatable :: rhores_zone
!     type(cparz), dimension(:), allocatable :: par_zone
! end type cparstruct

! !Cell parallel zone type
! type cparz
!     integer(in), dimension(:), allocatable :: parZcell,parZedge,parZcallow
!     integer(in), dimension(:,:), allocatable :: parZe2c
!     real(dp) :: parZ_cl_cd_cm(3),parZ_cx_cy(3)
! end type cparz

end module flow_data_mod