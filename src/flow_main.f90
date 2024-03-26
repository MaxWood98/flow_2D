!Flow Solver Main Program
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering
!Version: 2.5.2
!Updated: 23/03/24

!Main
program flow
use flow_io_mod
use flow_mesh_mod
use flow_general_mod
use flow_iterator_mod
implicit none 

!Variables
integer(in) :: invalid_mesh
type(options_data) :: options 
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars,flowvarsAV
type(coeficient_data) :: fcoefs

!Det default paths 
call set_default_paths(options)

!Set default options 
call set_default_options(options)

!Get command line arguments 
call get_process_arguments(options)

!Import options
call flow_import_options(options)

!Display
if (options%csdisp) then
    write(*,*)
    write(*,*)'+--------------------------------------------+'
    write(*,*)'|        Flow2D - 2D Euler Flow Solver       |'
    write(*,*)'|         Version 2.6.1 || 26/03/2024        |'
    write(*,*)'|                 Max Wood                   |'
    write(*,*)'|           University of Bristol            |'
    write(*,*)'|    Department of Aerospace Engineering     |'
    write(*,*)'+--------------------------------------------+'
    write(*,*)
end if 

!Mesh data import
call flow_import_mesh(mesh,options)

!Preprocess mesh
call flow_preprocess_mesh(mesh,options,invalid_mesh,1_in)
if (invalid_mesh == 1) then 
    options%status = 2
    call export_force_coefs('fcoefs.dat',options,fcoefs)
    flowvars%rho_res = 0.0d0 
    flowvarsAV%rho_res = 0.0d0 
    call export_status(options,flowvars,flowvarsAV)
    write(*,'(A)') '--> ** invalid mesh **'
    stop 
end if

!Solve flow field
call flow_solve(mesh,options,flowvars,flowvarsAV,fcoefs)

!Export restart file 
if (options%status == 0) then 
    call write_restart_file('flow2d_restart.dat',flowvars,mesh,options)
end if 

!Scale flow field
if (options%ffexport_scale == 'real') then 
    call set_flowvars_scale(flowvars,options)
    call set_flowvars_scale(flowvarsAV,options)
elseif (options%ffexport_scale == 'normalised') then 
    !do nothing 
else
    write(*,'(A)') '--> ** invalid export scaling selected - defaulting to normalised **'
end if 

!Evaluate flowfield properties 
call evaluate_flowfield_properties(flowvars,mesh)
! call evaluate_flowfield_properties(flowvarsAV,mesh)

!Export flow coeficients 
call export_force_coefs('fcoefs.dat',options,fcoefs)

!Export result TECPLOT file
if (options%pptoggle) then 
    call saveVolPLT('flow_volume.plt',mesh,options,flowvars)
    ! call saveVolPLT('flow_volumeAV.plt',mesh,options,flowvarsAV)
end if 

!Export flow field raw data 
if (options%ffexporttoggle) then 
    call export_flowfield('flowfield.dat',flowvars,mesh,options)
    ! call export_flowfield('flowfieldAV.dat',flowvarsAV,mesh,options)
end if 

!Export status 
call export_status(options,flowvars,flowvarsAV)

! !Exit display
if (options%csdisp) then
    print *, 'COMPLETE'
    stop
end if 
end program flow