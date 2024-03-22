!Flow solver io module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.5.1
!Updated: 20/03/24

!Module
module flow_io_mod
use io_utilities
use flow_general_mod
contains 


!Read command arguments subroutine ===========================
subroutine get_process_arguments(options)
implicit none

!Variables - Import
type(options_data) :: options
    
!Variables - Local 
integer(in32) :: ii,jj
integer(in32) :: nargs,pathpos
character(len=:), allocatable :: argcurr,pathtemp

!Check and process supplied command arguments
nargs = command_argument_count()
if (nargs == 0) then !Do nothing

else !Scan for and select options 

    !Scan for arguments 
    do ii=1,nargs-1

        !Current argument 
        argcurr = get_command_argument_n_str(ii)

        !If option tag
        if (argcurr == '-o') then !next is options filename with path 

            !Read options filepath and extract optpath
            options%optpath = get_command_argument_n_str(ii+1)
            pathpos = 0 
            do jj=len_trim(options%optpath),1,-1
                if ((options%optpath(jj:jj) == '/') .OR. (options%optpath(jj:jj) == '\')) then 
                    pathpos = jj 
                    exit 
                end if 
            end do 
            if (pathpos .NE. 0) then !path
                pathtemp = options%optpath
                options%opt_filename = pathtemp(pathpos+1:len_trim(pathtemp))
                options%optpath = pathtemp(1:pathpos)
            else !no path 
                options%opt_filename = options%optpath
                options%optpath = ''
            end if 
            ! print *, pathpos
            ! print *, options%optpath
            ! print *, options%opt_filename
        elseif (argcurr == '-m') then !next is mesh filename with path

            !Read mesh filepath and extract iopath 
            options%iopath = get_command_argument_n_str(ii+1)
            pathpos = 0
            do jj=len_trim(options%iopath),1,-1
                if ((options%iopath(jj:jj) == '/') .OR. (options%iopath(jj:jj) == '\')) then 
                    pathpos = jj 
                    exit 
                end if 
            end do 
            if (pathpos .NE. 0) then !path
                pathtemp = options%iopath
                options%mesh_filename = pathtemp(pathpos+1:len_trim(pathtemp))
                options%iopath = pathtemp(1:pathpos)
            else !no path 
                options%mesh_filename = options%iopath
                options%iopath = ''
            end if 
            ! print *, pathpos
            ! print *, options%iopath
            ! print *, options%mesh_filename
        end if 
    end do 
end if 
return 
end subroutine get_process_arguments



!Import mesh subroutine =========================
subroutine flow_import_mesh(mesh,options)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(options_data) :: options 

!Variables - Local 
integer(in) :: ii,vcurr 
integer(in) :: ncell,nedge,nvtx
integer(in32) :: fh

!Display
if (options%csdisp) then
    write(*,'(A)') '--> importing mesh'
end if

!Check if file exists
if (.NOT.file_exists(options%iopath//options%mesh_filename)) then 
    write(*,'(A)') '    ** cannot locate mesh file: '//options%iopath//options%mesh_filename
    options%status = 1
    stop
end if 

!Import 
open(newunit=fh,file=options%iopath//options%mesh_filename)
read(fh,*) ncell,nedge,nvtx
call allocate_mesh(mesh,ncell,nedge,nvtx)
do ii=1,nedge
    read(fh,*) mesh%edge_1(ii),mesh%edge_2(ii),mesh%cell_left(ii),mesh%cell_right(ii)
end do
do ii=1,nvtx
    read(fh,*) vcurr,mesh%vtx_x(ii),mesh%vtx_y(ii)
end do
close(fh)
if (options%csdisp) then
    write(*,'(A,I0)') '    cells: ',ncell
    write(*,'(A,I0)') '    edges: ',nedge
end if 
return 
end subroutine flow_import_mesh




!Set default paths subroutine ===========================
subroutine set_default_paths(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Set paths 
options%optpath = ''
options%iopath = ''
options%mesh_filename = 'grid'
options%opt_filename = 'FLOW_options'
return 
end subroutine set_default_paths




!Set default options subroutine =========================
subroutine set_default_options(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Flow conditions  
options%aoa = 0.0d0 
options%machinf = 0.5d0 
options%rhoinf = 1.225d0 
options%tinf = 288.0d0 
options%gamma = 1.403d0 
options%R = 287.058

!Numerical properties 
options%ittlim_min = 100
options%ittlim = 20000
options%term_res = -8.0d0 
options%cfl = 1.5d0 
options%k2 = 1.0d0 
options%k4 = 0.05d0 
options%c2 = 1.0d0 
options%sp = 0.0d0 
options%nRKstage = 4
allocate(options%RKstagediss_eval(4))
options%RKstagediss_eval(1) = 1 
options%RKstagediss_eval(2) = 0
options%RKstagediss_eval(3) = 0
options%RKstagediss_eval(4) = 0
allocate(options%RKstagediss_apply(4))
options%RKstagediss_apply(1) = 1
options%RKstagediss_apply(2) = 1
options%RKstagediss_apply(3) = 1
options%RKstagediss_apply(4) = 1
options%psNsmooth = 0
options%ts_mode = 'local'
options%damptsteps = .false.

!Boundary condition properties 
options%ff_bc_type = 'riemann'
options%wall_bc_type = 'constant'
options%mflux_bc_type = 'pressure'
options%Mflux_relaxiter = 500
options%Mflux_relax = 0.95d0 
options%Mflux_toll = 0.005d0 
options%backp_pratio = 0.8d0 
options%force_fixed_pratio = .false.
options%force_fixed_ff_ifof_state = .false.

!General 
options%init_state = 'freestream'
options%num_threads = 2
options%mpartition_type = 'prin_axis'
options%avstencil_size = 50 
options%csdisp = .true.
options%pptoggle = .true.
options%ffexporttoggle = .false.
options%ffexport_scale = 'normalised'

!Additional 
options%evalMassflux = .false.
options%massflux_niterav = 2000
options%export_mesh_partitions = .false.
return 
end subroutine set_default_options



!Import options subroutine =========================
subroutine flow_import_options(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Variables - Local
integer(in32) :: fh

!Check if file exists
if (.NOT.file_exists(options%optpath//options%opt_filename)) then 
    write(*,'(A)') '    ** cannot locate options file: '//options%optpath//options%opt_filename
    options%status = 1
    stop
end if 

!Open file 
open(newunit=fh,file=options%optpath//options%opt_filename)

!Flow conditions  
call set_real_opt(options%aoa,fh,'aoa')
call set_real_opt(options%machinf,fh,'mach')
call set_real_opt(options%rhoinf,fh,'rho')
call set_real_opt(options%tinf,fh,'temp')
call set_real_opt(options%gamma,fh,'gamma')
call set_real_opt(options%R,fh,'R')

!Numerical properties 
call set_int_opt(options%ittlim_min,fh,'iterlim_min')
call set_int_opt(options%ittlim,fh,'iterlim')
call set_real_opt(options%term_res,fh,'term_res')
call set_real_opt(options%cfl,fh,'cfl')
call set_real_opt(options%k2,fh,'k2')
call set_real_opt(options%k4,fh,'k4')
call set_real_opt(options%c2,fh,'c2')
call set_real_opt(options%sp,fh,'sp')
call set_int_opt(options%nRKstage,fh,'nRKstage')
call set_int_opt_arr1(options%RKstagediss_eval,fh,'RKstagediss_eval')
call set_int_opt_arr1(options%RKstagediss_apply,fh,'RKstagediss_apply')
call set_int_opt(options%psNsmooth,fh,'npssmooth')
call set_str_opt(options%ts_mode,fh,'ts_mode')
call set_log_opt(options%damptsteps,fh,'damp_ts')

!Boundary condition properties 
call set_str_opt(options%ff_bc_type,fh,'ff_bc_type')
if (options%ff_bc_type == 'riemann') then 
    options%ff_bc_type = 'r'
elseif (options%ff_bc_type == 'characteristic') then 
    options%ff_bc_type = 'c'
end if 
call set_str_opt(options%wall_bc_type,fh,'wall_bc_type')
if (options%wall_bc_type == 'constant') then 
    options%wall_bc_type = 'c'
elseif (options%wall_bc_type == 'linear') then 
    options%wall_bc_type = 'l'
end if 
call set_str_opt(options%mflux_bc_type,fh,'mflux_bc_type')
call set_int_opt(options%Mflux_relaxiter,fh,'mflux_bc_relaxiter')
call set_real_opt(options%Mflux_relax,fh,'mflux_bc_relaxpar')
call set_real_opt(options%Mflux_toll,fh,'mflux_bc_toll')
call set_real_opt(options%backp_pratio,fh,'backp_pratio')
call set_log_opt(options%force_fixed_pratio,fh,'force_fixed_pratio')
call set_log_opt(options%force_fixed_ff_ifof_state,fh,'force_fixed_ff_ifof_state')

!General 
call set_str_opt(options%init_state,fh,'init_state')
call set_int_opt(options%num_threads,fh,'num_threads')
call set_str_opt(options%mpartition_type,fh,'mpartition_type')
call set_int_opt(options%avstencil_size,fh,'avstencil_size')
call set_log_opt(options%csdisp,fh,'cdisplay')
call set_log_opt(options%pptoggle,fh,'postprocess')
call set_log_opt(options%ffexporttoggle,fh,'export_flowfield')
call set_str_opt(options%ffexport_scale,fh,'flowfield_export_scale')

!Additional 
call set_log_opt(options%evalMassflux,fh,'evalMassflux')
call set_int_opt(options%massflux_niterav,fh,'massflux_niterav')
call set_log_opt(options%export_mesh_partitions,fh,'export_mesh_partitions')

!Close file 
close(11)

!Set error status
options%status = 0
return 
end subroutine flow_import_options




!Export force coeficients =========================
subroutine export_force_coefs(filename,options,fcoefs)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(coeficient_data) :: fcoefs
type(options_data) :: options 

!Set error state
if (options%status .NE. 0) then 
    fcoefs%cl_av = 0.0d0 
    fcoefs%cd_av = 0.0d0 
    fcoefs%cm_av = 0.0d0 
    fcoefs%mflux_in = 0.0d0 
end if 

!Export
open(11,file=options%iopath//filename,status='unknown')
    write(11,'(A,A)') 'cl = ',real2F0_Xstring(fcoefs%cl_av,12_in)
    write(11,'(A,A)') 'cd = ',real2F0_Xstring(fcoefs%cd_av,12_in)
    write(11,'(A,A)') 'cm = ',real2F0_Xstring(fcoefs%cm_av,12_in)
    write(11,'(A,A)') 'mflux = ',real2F0_Xstring(fcoefs%mflux_in_av,12_in)
    write(11,'(A,I0)') 'status = ',options%status
close(11)
return 
end subroutine export_force_coefs




!Export TECPLOT solution file subroutine ========================= (from edge2d by Tom Rendall)
subroutine saveVolPLT(filename,mesh,options,flowvars)
use ieee_arithmetic
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars
type(options_data) :: options 

!Variables - Local 
integer(in) :: fh,i,nperline

!Set to zero any nan values of cp
do i=1,mesh%ncell
    if (isnan(flowvars%cp(i))) then 
        flowvars%cp(i) = 0.0d0 
    elseif (.NOT. IEEE_IS_FINITE(flowvars%cp(i))) then
        flowvars%cp(i) = -0.0d0 
    end if 
end do 

!Open file
fh = 11
open(fh,file=options%iopath//filename,status='unknown')

!TECPLOT formatting
write(fh,'(A)',advance="no") 'VARIABLES="X" "Y"'
write(fh,*) '"rho" "u" "v" "vel" "E" "p" "mach" "cp" "vorticity" "resid"'
write(fh,*) 'ZONE T="FlowField"'
write(fh,'(A)',advance="no") 'VARLOCATION=([1,2]=NODAL'
write(fh,*) ',[3,4,5,6,7,8,9,10,11,12]=CELLCENTERED)'
write(fh,*) 'ZONETYPE=FEPOLYGON'
write(fh,'(A,I8)') ' Nodes=',mesh%nvtx
write(fh,'(A,I8)') ' Elements=',mesh%ncell
write(fh,'(A,I8)') ' Faces=',mesh%nedge
write(fh,*) 'NumConnectedBoundaryFaces=0 '
write(fh,*) 'TotalNumBoundaryConnections=0 '

! These loops are because tecplot has a maximum number of characters per line
nperline = 100
write(fh,*) ( mesh%vtx_x(i:min(i+nperline-1,mesh%nvtx)),NEW_LINE('A') , i=1,mesh%nvtx,nperline )
write(fh,*) ( mesh%vtx_y(i:min(i+nperline-1,mesh%nvtx)),NEW_LINE('A') , i=1,mesh%nvtx,nperline )
write(fh,*) ( flowvars%rho(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%u(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%v(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%vel(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%e(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%p(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%mach(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%cp(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%vorticity(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
write(fh,*) ( flowvars%resid(i:min(i+nperline-1,mesh%ncell)),NEW_LINE('A') , i=1,mesh%ncell,nperline )
do i=1,mesh%nedge
    write(fh,*) mesh%edge_1(i),mesh%edge_2(i)  
end do
write(fh,*) ( max(0,mesh%cell_left(i:min(i+nperline-1,mesh%nedge))),NEW_LINE('A') , i=1,mesh%nedge,nperline )
write(fh,*) ( max(0,mesh%cell_right(i:min(i+nperline-1,mesh%nedge))),NEW_LINE('A') , i=1,mesh%nedge,nperline )


!Close file
close(fh)
return 
end subroutine saveVolPLT




!Export TECPLOT mesh file with cell data subroutine =========================  
subroutine write_cell_dataPLT(filename,volume_mesh,cell_datap)
use ieee_arithmetic
implicit none 


!Variables - Import
real(dp), dimension(:) :: cell_datap
character(*), intent(in) :: filename
type(mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: fh,i,nperline

!Open file
fh = 11
open(fh,file=filename//'.plt',status='unknown')

!TECPLOT formatting
write(fh,'(A)',advance="no") 'VARIABLES="X" "Y" "data"'
write(fh,*) 'ZONE T="CellData"'
write(fh,'(A)',advance="no") 'VARLOCATION=([1,2]=NODAL,[3]=CELLCENTERED)'
write(fh,*) 'ZONETYPE=FEPOLYGON'
write(fh,'(A,I8)') ' Nodes=',volume_mesh%nvtx
write(fh,'(A,I8)') ' Elements=',volume_mesh%ncell
write(fh,'(A,I8)') ' Faces=',volume_mesh%nedge
write(fh,*) 'NumConnectedBoundaryFaces=0 '
write(fh,*) 'TotalNumBoundaryConnections=0 '

! These loops are because tecplot has a maximum number of characters per line
nperline = 100
write(fh,*) ( volume_mesh%vtx_x(i:min(i+nperline-1,volume_mesh%nvtx)),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( volume_mesh%vtx_y(i:min(i+nperline-1,volume_mesh%nvtx)),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( cell_datap(i:min(i+nperline-1,volume_mesh%ncell)),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
do i=1,volume_mesh%nedge
    write(fh,*) volume_mesh%edge_1(i),volume_mesh%edge_2(i)  
end do
write(fh,*) ( max(0,volume_mesh%cell_left(i:min(i+nperline-1,volume_mesh%nedge))),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )
write(fh,*) ( max(0,volume_mesh%cell_right(i:min(i+nperline-1,volume_mesh%nedge))),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )

!Close file
close(fh)
return 
end subroutine write_cell_dataPLT




!Export flow field data surboutine  =========================
subroutine export_flowfield(filename,flowvars,mesh,options)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(flow_var_data) :: flowvars
type(options_data) :: options 
type(mesh_data) :: mesh

!Variables - Local 
integer(in) :: cc 

!Export
open(11,file=options%iopath//filename,status='unknown')
write(11,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)') '           ','u','                ','v','               ','mach','              ',&
                            'p','               ','rho','               ','cp','               ','T','            ','vorticity'
    do cc=1,mesh%ncell 
        write(11,'(A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8)') ' ',flowvars%u(cc),' ',flowvars%v(cc),' ',&
    flowvars%mach(cc),' ',flowvars%p(cc),' ',flowvars%rho(cc),' ',flowvars%cp(cc),' ',flowvars%temp(cc),' ',flowvars%vorticity(cc)
    end do 
close(11)
return 
end subroutine export_flowfield




!Write restart file =========================
subroutine write_restart_file(filename,flowvars,mesh,options)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(flow_var_data) :: flowvars
type(options_data) :: options 
type(mesh_data) :: mesh

!Variables - Local 
integer(in) :: cc 

!Write
open(11,file=options%iopath//filename,status='unknown')

!Write general data 
write(11,'(A,f16.8,A,f16.8)') ' ',flowvars%mflux_pset,' ',flowvars%vnset_of

!Write flowfield data
do cc=1,mesh%ncell
    write(11,'(A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8,A,f16.8)') ' ',flowvars%u(cc),' ',flowvars%v(cc),' ',&
    flowvars%mach(cc),' ',flowvars%p(cc),' ',flowvars%rho(cc),' ',flowvars%cp(cc),' ',flowvars%e(cc)
end do 
close(11)
return 
end subroutine write_restart_file




!Read restart file =========================
subroutine read_restart_file(filename,flowvars,mesh,options)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(flow_var_data) :: flowvars
type(options_data) :: options 
type(mesh_data) :: mesh

!Variables - Local 
integer(in) :: cc 

!Read
open(11,file=options%iopath//filename,status='unknown')

!Write general data 
read(11,*) flowvars%mflux_pset,flowvars%vnset_of

!Write flowfield data
do cc=1,mesh%ncell
    read(11,*) flowvars%u(cc),flowvars%v(cc),&
    flowvars%mach(cc),flowvars%p(cc),flowvars%rho(cc),flowvars%cp(cc),flowvars%e(cc)
end do 
close(11)
return 
end subroutine read_restart_file
    



!Write status subroutine ===========================
subroutine export_status(options,flowvars,flowvarsAV)
implicit none 

!Variables - Import
type(options_data) :: options
type(flow_var_data) :: flowvars,flowvarsAV

!Write
open(11,file=options%iopath//'flow2d_status') 
    write(11,'(I0)') options%status 
    write(11,'(A)') real2F0_Xstring(flowvars%rho_res,12_in)
    write(11,'(A)') real2F0_Xstring(flowvarsAV%rho_res,12_in) 
close(11)
return 
end subroutine export_status


end module flow_io_mod