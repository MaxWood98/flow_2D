!Flow solver io module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.5.9
!Updated: 23/08/23

!Module
module flow_io_mod
! use flow_data_mod
use flow_general_mod
contains 


!Read command arguments subroutine ===========================
subroutine get_process_arguments(options)
implicit none

!Variables - Import
type(options_data) :: options
    
!Variables - Local 
integer(in) :: nargs
integer(in32) :: arglen,argstat

!Check and process supplied command arguments
nargs = command_argument_count()
if (nargs == 0) then !Use default path
    allocate(character(len=3) :: options%optpath)
    options%optpath = 'io/'
    allocate(character(len=3) :: options%iopath)
    options%iopath = 'io/'
else !Use specified path 
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: options%optpath)
    call get_command_argument(number=1, value=options%optpath, status=argstat)
    call get_command_argument(number=2, length=arglen)
    allocate(character(len=arglen) :: options%iopath)
    call get_command_argument(number=2, value=options%iopath, status=argstat)
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

!Display
if (options%csdisp) then
    write(*,'(A)') '--> importing mesh'
end if

!Import 
open(11,file=options%iopath//'grid')
read(11,*) ncell,nedge,nvtx
call allocate_mesh(mesh,ncell,nedge,nvtx)
do ii=1,nedge
    read(11,*) mesh%edge_1(ii),mesh%edge_2(ii),mesh%cell_left(ii),mesh%cell_right(ii)
end do
do ii=1,nvtx
    read(11,*) vcurr,mesh%vtx_x(ii),mesh%vtx_y(ii)
end do
close(11)
if (options%csdisp) then
    write(*,*) '   cells: ',ncell
    write(*,*) '   edges: ',nedge
end if 
return 
end subroutine flow_import_mesh




!Import options subroutine =========================
subroutine flow_import_options(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Import 
open(11,file=options%optpath//'FLOW_options.dat')
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) options%aoa
read(11,*) !skip
read(11,*) !skip
read(11,*) options%machinf
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) options%rhoinf
read(11,*) !skip
read(11,*) !skip
read(11,*) options%tinf
read(11,*) !skip
read(11,*) !skip
read(11,*) options%gamma
read(11,*) !skip
read(11,*) !skip
read(11,*) options%R 
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) options%ittlim_min
read(11,*) !skip
read(11,*) !skip
read(11,*) options%ittlim
read(11,*) !skip
read(11,*) !skip
read(11,*) options%term_res
read(11,*) !skip
read(11,*) !skip
read(11,*) options%cfl
read(11,*) !skip
read(11,*) !skip
read(11,*) options%k2
read(11,*) options%k4
read(11,*) options%c2
read(11,*) options%sp
read(11,*) !skip
read(11,*) !skip
read(11,*) options%nRKstage
read(11,*) !skip
read(11,*) !skip
read(11,*) options%RKstagediss(1)
read(11,*) options%RKstagediss(2)
read(11,*) options%RKstagediss(3)
read(11,*) options%RKstagediss(4)
read(11,*) !skip
read(11,*) !skip
read(11,*) options%psNsmooth
read(11,*) !skip
read(11,*) !skip
read(11,*) options%ts_mode
read(11,*) !skip
read(11,*) !skip
read(11,*) options%damptsteps
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) options%ff_bc_type
read(11,*) !skip
read(11,*) !skip
read(11,*) options%wall_bc_type
read(11,*) !skip
read(11,*) !skip
read(11,*) options%mflux_bc_type
read(11,*) !skip
read(11,*) !skip
read(11,*) options%Mflux_relaxiter
read(11,*) options%Mflux_relax
read(11,*) options%Mflux_toll
read(11,*) !skip
read(11,*) !skip
read(11,*) options%backp_pratio
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) options%init_state
read(11,*) !skip
read(11,*) !skip
read(11,*) options%num_threads
read(11,*) !skip
read(11,*) !skip
read(11,*) options%avstencil_size
read(11,*) !skip
read(11,*) !skip
read(11,*) options%csdisp
read(11,*) !skip
read(11,*) !skip
read(11,*) options%pptoggle
read(11,*) !skip
read(11,*) !skip
read(11,*) options%ffexporttoggle
read(11,*) !skip
read(11,*) !skip
read(11,*) options%ffexport_scale
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) options%evalMassflux
read(11,*) options%massflux_niterav
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
    write(11,'(A,A)') 'mflux = ',real2F0_Xstring(fcoefs%mflux_in,12_in)
    write(11,'(A,I0)') 'status = ',options%status
close(11)
return 
end subroutine export_force_coefs




!F0.X format with leading zero function =========================
function real2F0_Xstring(val,X) result(str)

!Result 
character(len=:), allocatable :: str

!Variables - Import 
character(len=10) :: frmtI
character(len=:), allocatable :: frmt,str_I
integer(in) :: X,len_frmt,len_str
real(dp) :: val

!Set format descriptor
write(frmtI,'(I0)') X
len_frmt = len_trim(frmtI)
allocate(character(len=len_frmt) :: frmt)
frmt = frmtI(1:len_frmt)
frmt = frmt//')'
frmt = '(F0.'//frmt

!Allocate initial character
allocate(character(len=2*X) :: str_I)

!Write data to return charachter
write(str_I,frmt) val

!Allocate return character
len_str = len_trim(str_I)
allocate(character(len=len_str) :: str)
str = str_I(1:len_str)

!Assign leading zero if required
if (str(1:1) == '.') then 
    str = '0'//str
elseif (str(1:2) == '-.') then 
    str = '-0.'//str(3:len_str)
end if 
end function real2F0_Xstring




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