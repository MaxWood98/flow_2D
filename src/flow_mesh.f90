!Flow solver mesh subroutines module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.5.8
!Updated: 26/03/24

!Module
module flow_mesh_mod
use flow_io_mod
contains 

!Mesh preprocessing subroutine ========================= 
subroutine flow_preprocess_mesh(mesh,options,invalid_mesh,dispt)
implicit none 

!Variables - Import
integer(in) :: invalid_mesh,dispt
type(mesh_data) :: mesh
type(options_data) :: options 

!Variables - Local
integer(in) :: ii,neonw,nconw
integer(in) :: cell_on_wall(mesh%ncell),edge_on_wallc(mesh%nedge)
real(dp) :: vol_edge,minvol,maxvol,emx,emy

!Display
if ((options%csdisp) .AND. (dispt == 1)) then
    write(*,'(A)') '--> preprocessing mesh'
end if

!Calculate cell properties  
mesh%cell_elenint(:) = 0.0d0  
mesh%cell_elentot(:) = 0.0d0  
mesh%cellcenx(:) = 0.0d0 
mesh%cellceny(:) = 0.0d0 
mesh%cell_nedge(:) = 0
mesh%cell_nedgei(:) = 0
mesh%cell_vol(:) = 0.0d0 
invalid_mesh = 0
do ii=1,mesh%nedge

    !Check for invalid mesh configuration
    if (mesh%cell_right(ii) .LE. 0) then
        invalid_mesh = 1
        options%status = invalid_mesh
        write(*,*) '** warning -> invalid mesh configuration (cell right <= 0)'
        return 
    end if

    !Edge volume contribution 
    vol_edge = 0.5d0*(mesh%vtx_x(mesh%edge_1(ii))*mesh%vtx_y(mesh%edge_2(ii)) - &
                      mesh%vtx_x(mesh%edge_2(ii))*mesh%vtx_y(mesh%edge_1(ii)))
    if (mesh%cell_left(ii) .GT. 0) then
        mesh%cell_vol(mesh%cell_left(ii)) = mesh%cell_vol(mesh%cell_left(ii)) + vol_edge
    end if
    if (mesh%cell_right(ii) .GT. 0) then
        mesh%cell_vol(mesh%cell_right(ii)) = mesh%cell_vol(mesh%cell_right(ii)) - vol_edge
    end if

    !Edge dx and dy
    mesh%edge_dx(ii) = (mesh%vtx_x(mesh%edge_2(ii)) - mesh%vtx_x(mesh%edge_1(ii)))
    mesh%edge_dy(ii) = (mesh%vtx_y(mesh%edge_2(ii)) - mesh%vtx_y(mesh%edge_1(ii)))

    !Edge length 
    mesh%edgelen(ii) = sqrt(mesh%edge_dx(ii)**2 + mesh%edge_dy(ii)**2)

    !Edge midpoint 
    mesh%edge_mx(ii) = 0.5d0*(mesh%vtx_x(mesh%edge_1(ii)) + mesh%vtx_x(mesh%edge_2(ii)))
    mesh%edge_my(ii) = 0.5d0*(mesh%vtx_y(mesh%edge_1(ii)) + mesh%vtx_y(mesh%edge_2(ii)))

    !Accumulate cell centres
    emx = mesh%edge_mx(ii) 
    emy = mesh%edge_my(ii) 
    mesh%cellcenx(mesh%cell_right(ii)) = mesh%cellcenx(mesh%cell_right(ii)) + emx*mesh%edgelen(ii)
    mesh%cellceny(mesh%cell_right(ii)) = mesh%cellceny(mesh%cell_right(ii)) + emy*mesh%edgelen(ii)
    mesh%cell_elentot(mesh%cell_right(ii)) = mesh%cell_elentot(mesh%cell_right(ii)) + mesh%edgelen(ii)
    if (mesh%cell_left(ii) .GT. 0) then
        mesh%cellcenx(mesh%cell_left(ii)) = mesh%cellcenx(mesh%cell_left(ii)) + emx*mesh%edgelen(ii)
        mesh%cellceny(mesh%cell_left(ii)) = mesh%cellceny(mesh%cell_left(ii)) + emy*mesh%edgelen(ii)
        mesh%cell_elentot(mesh%cell_left(ii)) = mesh%cell_elentot(mesh%cell_left(ii)) + mesh%edgelen(ii)
    end if 
    
    !Edge unit normal vector 
    if (mesh%edgelen(ii) == 0.0d0) then 
        write(*,'(A,I8)') '** warning -> zero length mesh edge: ',ii
        mesh%edge_nx(ii) = 0.0d0 
        mesh%edge_ny(ii) = 0.0d0 
    else
        mesh%edge_nx(ii) = mesh%edge_dy(ii)/mesh%edgelen(ii)
        mesh%edge_ny(ii) = -mesh%edge_dx(ii)/mesh%edgelen(ii)
    end if

    !Accumulate number of edges on each cell and internal edge lengths
    if (mesh%cell_left(ii) .GT. 0) then
        mesh%cell_nedge(mesh%cell_left(ii)) = mesh%cell_nedge(mesh%cell_left(ii)) + 1
    end if
    if (mesh%cell_right(ii) .GT. 0) then
        mesh%cell_nedge(mesh%cell_right(ii)) = mesh%cell_nedge(mesh%cell_right(ii)) + 1
    end if
    ! if ((mesh%cell_left(ii) .GT. 0) .AND. (mesh%cell_right(ii) .GT. 0)) then
    !     mesh%cell_nedgei(mesh%cell_left(ii)) = mesh%cell_nedgei(mesh%cell_left(ii)) + 1
    !     mesh%cell_nedgei(mesh%cell_right(ii)) = mesh%cell_nedgei(mesh%cell_right(ii)) + 1
    !     mesh%cell_elenint(mesh%cell_left(ii)) = mesh%cell_elenint(mesh%cell_left(ii)) + mesh%edgelen(ii)
    !     mesh%cell_elenint(mesh%cell_right(ii)) = mesh%cell_elenint(mesh%cell_right(ii)) + mesh%edgelen(ii)
    ! end if 
    if ((mesh%cell_left(ii) .NE. -1) .AND. (mesh%cell_right(ii) .NE. -1)) then
        mesh%cell_nedgei(mesh%cell_left(ii)) = mesh%cell_nedgei(mesh%cell_left(ii)) + 1
        mesh%cell_nedgei(mesh%cell_right(ii)) = mesh%cell_nedgei(mesh%cell_right(ii)) + 1
        mesh%cell_elenint(mesh%cell_left(ii)) = mesh%cell_elenint(mesh%cell_left(ii)) + mesh%edgelen(ii)
        mesh%cell_elenint(mesh%cell_right(ii)) = mesh%cell_elenint(mesh%cell_right(ii)) + mesh%edgelen(ii)
    end if 
end do
mesh%cellcenx(:) = mesh%cellcenx(:)/mesh%cell_elentot(:)
mesh%cellceny(:) = mesh%cellceny(:)/mesh%cell_elentot(:)
minvol = minval(mesh%cell_vol(:))
maxvol = maxval(mesh%cell_vol(:))
mesh%cvol_tot = sum(mesh%cell_vol(:))
if ((options%csdisp) .AND. (dispt == 1)) then
    write(*,'(A,E12.5,A,E12.5)') '    cell volume (max/min): ',maxvol, ' / ',minvol
end if 

!Identify active boundary conditions 
mesh%bc_active(:) = 0 
do ii=1,mesh%nedge
    if (mesh%cell_left(ii) .LT. 0) then
        mesh%bc_active(abs(mesh%cell_left(ii))) = 1
    end if
end do 

!Initialise boundary condition states
mesh%bc_state(:) = 0 

!If pressure outflow boundary condition is active then tag and identify outflow zones 
if (mesh%bc_active(7) == 1) then 
    call identify_outflow_zones(mesh,invalid_mesh)
else
    mesh%nzone_outflow = 0
    mesh%outflow_zone(:) = 0
    mesh%nedge_of_zone(:) = 0
end if 
if (mesh%nzone_outflow .NE. 0) then 
    if ((options%csdisp) .AND. (dispt == 1)) then
        write(*,'(A,I0,A)') '    identifed ',mesh%nzone_outflow,' outflow zones'
    end if 
end if 

!Calculate global chord
call calculate_chord(mesh,invalid_mesh)
if ((options%csdisp) .AND. (dispt == 1)) then
    write(*,'(A,A)') '    chord (x): ',real2F0_Xstring(mesh%chordx,6_in) 
end if 

!Build list of cells adjacent to walls and the edges on these cells 
cell_on_wall(:) = 0 
do ii=1,mesh%nedge
    if (mesh%cell_left(ii) == -1) then
        cell_on_wall(mesh%cell_right(ii)) = 1
    end if
end do 
edge_on_wallc(:) = 0 
do ii=1,mesh%nedge
    if (cell_on_wall(mesh%cell_right(ii)) == 1) then 
        edge_on_wallc(ii) = 1
    end if 
    if (mesh%cell_left(ii) .GT. 0) then
        if (cell_on_wall(mesh%cell_left(ii)) == 1) then 
            edge_on_wallc(ii) = 1
        end if 
    end if 
end do 
mesh%npwallcell = sum(cell_on_wall(:))

!Set edge counts for damping calculations 
mesh%cell_nedgeiw(:) = mesh%cell_nedgei(:)
do ii=1,mesh%ncell 
    if (cell_on_wall(ii) == 1) then 
        mesh%cell_nedgeiw(ii) = mesh%cell_nedgeiw(ii) + 1
    end if 
end do 

!Find adjacent cells to cells on walls 
do ii=1,mesh%nedge
    if (cell_on_wall(mesh%cell_right(ii)) == 1) then 
        if (mesh%cell_left(ii) .GT. 0) then
            if (cell_on_wall(mesh%cell_left(ii)) == 0) then 
                cell_on_wall(mesh%cell_left(ii)) = 2
            end if 
        end if 
    end if 
    if (mesh%cell_left(ii) .GT. 0) then
        if (cell_on_wall(mesh%cell_left(ii)) == 1) then 
            if (cell_on_wall(mesh%cell_right(ii)) == 0) then 
                cell_on_wall(mesh%cell_right(ii)) = 2
            end if 
        end if 
    end if 
end do 
do ii=1,mesh%ncell
    if (cell_on_wall(ii) == 2) then 
        cell_on_wall(ii) = 1
    end if 
end do 
mesh%npwallcell = sum(cell_on_wall(:))

!Build list of cells for wall pressure gradient evaluation 
allocate(mesh%cellonwall(mesh%npwallcell))
nconw = 0 
do ii=1,mesh%ncell
    if (cell_on_wall(ii) == 1) then 
        nconw = nconw + 1
        mesh%cellonwall(nconw) = ii
    end if 
end do 

!Build list of edges for wall pressure gradient evaluation 
mesh%npwalledge = sum(edge_on_wallc(:))
neonw = 0 
allocate(mesh%edgeonwcell(mesh%npwalledge))
do ii=1,mesh%nedge
    if (edge_on_wallc(ii) == 1) then 
        neonw = neonw + 1
        mesh%edgeonwcell(neonw) = ii
    end if 
end do 

!Set status
options%status = invalid_mesh
return 
end subroutine flow_preprocess_mesh




!Outflow zone identification subroutine =========================
subroutine identify_outflow_zones(mesh,invalid_mesh)
implicit none 

!Variables - Import
integer(in) :: invalid_mesh
type(mesh_data) :: mesh

!Variables - Local
integer(in) :: ii,ee,ff
integer(in) :: Nedge_outflow,etgt,eadj,ebase,v1,v2,Nzone_outflow,nupdate
integer(in) :: Nedge_of_zone(mesh%nedge)
integer(in) :: V2E_of(mesh%nvtx,2)
integer(in), dimension(:), allocatable :: edges_outflow,of_zone

!Build list of outflow edges 
Nedge_outflow = 0 
do ee=1,mesh%nedge
    if (mesh%cell_left(ee) == -7) then
        Nedge_outflow = Nedge_outflow + 1
    end if 
end do 
allocate(edges_outflow(Nedge_outflow))
Nedge_outflow = 0 
do ee=1,mesh%nedge
    if (mesh%cell_left(ee) == -7) then
        Nedge_outflow = Nedge_outflow + 1
        edges_outflow(Nedge_outflow) = ee 
    end if 
end do 

!Build V2E for outflow edges 
V2E_of(:,:) = 0 
do ee=1,Nedge_outflow
    etgt = edges_outflow(ee)
    v1 = mesh%edge_1(etgt)
    v2 = mesh%edge_2(etgt)
    if (V2E_of(v1,1) == 0) then 
        V2E_of(v1,1) = etgt 
    elseif (V2E_of(v1,2) == 0) then 
        V2E_of(v1,2) = etgt 
    else 
        print *, '** greater than valence 2 outflow vertex: ',v1
        invalid_mesh = 1
        return 
    end if 
    if (V2E_of(v2,1) == 0) then 
        V2E_of(v2,1) = etgt 
    elseif (V2E_of(v2,2) == 0) then 
        V2E_of(v2,2) = etgt 
    else 
        print *, '** greater than valence 2 outflow vertex: ',v2
        invalid_mesh = 1
        return 
    end if 
end do 

!Flood outflow zones 
Nzone_outflow = 0 
allocate(of_zone(mesh%nedge))
of_zone(:) = 0 
Nedge_of_zone(:) = 0 
do ff=1,Nedge_outflow

    !Find zero tagged edge 
    ebase = 0 
    do ee=1,Nedge_outflow
        etgt = edges_outflow(ee)
        if (of_zone(etgt) == 0) then 
            ebase = etgt 
            exit 
        end if 
    end do  

    !Exit at complete flood of all zones 
    if (ebase == 0) then 
        exit 
    end if 

    !Increment zone count 
    Nzone_outflow = Nzone_outflow + 1

    !Flood this outflow zone 
    of_zone(ebase) = Nzone_outflow
    Nedge_of_zone(Nzone_outflow) = 1
    do ii=1,2*Nedge_outflow
        nupdate = 0 
        do ee=1,Nedge_outflow
            etgt = edges_outflow(ee)
            v1 = mesh%edge_1(etgt)
            v2 = mesh%edge_2(etgt)
            if (of_zone(etgt) .NE. 0) then !flood to adjacent edges 
                
                !Flood across v1
                eadj = 0 
                if (V2E_of(v1,1) == etgt) then 
                    eadj = V2E_of(v1,2) 
                elseif (V2E_of(v1,2) == etgt) then 
                    eadj = V2E_of(v1,1) 
                end if 
                if (eadj .NE. 0) then 
                    if (of_zone(eadj) == 0) then 
                        of_zone(eadj) = of_zone(etgt)
                        Nedge_of_zone(Nzone_outflow) = Nedge_of_zone(Nzone_outflow) + 1
                        nupdate = nupdate + 1
                    end if 
                end if 
                
                !Flood across v2
                eadj = 0 
                if (V2E_of(v2,1) == etgt) then 
                    eadj = V2E_of(v2,2) 
                elseif (V2E_of(v2,2) == etgt) then 
                    eadj = V2E_of(v2,1) 
                end if 
                if (eadj .NE. 0) then 
                    if (of_zone(eadj) == 0) then 
                        of_zone(eadj) = of_zone(etgt)
                        Nedge_of_zone(Nzone_outflow) = Nedge_of_zone(Nzone_outflow) + 1
                        nupdate = nupdate + 1
                    end if 
                end if
            end if 
        end do
        if (nupdate == 0) then 
            exit 
        end if 
    end do 
end do 

!Tag edges with their outflow zone 
mesh%nzone_outflow = Nzone_outflow
mesh%outflow_zone(:) = of_zone(:)
mesh%nedge_of_zone(:) = Nedge_of_zone(:)
return 
end subroutine identify_outflow_zones




!Chord calculation subroutine ========================= 
subroutine calculate_chord(mesh,invalid_mesh)
implicit none 

!Variables - Import
integer(in) :: invalid_mesh
type(mesh_data) :: mesh

!Variables - Local
integer(in) :: ii,vtgt
real(dp) :: xmax,xmin 

!Determine xmax and xmin surface points
xmax = minval(mesh%vtx_x(:))
xmin = maxval(mesh%vtx_x(:))
do ii=1,mesh%nedge
    if (mesh%cell_left(ii) == -1) then !Wall segment
        vtgt = mesh%edge_1(ii)
        if (mesh%vtx_x(vtgt) .GT. xmax) then 
            xmax = mesh%vtx_x(vtgt)
        end if
        if (mesh%vtx_x(vtgt) .LT. xmin) then 
            xmin = mesh%vtx_x(vtgt)
        end if
        vtgt = mesh%edge_2(ii)
        if (mesh%vtx_x(vtgt) .GT. xmax) then 
            xmax = mesh%vtx_x(vtgt)
        end if
        if (mesh%vtx_x(vtgt) .LT. xmin) then 
            xmin = mesh%vtx_x(vtgt)
        end if
    end if
end do 

!Calculate chord 
mesh%chordx = abs(xmax - xmin)

!Check for invalid chord
if (mesh%chordx == 0.0d0) then
    invalid_mesh = 1
    write(*,*) '** warning -> invalid mesh configuration (zero chord)'
end if
return 
end subroutine calculate_chord




!Mesh parallel zone construction subroutine ========================= 
subroutine flow_par_segment_mesh(mesh_par,mesh,options,nanflag)
implicit none 

!Variables - Import
integer(in) :: nanflag
type(mesh_data) :: mesh
type(options_data) :: options
type(mesh_data), dimension(:), allocatable :: mesh_par

!Variables - Local 
integer(in) :: ii,zz,ff,cc,ee,tt,cell_maxne,Nminzone,Nremain,zoneidx,nzonec,nzonec0,NfrontC,NfrontN,Nflooditer,cellC
integer(in) :: maxcellzone,cl,cr,ncell,nedge,nvtx,v1,v2,invalid_mesh,zl,zr,etgt,Npartition
integer(in) :: cellnedge(mesh%ncell),cellzone(mesh%ncell),ncellzone(options%num_threads),nedgezone(options%num_threads)
integer(in) :: frontC(mesh%ncell),frontN(mesh%ncell),edge_zone(mesh%nedge),vtxtag(mesh%nvtx),cell_zone_idx(mesh%ncell)
integer(in) :: cell_mult_zone(mesh%ncell),zone_tag(options%num_threads)
integer(in), dimension(:,:), allocatable :: cell2cell,cell2edge
real(dp) :: ztol,mxmin,dx1,dx2 

!Set zero tollerance
ztol = 0.000001d0 

!Allocate parallel mesh structure 
allocate(mesh_par(options%num_threads))

!Build cell2cell 
cellnedge(:) = 0 
do ii=1,mesh%nedge 
    if (mesh%cell_left(ii) .GT. 0) then
        cellnedge(mesh%cell_left(ii)) = cellnedge(mesh%cell_left(ii)) + 1
    end if
    if (mesh%cell_right(ii) .GT. 0) then
        cellnedge(mesh%cell_right(ii)) = cellnedge(mesh%cell_right(ii)) + 1
    end if
end do 
cell_maxne = maxval(cellnedge(:))
cellnedge(:) = 0 
allocate(cell2cell(mesh%ncell,cell_maxne)) 
allocate(cell2edge(mesh%ncell,cell_maxne))
cell2cell(:,:) = 0 
cell2edge(:,:) = 0 
do ii=1,mesh%nedge 
    if (mesh%cell_left(ii) .GT. 0) then
        cellnedge(mesh%cell_left(ii)) = cellnedge(mesh%cell_left(ii)) + 1
        cell2cell(mesh%cell_left(ii),cellnedge(mesh%cell_left(ii))) = mesh%cell_right(ii)
        cell2edge(mesh%cell_left(ii),cellnedge(mesh%cell_left(ii))) = ii
    end if
    if (mesh%cell_right(ii) .GT. 0) then
        cellnedge(mesh%cell_right(ii)) = cellnedge(mesh%cell_right(ii)) + 1
        cell2cell(mesh%cell_right(ii),cellnedge(mesh%cell_right(ii))) = mesh%cell_left(ii)
        cell2edge(mesh%cell_right(ii),cellnedge(mesh%cell_right(ii))) = ii
    end if
end do  

!Initialise cell zones 
nzonec = 1 
cellzone(:) = 1 

!Divide mesh into paralell zones 
if (options%mpartition_type == 'prin_axis') then 
    Npartition = nint(log(real(options%num_threads,dp))/log(2.0d0))
    do ii=1,Npartition
        nzonec0 = nzonec
        if (mod(ii,2) == 0) then 
            do zz=1,nzonec
                call split_half_zone(cellzone,nzonec,mesh,mesh%cellcenx,mesh%cellceny,zz)
            end do 
        else
            do zz=1,nzonec
                call split_half_zone(cellzone,nzonec,mesh,mesh%cellceny,mesh%cellcenx,zz)
            end do 
        end if 
    end do 
elseif (options%mpartition_type == 'flood') then 

    !Set number of cells in each parallel zone 
    ncellzone(:) = 0 
    Nminzone = floor(real(mesh%ncell,dp)/real(options%num_threads,dp))
    Nremain = mesh%ncell - options%num_threads*Nminzone
    ncellzone(:) = Nminzone
    do ii=1,Nremain
        ncellzone(ii) = ncellzone(ii) + 1
    end do 
    maxcellzone = maxval(ncellzone(:))

    !Initialise zone flooding state
    zoneidx = 1
    nzonec = 0

    !Identify flood base cells with xmin edge at xmin far field of mesh
    mxmin = ieee_value(1.0d0,ieee_positive_inf) !minval(mesh%vtx_x(:))
    do ii=1,mesh%nedge 
        if (mesh%vtx_x(mesh%edge_1(ii)) .LE. mxmin) then 
            mxmin = mesh%vtx_x(mesh%edge_1(ii)) 
        end if 
        if (mesh%vtx_x(mesh%edge_2(ii)) .LE. mxmin) then 
            mxmin = mesh%vtx_x(mesh%edge_2(ii)) 
        end if 
    end do 
    cellzone(:) = 0
    do ii=1,mesh%nedge 
        dx1 = sqrt((mxmin - mesh%vtx_x(mesh%edge_1(ii)))**2)
        dx2 = sqrt((mxmin - mesh%vtx_x(mesh%edge_2(ii)))**2) 
        if ((dx1 .LE. ztol) .OR. (dx2 .LE. ztol)) then 
            if (mesh%cell_left(ii) .GT. 0) then
                if (cellzone(mesh%cell_left(ii)) == 0) then 
                    cellzone(mesh%cell_left(ii)) = zoneidx
                    nzonec = nzonec + 1
                    if (nzonec .GE. ncellzone(zoneidx)) then 
                        nzonec = 0 
                        zoneidx = zoneidx + 1
                    end if
                end if 
            end if
            if (mesh%cell_right(ii) .GT. 0) then
                if (cellzone(mesh%cell_right(ii)) == 0) then 
                    cellzone(mesh%cell_right(ii)) = zoneidx
                    nzonec = nzonec + 1
                    if (nzonec .GE. ncellzone(zoneidx)) then 
                        nzonec = 0 
                        zoneidx = zoneidx + 1
                    end if
                end if 
            end if 
        end if 
    end do 

    !Construct initial flood front 
    NfrontC = 0 
    NfrontN = 0 
    frontC(:) = 0
    frontN(:) = 0
    do ii=1,mesh%ncell
        if (cellzone(ii) .NE. 0) then 
            NfrontC = NfrontC + 1
            frontC(NfrontC) = ii
        end if
    end do 

    !Flood mesh to build parallel zones 
    Nflooditer = 0 
    do ff=1,mesh%ncell !Flood iteration
        NfrontN = 0
        do cc=1,NfrontC !Each cell in front 
            cellC = frontC(cc) !Cell 
            do ee=1,cellnedge(cellC) !Adjacent cells
                if (cell2cell(cellC,ee) .GT. 0) then 
                    if (cellzone(cell2cell(cellC,ee)) == 0) then 

                        !Add to current zone 
                        cellzone(cell2cell(cellC,ee)) = zoneidx
                        nzonec = nzonec + 1

                        !Increment zone
                        if (nzonec .GE. ncellzone(zoneidx)) then 
                            nzonec = 0 
                            zoneidx = zoneidx + 1
                        end if 
                            
                        !Add to new front
                        frontN(NfrontN+1) = cell2cell(cellC,ee)
                        NfrontN = NfrontN + 1
                    end if 
                end if 
            end do 
        end do 

        !Exit at complete flood
        if (NfrontN == 0) then 
            exit 
        end if 

        !Update front
        NfrontC = NfrontN
        frontC(1:NfrontN) = frontN(1:NfrontN)
        Nflooditer = Nflooditer + 1
    end do 

    !Check for completed flood
    if (minval(cellzone) .NE. 0) then 
        if (options%csdisp) then
            write(*,'(A,I4,A)') '    {mesh partitioning complete in ',Nflooditer,' iterations}'
        end if
    else    
        nanflag = 1
        write(*,'(A)') '     *** mesh partitioning error ***'
        return 
    end if 
else
    nanflag = 1
    write(*,'(A)') '--> ** invalid mesh partitioning method selected **'
    stop 
end if 

!Assign edges to the zone of their cell right cell 
nedgezone(:) = 0
edge_zone(:) = 0 
do ee=1,mesh%nedge
    edge_zone(ee) = cellzone(mesh%cell_right(ee))
end do 
do ee=1,mesh%nedge
    nedgezone(edge_zone(ee)) = nedgezone(edge_zone(ee)) + 1
end do 

!Index cells in each zone 
ncellzone(:) = 0 
cell_zone_idx(:) = 0 
do cc=1,mesh%ncell 
    ncellzone(cellzone(cc)) = ncellzone(cellzone(cc)) + 1
    cell_zone_idx(cc) = ncellzone(cellzone(cc)) 
end do 
! do cc=1,options%num_threads
!     print *,'zone = ',cc,' Ncell = ',ncellzone(cc)
! end do 

!Tag each cell that contains edges of multiple zones 
cell_mult_zone(:) = 0 
do cc=1,mesh%ncell 
    zone_tag(:) = 0 
    do ee=1,cell_maxne
        etgt = cell2edge(cc,ee)
        if (etgt .GT. 0) then 
            cl = mesh%cell_left(etgt)
            cr = mesh%cell_right(etgt)
            if (cl .GT. 0) then
                zl = cellzone(cl)
            else
                zl = 0 
            end if 
            if (cr .GT. 0) then
                zr = cellzone(cr)
            else 
                zr = 0 
            end if 
            if (zl .GT. 0) then  
                zone_tag(zl) = 1
            end if 
            if (zr .GT. 0) then  
                zone_tag(zr) = 1
            end if 
        end if 
    end do 
    if (sum(zone_tag(:)) .GT. 1) then 
        cell_mult_zone(cc) = 1
    end if 
end do 
! cell_mult_zone(:) = 1 !force atomic accumulation on all cells

!Construct local mesh in each zone 
invalid_mesh = 0 
do tt=1,options%num_threads

    !Count vertices in this zone
    vtxtag(:) = 0 
    do ii=1,mesh%nedge
        if (edge_zone(ii) == tt) then 
            vtxtag(mesh%edge_1(ii)) = 1
            vtxtag(mesh%edge_2(ii)) = 1
        end if 
    end do 

    !Set quantities
    ncell = ncellzone(tt)
    nedge = nedgezone(tt)
    nvtx = sum(vtxtag(:))

    !Allocate mesh in this zone
    call allocate_mesh(mesh_par(tt),ncell,nedge,nvtx)

    !Index vertices in this zone 
    vtxtag(:) = 0 
    nvtx = 0 
    do ii=1,mesh%nedge
        if (edge_zone(ii) == tt) then 
            if (vtxtag(mesh%edge_1(ii)) == 0) then 
                nvtx = nvtx + 1
                vtxtag(mesh%edge_1(ii)) = nvtx
            end if
            if (vtxtag(mesh%edge_2(ii)) == 0) then 
                nvtx = nvtx + 1
                vtxtag(mesh%edge_2(ii)) = nvtx
            end if 
        end if 
    end do 

    !Assign base edge mesh items in this zone 
    nedge = 0 
    do ii=1,mesh%nedge
        if (edge_zone(ii) == tt) then 

            !Increment edge 
            nedge = nedge + 1

            !Edge cells and vertices
            v1 = mesh%edge_1(ii)
            v2 = mesh%edge_2(ii)
            cl = mesh%cell_left(ii)
            cr = mesh%cell_right(ii)

            !Assign this edge 
            mesh_par(tt)%edge_1(nedge) = vtxtag(v1)
            mesh_par(tt)%edge_2(nedge) = vtxtag(v2)
            if (cl .GT. 0) then 
                mesh_par(tt)%cell_left(nedge) = cell_zone_idx(cl)
                mesh_par(tt)%zone_left(nedge) = cellzone(cl)
            else
                mesh_par(tt)%cell_left(nedge) = cl
                mesh_par(tt)%zone_left(nedge) = 0 
            end if
            if (cr .GT. 0) then 
                mesh_par(tt)%cell_right(nedge) = cell_zone_idx(cr)
                mesh_par(tt)%zone_right(nedge) = cellzone(cr)
            else
                mesh_par(tt)%cell_right(nedge) = cr
                mesh_par(tt)%cell_right(nedge) = 0 
            end if
            
            !Assign vertices on this edge
            mesh_par(tt)%vtx_x(vtxtag(v1)) = mesh%vtx_x(v1)
            mesh_par(tt)%vtx_y(vtxtag(v1)) = mesh%vtx_y(v1)
            mesh_par(tt)%vtx_x(vtxtag(v2)) = mesh%vtx_x(v2)
            mesh_par(tt)%vtx_y(vtxtag(v2)) = mesh%vtx_y(v2)

            !Add properties to this edge 
            mesh_par(tt)%edge_dx(nedge) = mesh%edge_dx(ii)
            mesh_par(tt)%edge_dy(nedge) = mesh%edge_dy(ii)
            mesh_par(tt)%edgelen(nedge) = mesh%edgelen(ii)
            mesh_par(tt)%edge_mx(nedge) = mesh%edge_mx(ii)
            mesh_par(tt)%edge_my(nedge) = mesh%edge_my(ii)
            mesh_par(tt)%edge_nx(nedge) = mesh%edge_nx(ii)
            mesh_par(tt)%edge_ny(nedge) = mesh%edge_ny(ii)
            mesh_par(tt)%outflow_zone(nedge) = mesh%outflow_zone(ii)
        end if 
    end do 

    !Assign cell properties in this zone 
    do cc=1,mesh%ncell
        if (cellzone(cc) == tt) then 
            mesh_par(tt)%cell_mz(cell_zone_idx(cc)) = cell_mult_zone(cc)
            mesh_par(tt)%cellcenx(cell_zone_idx(cc)) = mesh%cellcenx(cc) 
            mesh_par(tt)%cellceny(cell_zone_idx(cc)) = mesh%cellceny(cc) 
            mesh_par(tt)%cell_nedge(cell_zone_idx(cc)) = mesh%cell_nedge(cc) 
            mesh_par(tt)%cell_nedgei(cell_zone_idx(cc)) = mesh%cell_nedgei(cc)
            mesh_par(tt)%cell_nedgeiw(cell_zone_idx(cc)) = mesh%cell_nedgeiw(cc)
            mesh_par(tt)%cell_vol(cell_zone_idx(cc)) = mesh%cell_vol(cc) 
            mesh_par(tt)%cell_elenint(cell_zone_idx(cc)) = mesh%cell_elenint(cc)
            mesh_par(tt)%cell_elentot(cell_zone_idx(cc)) = mesh%cell_elentot(cc)
            mesh_par(tt)%cell_link(cell_zone_idx(cc)) = cc 
        end if 
    end do 

    !Set outflow properties
    mesh_par(tt)%nzone_outflow = mesh%nzone_outflow
    mesh_par(tt)%nedge_of_zone(1:mesh%nzone_outflow) = mesh%nedge_of_zone(1:mesh%nzone_outflow)

    !Set active boundary conditions
    mesh_par(tt)%bc_active(:) = mesh%bc_active(:) 
end do 

!Export mesh cell zones if requested
if (options%export_mesh_partitions) then 
    call write_cell_dataPLT(options%iopath//'mesh_partitions',mesh,real(cellzone,dp))
end if
return 
end subroutine flow_par_segment_mesh




!Subroutine to split mesh zone into two approximatly equal zones =========================
subroutine split_half_zone(cellzone,Nzone,mesh,cellcenx,cellceny,zone2half)
implicit none 

!Variables - Import
integer(in) :: zone2half,Nzone
integer(in), dimension(:) :: cellzone
real(dp), dimension(:) :: cellcenx,cellceny
type(mesh_data) :: mesh

!Variables - Local 
integer(in) :: cc 
integer(in) :: Ncell_zone,cidx
integer(in), dimension(:), allocatable :: cellinzone
real(dp) :: xc,yc,qyy,qxx,azone,Ixx,Iyy,Ixy,thetap,cpval
real(dp) :: pv1(2),pv2(2),ccen(2)
real(dp) :: R(2,2)

!List cells in this zone 
Ncell_zone = 0 
do cc=1,mesh%ncell
    if (cellzone(cc) == zone2half) then 
        Ncell_zone = Ncell_zone + 1
    end if 
end do 
allocate(cellinzone(Ncell_zone))
cellinzone(:) = 0 
Ncell_zone = 0 
do cc=1,mesh%ncell
    if (cellzone(cc) == zone2half) then 
        Ncell_zone = Ncell_zone + 1
        cellinzone(Ncell_zone) = cc  
    end if 
end do 

!Evaluate centroid location 
qxx = 0.0d0 
qyy = 0.0d0 
azone = 0.0d0 
do cc=1,Ncell_zone
    cidx = cellinzone(cc)
    qyy = qyy + cellcenx(cidx)!*mesh%cell_vol(cidx)
    qxx = qxx + cellceny(cidx)!*mesh%cell_vol(cidx)
    azone = azone + 1.0d0 !mesh%cell_vol(cidx)
end do 
xc = qyy/azone
yc = qxx/azone

!Evaluate moments of area
Ixx = 0.0d0 
Iyy = 0.0d0 
Ixy = 0.0d0 
do cc=1,Ncell_zone
    cidx = cellinzone(cc)
    Ixx = Ixx + ((cellceny(cidx) - yc)**2)!*mesh%cell_vol(cidx)
    Iyy = Iyy + ((cellcenx(cidx) - xc)**2)!*mesh%cell_vol(cidx)
    Ixy = Ixy + ((cellcenx(cidx) - xc)*(cellceny(cidx) - yc))!*mesh%cell_vol(cidx)
end do 

!Find angle of principal axis 
thetap = 0.5d0*atan(2.0d0*Ixy/(Iyy - Ixx))

!Construct line along principal axis direction 
pv1(1) = -1.0d0 
pv1(2) = 0.0d0 
pv2(1) = 1.0d0 
pv2(2) = 0.0d0 
R(1,1) = cos(thetap)
R(1,2) = -sin(thetap)
R(2,1) = sin(thetap)
R(2,2) = cos(thetap)
pv1(1) = R(1,1)*pv1(1) + R(1,2)*pv1(2)
pv1(2) = R(2,1)*pv1(1) + R(2,2)*pv1(2)
pv2(1) = R(1,1)*pv2(1) + R(1,2)*pv2(2)
pv2(2) = R(2,1)*pv2(1) + R(2,2)*pv2(2)
pv1(1) = pv1(1) + xc
pv1(2) = pv1(2) + yc
pv2(1) = pv2(1) + xc
pv2(2) = pv2(2) + yc

!Divide cells by adding new zone 
Nzone = Nzone + 1
do cc=1,Ncell_zone
    cidx = cellinzone(cc)
    ccen(1) = cellcenx(cidx)
    ccen(2) = cellceny(cidx)
    cpval = cp2d(pv1,pv2,ccen)
    if (cpval .GT. 0.0d0) then 
        cellzone(cidx) = Nzone
    end if 
end do 
return 
end subroutine split_half_zone




!2D Cross product functions ===========================
function cp2d(A,B,C) result(val) !vectors (A->B) and (A->C)
implicit none 

!Variables - Import
real(dp) :: val
real(dp) :: A(2),B(2),C(2)

!Evaluate
val = (B(1) - A(1))*(C(2) - A(2)) - (B(2) - A(2))*(C(1) - A(1))
return 
end function cp2d


end module flow_mesh_mod