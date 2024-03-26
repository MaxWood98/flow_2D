!Flow solver boundary conditions subroutines module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.5.8
!Updated: 26/03/24

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
module flow_boundary_cond_mod
use flow_general_mod
contains


!================================================================================
!Wall and Internal Edge Conditions ==============================================
!================================================================================


!Calculate cell centre flux vectors subroutine =========================
subroutine cell_flux(F,G,flowvars,cellidx,zidx)
implicit none 

!Variables - Import
integer(in) :: cellidx,zidx
real(dp) :: F(4),G(4)
type(flow_var_data), dimension(:) :: flowvars

!Find flux vectors at cell centre 
F(1) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%u(cellidx)
F(2) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%u(cellidx)*flowvars(zidx)%u(cellidx) + flowvars(zidx)%p(cellidx)
F(3) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%u(cellidx)*flowvars(zidx)%v(cellidx)
F(4) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%u(cellidx)*flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx)*&
flowvars(zidx)%u(cellidx)
! F(4) = flowvars(zidx)%u(cellidx)*(flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx))

G(1) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)
G(2) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)*flowvars(zidx)%u(cellidx)
G(3) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)*flowvars(zidx)%v(cellidx) + flowvars(zidx)%p(cellidx)
G(4) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)*flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx)*&
flowvars(zidx)%v(cellidx)
! G(4) = flowvars(zidx)%v(cellidx)*(flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx))
return 
end subroutine cell_flux




!Constant extrapolation wall flux boundary condition subroutine =========================
subroutine wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
implicit none 

!Variables - Import
integer(in) :: cr,zr 
real(dp) :: pwall
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data), dimension(:) :: flowvars

!Set pressure at wall 
pwall = flowvars(zr)%p(cr)

!Set wall boundary condition
Fleft(1) = 0.0d0
Fleft(2) = flowvars(zr)%p(cr)
Fleft(3) = 0.0d0
Fleft(4) = 0.0d0

Gleft(1) = 0.0d0
Gleft(2) = 0.0d0
Gleft(3) = flowvars(zr)%p(cr)
Gleft(4) = 0.0d0
return 
end subroutine wall_flux_bc_Cextrp




!Linear extrapolation wall flux boundary condition subroutine =========================
subroutine wall_flux_bc_Lextrp(pwall,Fleft,Gleft,flowvars,mesh,edg_idx)
implicit none 

!Variables - Import
integer(in) :: edg_idx 
real(dp) :: pwall
real(dp), dimension(:) :: Fleft,Gleft
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Variables - Local
integer(in) :: cr 
real(dp) :: nx,ny,pD,nx1,ny1,nx2,ny2,lprj

!Cell right 
cr = mesh%cell_right(edg_idx)

!Vector from cell centre to wall magnitude 
nx1 = mesh%vtx_x(mesh%edge_1(edg_idx)) - mesh%cellcenx(cr)
ny1 = mesh%vtx_y(mesh%edge_1(edg_idx)) - mesh%cellceny(cr)
nx2 = mesh%vtx_x(mesh%edge_2(edg_idx)) - mesh%cellcenx(cr)
ny2 = mesh%vtx_y(mesh%edge_2(edg_idx)) - mesh%cellceny(cr)
nx = (nx1 + nx2)/2.0d0 
ny = (ny1 + ny2)/2.0d0
lprj = sqrt(nx**2 + ny**2)

!Normal vector of edge pointing out from the wall 
nx = mesh%edge_nx(edg_idx)*lprj
ny = mesh%edge_ny(edg_idx)*lprj

!Pressure directional derivative away from the wall 
pD = flowvars%pgrad(cr,1)*nx + flowvars%pgrad(cr,2)*ny

!Set pressure at wall 
pwall = flowvars%p(cr) - pD

!Set wall boundary condition
Fleft(1) = 0.0d0
Fleft(2) = pwall
Fleft(3) = 0.0d0
Fleft(4) = 0.0d0

Gleft(1) = 0.0d0
Gleft(2) = 0.0d0
Gleft(3) = pwall
Gleft(4) = 0.0d0
return 
end subroutine wall_flux_bc_Lextrp


!================================================================================
!Riemann Invarient Farfield Conditions ==========================================
!================================================================================


!Supersonic inflow boundary condition subroutine =========================
subroutine supersonic_inflow_bc_ri(Fleft,Gleft,cr,zr,edgidx,flowvars,mesh) 
implicit none 

!Variables - Import
integer(in) :: cr,zr,edgidx
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data),dimension(:) :: flowvars
type(mesh_data),dimension(:) :: mesh

!Variables - Local
real(dp) :: Rp,Rm,nx,ny
real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,vnorm_in,vnorm_inf,cin,cinf
real(dp) :: vbound,cbound
real(dp) :: ub,vb,rhob,pb,sb,eb

!Cell unit normal vector
nx = -mesh(zr)%edge_nx(edgidx)
ny = -mesh(zr)%edge_ny(edgidx)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Freestream primative vector
rhoinf = flowvars(zr)%rhoinf
uinf = flowvars(zr)%uinf
vinf = flowvars(zr)%vinf
pinf = flowvars(zr)%pinf
cinf = flowvars(zr)%cinf

!Boundary normal velocity 
vnorm_in = nx*uin + ny*vin
vnorm_inf = nx*uinf + ny*vinf

!R values
Rp = vnorm_inf + (2.0d0*cinf)/(flowvars(zr)%gam - 1.0d0)
Rm = vnorm_inf - (2.0d0*cinf)/(flowvars(zr)%gam - 1.0d0)

!Properties at boundary
vbound = 0.5d0*(Rp + Rm)
cbound = 0.25d0*(flowvars(zr)%gam - 1.0d0)*(Rp - Rm)

!Primatives at boundary
ub = uinf + (vbound - vnorm_inf)*nx
vb = vinf + (vbound - vnorm_inf)*ny
sb = (cinf**2)/(flowvars(zr)%gam*rhoinf**(flowvars(zr)%gam - 1.0d0))
rhob = ((cbound**2)/(flowvars(zr)%gam*sb))**(1.0d0/(flowvars(zr)%gam - 1.0d0))
pb = (rhob*cbound**2)/(flowvars(zr)%gam)

!Boundary flow energy
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return
end subroutine supersonic_inflow_bc_ri




!Supersonic outflow boundary condition subroutines =========================
subroutine supersonic_outflow_bc_ri(Fleft,Gleft,cr,zr,edgidx,flowvars,mesh)
implicit none 

!Variables - Import
integer(in) :: cr,zr,edgidx
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh

!Variables - Local
real(dp) :: Rp,Rm,nx,ny
real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,vnorm_in,vnorm_inf,cin,cinf
real(dp) :: vbound,cbound
real(dp) :: ub,vb,rhob,pb,sb,eb

!Cell normal vector
nx = -mesh(zr)%edge_nx(edgidx)
ny = -mesh(zr)%edge_ny(edgidx)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Freestream primative vector
rhoinf = flowvars(zr)%rhoinf
uinf = flowvars(zr)%uinf
vinf = flowvars(zr)%vinf
pinf = flowvars(zr)%pinf
cinf = flowvars(zr)%cinf

!Surface normal velocity
vnorm_in = nx*uin + ny*vin
vnorm_inf = nx*uinf + ny*vinf

!R values
Rp = vnorm_in + (2.0d0*cin)/(flowvars(zr)%gam - 1.0d0)
Rm = vnorm_in - (2.0d0*cin)/(flowvars(zr)%gam - 1.0d0)

!Properties at boundary
vbound = 0.5d0*(Rp + Rm)
cbound = 0.25d0*(flowvars(zr)%gam - 1.0d0)*(Rp - Rm)

!Primatives at boundary
ub = (uin + (vbound - vnorm_in)*nx)
vb = (vin + (vbound - vnorm_in)*ny)
sb = (cin**2)/(flowvars(zr)%gam*rhoin**(flowvars(zr)%gam - 1.0d0))
rhob = ((cbound**2)/(flowvars(zr)%gam*sb))**(1.0d0/(flowvars(zr)%gam - 1.0d0))
pb = (rhob*cbound**2)/(flowvars(zr)%gam)

!Boundary flow energy
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return
end subroutine supersonic_outflow_bc_ri




!Subsonic inflow boundary condition subroutine =========================
subroutine subsonic_inflow_bc_ri(Fleft,Gleft,cr,zr,edgidx,flowvars,mesh)
implicit none 

!Variables - Import
integer(in) :: cr,zr,edgidx
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh

!Variables - Local
real(dp) :: Rp,Rm,nx,ny
real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,vnorm_in,vnorm_inf,cin,cinf
real(dp) :: vbound,cbound
real(dp) :: ub,vb,rhob,pb,sb,eb

!Cell normal vector
nx = -mesh(zr)%edge_nx(edgidx)
ny = -mesh(zr)%edge_ny(edgidx)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Freestream primative vector
rhoinf = flowvars(zr)%rhoinf
uinf = flowvars(zr)%uinf
vinf = flowvars(zr)%vinf
pinf = flowvars(zr)%pinf
cinf = flowvars(zr)%cinf

!Surface normal velocity
vnorm_in = nx*uin + ny*vin
vnorm_inf = nx*uinf + ny*vinf

!R values
Rp = vnorm_in + (2.0d0*cin)/(flowvars(zr)%gam - 1.0d0)
Rm = vnorm_inf - (2.0d0*cinf)/(flowvars(zr)%gam - 1.0d0)

!Properties at boundary
vbound = 0.5d0*(Rp + Rm)
cbound = 0.25d0*(flowvars(zr)%gam - 1.0d0)*(Rp - Rm)

!Primatives at boundary
ub = uinf + (vbound - vnorm_inf)*nx
vb = vinf + (vbound - vnorm_inf)*ny
sb = (cinf**2)/(flowvars(zr)%gam*rhoinf**(flowvars(zr)%gam - 1.0d0))
rhob = ((cbound**2)/(flowvars(zr)%gam*sb))**(1.0d0/(flowvars(zr)%gam - 1.0d0))
pb = (rhob*cbound**2)/(flowvars(zr)%gam)

!Boundary flow energy
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return
end subroutine subsonic_inflow_bc_ri




!Subsonic outflow boundary condition subroutine =========================
subroutine subsonic_outflow_bc_ri(Fleft,Gleft,cr,zr,edgidx,flowvars,mesh)
implicit none 

!Variables - Import
integer(in) :: cr,zr,edgidx
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh

!Variables - Local
real(dp) :: Rp,Rm,nx,ny
real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,vnorm_in,vnorm_inf,cin,cinf
real(dp) :: vbound,cbound
real(dp) :: ub,vb,rhob,pb,sb,eb

!Cell normal vector
nx = -mesh(zr)%edge_nx(edgidx)
ny = -mesh(zr)%edge_ny(edgidx)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Freestream primative vector
rhoinf = flowvars(zr)%rhoinf
uinf = flowvars(zr)%uinf
vinf = flowvars(zr)%vinf
pinf = flowvars(zr)%pinf
cinf = flowvars(zr)%cinf

!Surface normal velocity
vnorm_in = nx*uin + ny*vin
vnorm_inf = nx*uinf + ny*vinf

!R values
Rp = vnorm_in + (2.0d0*cin)/(flowvars(zr)%gam - 1.0d0)
Rm = vnorm_inf - (2.0d0*cinf)/(flowvars(zr)%gam - 1.0d0)

!Properties at boundary
vbound = 0.5d0*(Rp + Rm)
cbound = 0.25d0*(flowvars(zr)%gam - 1.0d0)*(Rp - Rm)

!Primatives at boundary 
ub = (uin + (vbound - vnorm_in)*nx)
vb = (vin + (vbound - vnorm_in)*ny)
sb = (cin**2)/(flowvars(zr)%gam*rhoin**(flowvars(zr)%gam - 1.0d0))
rhob = ((cbound**2)/(flowvars(zr)%gam*sb))**(1.0d0/(flowvars(zr)%gam - 1.0d0))
pb = (rhob*cbound**2)/(flowvars(zr)%gam)

!Boundary flow energy
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return
end subroutine subsonic_outflow_bc_ri


!================================================================================
!Characteristic Farfield Conditions =============================================
!================================================================================


!Subsonic inflow/outflow boundary condition subroutine =========================
subroutine subsonic_bc_char(Fleft,Gleft,cr,zr,flowvars,if_of,mesh,edge_idx)
implicit none 

!Variables - Import
integer(in) :: cr,zr,if_of,edge_idx
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh

!Variables - Local
real(dp) :: rhoin,pin,rhoinf,pinf,cin,cinf
real(dp) :: ub,vb,rhob,pb,eb,nx,ny
real(dp) :: vel_in(2),vel_inf(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2),basis_ax(2),basis_ay(2)
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4),Mb2a(2,2),Ma2b(2,2)

!Cell normal vector with positive dot product with the velocity vector 
if (if_of == 1) then !inflow
    nx = mesh(zr)%edge_nx(edge_idx)
    ny = mesh(zr)%edge_ny(edge_idx)
elseif (if_of == -1) then !outflow
    nx = -mesh(zr)%edge_nx(edge_idx)
    ny = -mesh(zr)%edge_ny(edge_idx)
end if 

!Define global coordiante system basis 
basis_bx(1) = 1.0d0 
basis_bx(2) = 0.0d0 
basis_by(1) = 0.0d0 
basis_by(2) = 1.0d0 

!Define edge normal coordiante system basis 
basis_ax(1) = nx
basis_ax(2) = ny
basis_ay(1) = mesh(zr)%edge_dx(edge_idx)
basis_ay(2) = mesh(zr)%edge_dy(edge_idx)
basis_ay(:) = basis_ay(:)/norm2(basis_ay(:))

!Get change of basis from global to edge normal coordinate system 
call get_basis_change(Mb2a,Ma2b,basis_bx,basis_by,basis_ax,basis_ay)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
vel_in(1) = flowvars(zr)%u(cr)
vel_in(2) = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Freestream primative vector
rhoinf = flowvars(zr)%rhoinf
vel_inf(1) = flowvars(zr)%uinf
vel_inf(2) = flowvars(zr)%vinf
pinf = flowvars(zr)%pinf
cinf = flowvars(zr)%cinf

!Change velocities to edge normal basis
vel_in_n = change_basis(Mb2a,vel_in)
vel_inf_n = change_basis(Mb2a,vel_inf)

!Evaluate characteristic matrix
! call charmat(chm,chmi,rhoinf,cinf)
call charmat(chm,chmi,rhoin,cin)

!Evaluate flow deltas from freestream
fdeltasi(1,1) = rhoin - rhoinf
fdeltasi(2,1) = vel_in_n(1) - vel_inf_n(1) 
fdeltasi(3,1) = vel_in_n(2) - vel_inf_n(2)  
fdeltasi(4,1) = pin - pinf

!Evaluate characteristics
fchars = matmul(chm,fdeltasi)

!Apply boundary condition 
if (if_of == 1) then !inflow
    fchars(1:3,1) = 0.0d0 
elseif (if_of == -1) then !outflow
    fchars(4,1) = 0.0d0 
end if 

!Calculate boundary deltas
fdeltasb = matmul(chmi,fchars)

!Evaluate boundary velocity
vel_b_n(1) = fdeltasb(2,1) + vel_inf_n(1) 
vel_b_n(2) = fdeltasb(3,1) + vel_inf_n(2) 

!Change boundary velocity to the global basis
vel_b = change_basis(Ma2b,vel_b_n)
ub = vel_b(1)
vb = vel_b(2)

!Evaluate boundary primatives 
rhob = fdeltasb(1,1) + rhoinf
pb = fdeltasb(4,1) + pinf
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return 
end subroutine subsonic_bc_char




!Characteristic matrix and inverse subroutines =========================
subroutine charmat(chm,chmi,rho,c)
implicit none 

!Variables - Import 
real(dp) :: rho,c
real(dp) :: chm(4,4),chmi(4,4)

!Base matrix 
chm(1,1) = -c*c
chm(1,2) = 0.0d0 
chm(1,3) = 0.0d0 
chm(1,4) = 1.0d0
chm(2,1) = 0.0d0 
chm(2,2) = 0.0d0 
chm(2,3) = rho*c
chm(2,4) = 0.0d0 
chm(3,1) = 0.0d0 
chm(3,2) = rho*c
chm(3,3) = 0.0d0 
chm(3,4) = 1.0d0 
chm(4,1) = 0.0d0 
chm(4,2) = -rho*c
chm(4,3) = 0.0d0 
chm(4,4) = 1.0d0

!Inverse matrix 
chmi(1,1) = -1.0d0/(c*c)
chmi(1,2) = 0.0d0 
chmi(1,3) = 1.0d0/(2.0d0*c*c)
chmi(1,4) = 1.0d0/(2.0d0*c*c)
chmi(2,1) = 0.0d0 
chmi(2,2) = 0.0d0 
chmi(2,3) = 1.0d0/(2.0d0*rho*c)
chmi(2,4) = -1.0d0/(2.0d0*rho*c)
chmi(3,1) = 0.0d0 
chmi(3,2) = 1.0d0/(rho*C)
chmi(3,3) = 0.0d0 
chmi(3,4) = 0.0d0 
chmi(4,1) = 0.0d0 
chmi(4,2) = 0.0d0 
chmi(4,3) = 0.5d0
chmi(4,4) = 0.5d0
return 
end subroutine charmat


!================================================================================
!Internal Inflow/Outflow Conditions ==============================================
!================================================================================


!Stagnation state inflow boundary condition subroutine =========================
subroutine stagnation_inflow_bc(Fleft,Gleft,flowvars,mesh,edge_idx,cr,zr)
implicit none 

!Variables - Import
integer(in) :: edge_idx,cr,zr
real(dp), dimension(:) :: Fleft,Gleft 
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh

!Variables - Local
real(dp) :: rhoin,pin,rhoinf,pinf,cin,tinf,vnorm_in
real(dp) :: ub,vb,rhob,pb,eb,Mb,nx,ny,H,Rm
real(dp) :: a,b,c,cbp,cbm,cb,vninf,velinf
real(dp) :: vel_in(2),vel_in_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2),basis_ax(2),basis_ay(2)
real(dp) :: Mb2a(2,2),Ma2b(2,2)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
vel_in(1) = flowvars(zr)%u(cr)
vel_in(2) = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Boundary normal velocity
vnorm_in = vel_in(1)*mesh(zr)%edge_nx(edge_idx) + vel_in(2)*mesh(zr)%edge_ny(edge_idx)

!Cell normal vector with positive dot product with the velocity vector 
if (vnorm_in .GE. 0.0d0) then !inflow
    nx = mesh(zr)%edge_nx(edge_idx)
    ny = mesh(zr)%edge_ny(edge_idx)
else !backflow 
    nx = -mesh(zr)%edge_nx(edge_idx)
    ny = -mesh(zr)%edge_ny(edge_idx)
end if 

!Define global coordiante system basis 
basis_bx(1) = 1.0d0 
basis_bx(2) = 0.0d0 
basis_by(1) = 0.0d0 
basis_by(2) = 1.0d0 

!Define edge normal coordiante system basis 
basis_ax(1) = nx
basis_ax(2) = ny
basis_ay(1) = mesh(zr)%edge_dx(edge_idx)
basis_ay(2) = mesh(zr)%edge_dy(edge_idx)
basis_ay(:) = basis_ay(:)/norm2(basis_ay(:))

!Get change of basis from global to edge normal coordinate system 
call get_basis_change(Mb2a,Ma2b,basis_bx,basis_by,basis_ax,basis_ay)

!Change internal velocity to edge normal basis
vel_in_n = change_basis(Mb2a,vel_in)

!Enthalpy (from interior properties)
H = ((cin*cin)/(flowvars(zr)%gam - 1.0d0)) + 0.5d0*(vel_in_n(1)*vel_in_n(1) + vel_in_n(2)*vel_in_n(2)) 

!Reimann invariant - 
Rm = -vel_in_n(1) + (2.0d0*cin/(flowvars(zr)%gam - 1.0d0))

!Solve for boundary speed of sound 
a = 1.0d0 + (2.0d0/(flowvars(zr)%gam - 1.0d0))
b = -2.0d0*Rm
c = 0.5d0*(flowvars(zr)%gam - 1.0d0)*(Rm*Rm - 2.0d0*H)
cbp = (-b + sqrt(b*b - 4.0d0*a*c))/(2.0d0*a)
cbm = (-b - sqrt(b*b - 4.0d0*a*c))/(2.0d0*a)
cb = max(cbp,cbm)

!Solve for boundary normal velocity and mach number 
vninf = (2.0d0*cb/(flowvars(zr)%gam - 1.0d0)) - Rm
velinf = sqrt(vninf**2 + vel_in_n(2)**2)
Mb = velinf/cb

!Boundary pressure temperature and density
pinf = flowvars(zr)%p0inf/(((1.0d0 + 0.5d0*(flowvars(zr)%gam - 1.0d0)*Mb*Mb)**&
(flowvars(zr)%gam/(flowvars(zr)%gam - 1.0d0))))
tinf = flowvars(zr)%t0inf/(1.0d0 + 0.5d0*(flowvars(zr)%gam - 1.0d0)*Mb*Mb)
rhoinf = pinf/(tinf*flowvars(zr)%R) 

!Evaluate boundary velocity
vel_b_n(1) = vninf
vel_b_n(2) = vel_in_n(2)

!Change boundary velocity to the global basis
vel_b = change_basis(Ma2b,vel_b_n)
ub = vel_b(1)
vb = vel_b(2)

!Evaluate boundary primatives 
rhob = rhoinf
pb = pinf
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)

!Accumulate inflow mass flux properties
flowvars(zr)%mflux_in = flowvars(zr)%mflux_in + rhob*vel_b_n(1)*mesh(zr)%edgelen(edge_idx)
flowvars(zr)%time_in = flowvars(zr)%time_in + mesh(zr)%deltaT(cr)!*mesh(zr)%cell_vol(cr)

!Accumulate number of cells 
flowvars(zr)%cvol_in = flowvars(zr)%cvol_in + 1.0d0 !mesh(zr)%cell_vol(cr)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return 
end subroutine stagnation_inflow_bc




!Freestream inflow boundary condition subroutine =========================
subroutine freestream_inflow_c(Fleft,Gleft,cr,zr,flowvars,mesh,options,edge_idx)
implicit none 

!Variables - Import
integer(in) :: cr,zr,edge_idx
real(dp), dimension(:) :: Fleft,Gleft
type(mesh_data), dimension(:) :: mesh
type(flow_var_data), dimension(:) :: flowvars
type(options_data) :: options 

!Variables - Local
real(dp) :: rhoin,pin,rhoinf,pinf,cin,cinf,vnorm_in
real(dp) :: ub,vb,rhob,pb,eb,pwall,nx,ny
real(dp) :: vel_in(2),vel_inf(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2),basis_ax(2),basis_ay(2)
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4),Mb2a(2,2),Ma2b(2,2)

!Cell normal vector with positive dot product with the velocity vector 
nx = mesh(zr)%edge_nx(edge_idx)
ny = mesh(zr)%edge_ny(edge_idx)

!Define global coordiante system basis 
basis_bx(1) = 1.0d0 
basis_bx(2) = 0.0d0 
basis_by(1) = 0.0d0 
basis_by(2) = 1.0d0 

!Define edge normal coordiante system basis 
basis_ax(1) = nx
basis_ax(2) = ny
basis_ay(1) = mesh(zr)%edge_dx(edge_idx)
basis_ay(2) = mesh(zr)%edge_dy(edge_idx)
basis_ay(:) = basis_ay(:)/norm2(basis_ay(:))

!Get change of basis from global to edge normal coordinate system 
call get_basis_change(Mb2a,Ma2b,basis_bx,basis_by,basis_ax,basis_ay)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
vel_in(1) = flowvars(zr)%u(cr)
vel_in(2) = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Freestream primative vector
rhoinf = flowvars(zr)%rhoinf
vel_inf(1) = flowvars(zr)%uinf
vel_inf(2) = flowvars(zr)%vinf
pinf = flowvars(zr)%pinf
cinf = flowvars(zr)%cinf

!Change velocities to edge normal basis
vel_in_n = change_basis(Mb2a,vel_in)
vel_inf_n = change_basis(Mb2a,vel_inf)

!Boundary normal velocity
vnorm_in = vel_in(1)*mesh(zr)%edge_nx(edge_idx) + vel_in(2)*mesh(zr)%edge_ny(edge_idx)

!Switch on sign of normal velocity
if (vnorm_in .LT. 0.0d0) then !backflow so default to wall condition
    call wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
else

    !Evaluate characteristic matrix
    call charmat(chm,chmi,rhoinf,cinf)
    ! call charmat(chm,chmi,rhoin,cin)

    !Evaluate flow deltas from freestream
    fdeltasi(1,1) = rhoin - rhoinf
    fdeltasi(2,1) = vel_in_n(1) - vel_inf_n(1) 
    fdeltasi(3,1) = vel_in_n(2) - vel_inf_n(2)  
    fdeltasi(4,1) = pin - pinf

    !Evaluate characteristics
    fchars = matmul(chm,fdeltasi)

    !Apply boundary condition (inflow)
    fchars(1:3,1) = 0.0d0
    ! fchars(1:3,1) = 0.1d0*fchars(1:3,1) 
    
    !Calculate boundary deltas
    fdeltasb = matmul(chmi,fchars)

    !Evaluate boundary velocity
    vel_b_n(1) = fdeltasb(2,1) + vel_inf_n(1) 
    vel_b_n(2) = fdeltasb(3,1) + vel_inf_n(2) 

    !Change boundary velocity to the global basis
    vel_b = change_basis(Ma2b,vel_b_n)
    ub = vel_b(1)
    vb = vel_b(2)

    !Evaluate boundary primatives 
    rhob = fdeltasb(1,1) + rhoinf
    if (options%force_fixed_pratio) then 
        pb = pinf
    else
        pb = fdeltasb(4,1) + pinf
    end if 
    eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
    ! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

    !Accumulate inflow mass flux properties
    flowvars(zr)%mflux_in = flowvars(zr)%mflux_in + rhob*vel_b_n(1)*mesh(zr)%edgelen(edge_idx)
    flowvars(zr)%time_in = flowvars(zr)%time_in + mesh(zr)%deltaT(cr)!*mesh(zr)%cell_vol(cr)

    !Construct boundary fluxes
    Fleft(1) = rhob*ub
    Fleft(2) = rhob*ub*ub + pb
    Fleft(3) = rhob*ub*vb
    Fleft(4) = rhob*ub*eb + pb*ub
    ! Fleft(4) = ub*(eb + pb)

    Gleft(1) = rhob*vb
    Gleft(2) = rhob*vb*ub
    Gleft(3) = rhob*vb*vb + pb
    Gleft(4) = rhob*vb*eb + pb*vb
    ! Gleft(4) = vb*(eb + pb)
end if 

!Accumulate number of cells 
flowvars(zr)%cvol_in = flowvars(zr)%cvol_in + 1.0d0 !mesh(zr)%cell_vol(cr)
return 
end subroutine freestream_inflow_c




!Mass flux outflow boundary condition (velocity) subroutine =========================
subroutine mass_flux_outflow_bc_V(R,mesh,flowvars,options,iter,rk_stage,upd_allow,zone)
implicit none 

!Variables - Import 
integer(in) :: iter,rk_stage,upd_allow,zone
real(dp), dimension(:,:) :: R
type(mesh_data), dimension(:) :: mesh
type(flow_var_data), dimension(:) :: flowvars
type(options_data) :: options

!Variables - Local 
integer(in) :: ii,tt,cl,cr,zl,zr
real(dp) :: nx,ny,rhoin,uin,vin,pin,velNi,vel_in,Vn_new
real(dp) :: A_rho_out,mflux_out,mflux_outF,mflux_error,vnset_of
real(dp) :: Sb,Jm,Jp,cin,M_in,M_b
real(dp) :: cb,rhob,pb,ub,vb,eb,velb,un,vn,ut,vt
real(dp) :: Fleft(4),Fright(4),Gleft(4),Gright(4)
real(dp) :: Fedge(4),Gedge(4),Flxe(4)

!Check and update outflow velocity if correct iteration 
if (upd_allow == 1) then 
    if ((mod(iter,options%Mflux_relaxiter) == 0) .AND. (rk_stage == 1)) then 

        !Find outflow mass flux in this zone 
        A_rho_out = 0.0d0 
        do ii=1,mesh(zone)%nedge

            !Left and right cells for this edge
            cl = mesh(zone)%cell_left(ii)
            cr = mesh(zone)%cell_right(ii)

            !Accumulate outflow density*area
            if (cl == -4) then 
                A_rho_out = A_rho_out + mesh(zone)%edgelen(ii)*flowvars(zone)%rho(cr)
            end if 
        end do 
        mflux_out = A_rho_out*flowvars(zone)%vnset_of

        !Evaluate mass flux error and set outflow velocity 
        !$OMP barrier
        !$OMP single

        !Accumulate flux between zones 
        mflux_outF = 0.0d0 
        do tt=1,options%num_threads
            mflux_outF = mflux_outF + mflux_out
        end do 
        
        !Evaluate mass flux error and push to each zone 
        mflux_out = mflux_outF
        mflux_error = abs(mflux_out - flowvars(zone)%mflow_in)/flowvars(zone)%mflow_in
        do tt=1,options%num_threads
            flowvars(tt)%mflux_error = mflux_error
        end do 

        !Update outflow velocity if outside tolerance on valid iterations 
        if (mflux_error .GT. options%Mflux_toll) then !if outside tollerance -> update outflow velocity
            Vn_new = flowvars(zone)%mflow_in/A_rho_out
            vnset_of = flowvars(zone)%vnset_of + options%Mflux_relax*(Vn_new - flowvars(zone)%vnset_of)
            do tt=1,options%num_threads
                flowvars(tt)%vnset_of = vnset_of
            end do 
            !print *,flowvars%vnset_of,' ========================='
        end if
        !$OMP end single 
    end if 
end if 
!$OMP barrier
! print *,flowvars%vnset_of

!Set downstream properties
do ii=1,mesh(zone)%nedge

    !Left and right cells for this edge
    cl = mesh(zone)%cell_left(ii)
    cr = mesh(zone)%cell_right(ii)
    zl = mesh(zone)%zone_left(ii)
    zr = mesh(zone)%zone_right(ii) 

    !If outflow face
    if (cl == -4) then 

        !Edge unit normal vector
        nx = -mesh(zone)%edge_nx(ii)
        ny = -mesh(zone)%edge_ny(ii)

        !Local primative vector
        rhoin = flowvars(zone)%rho(cr) 
        pin = flowvars(zone)%p(cr)
        uin = flowvars(zone)%u(cr)
        vin = flowvars(zone)%v(cr) 
        cin = sqrt(flowvars(zone)%gam*(pin/rhoin))
        vel_in = sqrt(uin*uin + vin*vin)
        velNi = uin*nx + vin*ny

        !Local boundary normal and tangent velocity components 
        un = uin*nx
        vn = vin*ny
        ut = uin - un
        vt = vin - vn 

        !Boundary entropy
        ! Sb = pin/(rhoin**flowvars(zone)%gam)
        Sb = (cin**2)/(flowvars(zone)%gam*rhoin**(flowvars(zone)%gam - 1.0d0))

        !Outgoing and incoming invariant 
        Jm = velNi + ((2.0d0*cin)/(flowvars(zone)%gam - 1.0d0))
        Jp = 2.0d0*flowvars(zone)%vnset_of - Jm

        !Boundary speed of sound
        cb = 0.25d0*(flowvars(zone)%gam - 1.0d0)*(Jm - Jp)

        !Internal and boundary mach number 
        M_in = velNi/cin
        M_b = flowvars(zone)%vnset_of/cb

        !Cases to set cell left flux 
        if ((M_in .LT. 1.0d0) .AND. (M_b .LT. 1.0d0)) then !both subsonic -> normal state, use calculated subsonic outflow velocity 

            !Use set boundary normal velocity
            velb = flowvars(zone)%vnset_of 

            !Set boundary properties
            rhob = ((cb*cb)/(flowvars(zone)%gam*Sb))**(1.0d0/(flowvars(zone)%gam - 1.0d0))
            pb = (rhob*cb*cb)/flowvars(zone)%gam 
            ub = velb*nx + ut 
            vb = velb*ny + vt 
            eb = (pb/((options%gamma - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
            ! eb = (pb/(options%gamma - 1.0d0)) + 0.5d0*rhob*(ub**2 + vb**2)

            !Construct boundary fluxes 
            Fleft(1) = rhob*ub
            Fleft(2) = rhob*ub*ub + pb
            Fleft(3) = rhob*ub*vb
            Fleft(4) = rhob*ub*eb + pb*ub
            ! Fleft(4) = ub*(eb + pb)

            Gleft(1) = rhob*vb
            Gleft(2) = rhob*vb*ub
            Gleft(3) = rhob*vb*vb + pb
            Gleft(4) = rhob*vb*eb + pb*vb
            ! Gleft(4) = vb*(eb + pb)
        elseif ((M_in .LT. 1.0d0) .AND. (M_b .GE. 1.0d0)) then !internal subsonic / boundary supersonic -> set boundary to sonic to choke 

            !Force boundary normal velocity to sonic 
            velb = cb

            !Set boundary properties
            rhob = ((cb*cb)/(flowvars(zone)%gam*Sb))**(1.0d0/(flowvars(zone)%gam - 1.0d0))
            pb = (rhob*cb*cb)/flowvars(zone)%gam 
            ub = velb*nx + ut 
            vb = velb*ny + vt
            eb = (pb/((options%gamma - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
            ! eb = (pb/(options%gamma - 1.0d0)) + 0.5d0*rhob*(ub**2 + vb**2)

            !Construct boundary fluxes 
            Fleft(1) = rhob*ub
            Fleft(2) = rhob*ub*ub + pb
            Fleft(3) = rhob*ub*vb
            Fleft(4) = rhob*ub*eb + pb*ub
            ! Fleft(4) = ub*(eb + pb)

            Gleft(1) = rhob*vb
            Gleft(2) = rhob*vb*ub
            Gleft(3) = rhob*vb*vb + pb
            Gleft(4) = rhob*vb*eb + pb*vb
            ! Gleft(4) = vb*(eb + pb)
        elseif ((M_in .GE. 1.0d0) .AND. (M_b .LT. 1.0d0)) then !internal supersonic / boundary subsonic -> use supersonic boundary
            ! call supersonic_internal_outflow(Fleft,Gleft,flowvars,cr)
            call supersonic_inflow_bc_ri(Fleft,Gleft,cr,zr,ii,flowvars,mesh) 
        elseif ((M_in .GT. 1.0d0) .AND. (M_b .GT. 1.0d0)) then !both supersonic -> use supersonic boundary
            ! call supersonic_internal_outflow(Fleft,Gleft,flowvars,cr)
            call supersonic_inflow_bc_ri(Fleft,Gleft,cr,zr,ii,flowvars,mesh) 
        end if 

        !Cell right flux
        call cell_flux(Fright,Gright,flowvars,cr,zr)

        !Edge flux value ([dy,-dx] normal)
        Fedge(:) = 0.5d0*(Fleft(:) + Fright(:))*mesh(zone)%edge_dy(ii)
        Gedge(:) = 0.5d0*(Gleft(:) + Gright(:))*mesh(zone)%edge_dx(ii)
        Flxe(:) =  Fedge(:) - Gedge(:)

        !Accumulate to cell right 
        if (mesh(zr)%cell_mz(cr) == 0) then
            R(cr,:) = R(cr,:) - Flxe(:)
        else
            !$OMP atomic 
            R(cr,1) = R(cr,1) - Flxe(1)
            !$OMP atomic 
            R(cr,2) = R(cr,2) - Flxe(2)
            !$OMP atomic 
            R(cr,3) = R(cr,3) - Flxe(3)
            !$OMP atomic 
            R(cr,4) = R(cr,4) - Flxe(4)
        end if 
    end if
end do 
return 
end subroutine mass_flux_outflow_bc_V




!Subroutine to set initial mass flux properties =========================
subroutine set_initial_massflux(flowvars,mesh,options)
implicit none 

!Variables - Import 
type(mesh_data) :: mesh
type(options_data) :: options
type(flow_var_data) :: flowvars

!Variables - Local 
integer(in) :: ii,cl,cr

!Evaluate inflow and outflow areas
flowvars%A_in = 0.0d0 
flowvars%A_out = 0.0d0 
do ii=1,mesh%nedge

    !Left and right cells for this edge
    cl = mesh%cell_left(ii)
    cr = mesh%cell_right(ii)

    !Accumulate inflow area
    if ((cl == -3) .OR. (cl == -5)) then !If inflow face
        flowvars%A_in = flowvars%A_in + mesh%edgelen(ii)
    end if 

    !Accumulate outflow area
    if ((cl == -4) .OR. (cl == -6)) then !If outflow face
        flowvars%A_out = flowvars%A_out + mesh%edgelen(ii)
    end if 
end do 

!Set inflow mass flux 
flowvars%mflow_in = flowvars%rhoinf*flowvars%A_in*flowvars%velinf

!Set outflow velocity for equal outflow mass flux
if (options%init_state .NE. 'r') then 
    flowvars%vnset_of = flowvars%mflow_in/(flowvars%A_out*flowvars%rhoinf)
end if 

!Set initial outflow pressure  
if (options%init_state .NE. 'r') then 
    flowvars%mflux_pset = flowvars%pinf
end if  
return 
end subroutine set_initial_massflux




!Subroutine to evaluate mass flux error =========================
subroutine mass_flux_error(mesh,flowvars)
implicit none 

!Variables - Import 
type(mesh_data) :: mesh
type(flow_var_data) :: flowvars

!Variables - Local 
integer(in) :: ii,cl,cr
real(dp) :: A_rho_out,mflux_out

!Find outflow mass flux 
A_rho_out = 0.0d0 
do ii=1,mesh%nedge

    !Left and right cells for this edge
    cl = mesh%cell_left(ii)
    cr = mesh%cell_right(ii)

    !Accumulate outflow density*area and average velocity
    if (cl == -4) then 
        A_rho_out = A_rho_out + mesh%edgelen(ii)*flowvars%rho(cr)
    end if 
end do 
mflux_out = A_rho_out*flowvars%vnset_of

!Evaluate mass flux error
flowvars%mflux_error = abs(mflux_out - flowvars%mflow_in)/flowvars%mflow_in
return 
end subroutine mass_flux_error




!Outflow back pressure boundary condition (Riemann) ========================= (not working ******)
subroutine backpressure_outflow_bc_ri(Fleft,Gleft,flowvars,mesh,edge_idx,cr,zr) 
implicit none 

!Variables - Import
integer(in) :: cr,zr,edge_idx
real(dp), dimension(:) :: Fleft,Gleft 
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh

!Variables - Local

real(dp) :: rhoin,uin,vin,pin,cin,Machin,velin,vnorm_in,unin,vnin,utin,vtin
real(dp) :: Sb,pb,rhob,ub,vb,eb,cb,nx,ny
real(dp) :: Jm,Jp,velen,pinf,Vnb 


!Edge normal 
nx = mesh(zr)%edge_nx(edge_idx)
ny = mesh(zr)%edge_ny(edge_idx)

!Internal primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
velin = sqrt(uin*uin + vin*vin)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))
Machin = velin/cin
velen = nx*uin + ny*vin

!Internal boundary normal and tangential velocity
vnorm_in = uin*nx + vin*ny
unin = nx*vnorm_in
vnin = ny*vnorm_in
utin = uin - unin
vtin = vin - vnin

!Set downstream pressure
if (Machin .GE. 1.0d0) then 
    pinf = pin
else
    pinf = flowvars(zr)%backp_pset
end if 
pb = pinf 

!Boundary entropy 
Sb = pin/(rhoin**flowvars(zr)%gam)

!Exiting invarient 
Jm = -velen + (2.0d0*cin/(flowvars(zr)%gam - 1.0d0))

!Boundary speed of sound 
cb = sqrt(flowvars(zr)%gam*(pinf**((flowvars(zr)%gam - 1.0d0)/flowvars(zr)%gam))*(Sb**(1.0d0/flowvars(zr)%gam)))

!Entering invarient 
Jp = Jm - (4.0d0/(flowvars(zr)%gam - 1.0d0))*cb

!Boundary density
rhob = (cb*cb/(flowvars(zr)%gam*Sb))**(1.0d0/(flowvars(zr)%gam - 1.0d0))

!Boundary normal velocity 
Vnb = -0.5d0*(Jm + Jp)

!Boundary velocity 
ub = Vnb*nx + utin 
vb = Vnb*ny + vtin

!Boundary energy
eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return 
end subroutine backpressure_outflow_bc_ri




!Outflow back pressure boundary condition (Characteristic) =========================
subroutine backpressure_outflow_bc_c(Fleft,Gleft,flowvars,mesh,options,edge_idx,cr,zr) 
implicit none 

!Variables - Import
integer(in) :: cr,zr,edge_idx
real(dp), dimension(:) :: Fleft,Gleft 
type(flow_var_data), dimension(:) :: flowvars
type(mesh_data), dimension(:) :: mesh
type(options_data) :: options 

!Variables - Local
integer(in) :: zof
real(dp) :: rhoin,pin,rhoinf,pinf,tin,cin,Machin,cinf,vnorm_in
real(dp) :: ub,vb,rhob,pb,eb,nx,ny
real(dp) :: vel_in(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2),basis_ax(2),basis_ay(2)
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4),Mb2a(2,2),Ma2b(2,2)

!Outflow zone of this edge
zof = mesh(zr)%outflow_zone(edge_idx)

!Local primative vector
rhoin = flowvars(zr)%rho(cr) 
vel_in(1) = flowvars(zr)%u(cr)
vel_in(2) = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))
tin = pin/(rhoin*flowvars(zr)%R)

!Boundary normal velocity
vnorm_in = vel_in(1)*mesh(zr)%edge_nx(edge_idx) + vel_in(2)*mesh(zr)%edge_ny(edge_idx)

!Cell normal vector with positive dot product with the velocity vector 
if (vnorm_in .GE. 0.0d0) then !backflow
    nx = mesh(zr)%edge_nx(edge_idx)
    ny = mesh(zr)%edge_ny(edge_idx)
else !outflow 
    nx = -mesh(zr)%edge_nx(edge_idx)
    ny = -mesh(zr)%edge_ny(edge_idx)
end if 

!Define global coordiante system basis 
basis_bx(1) = 1.0d0 
basis_bx(2) = 0.0d0 
basis_by(1) = 0.0d0 
basis_by(2) = 1.0d0 

!Define edge normal coordiante system basis 
basis_ax(1) = nx
basis_ax(2) = ny
basis_ay(1) = mesh(zr)%edge_dx(edge_idx)
basis_ay(2) = mesh(zr)%edge_dy(edge_idx)
basis_ay(:) = basis_ay(:)/norm2(basis_ay(:))

!Get change of basis from global to edge normal coordinate system 
call get_basis_change(Mb2a,Ma2b,basis_bx,basis_by,basis_ax,basis_ay)

!Change internal velocity to edge normal basis
vel_in_n = change_basis(Mb2a,vel_in)
Machin = norm2(vel_in_n)/cin

!Set downstream pressure
if (Machin .GE. 1.0d0) then 
    pinf = pin
else
    pinf = flowvars(zr)%backp_pset
end if 

!Apply boundary condition 
if (Machin .LT. 1.0d0) then !subsonic 

    !Freestream primative vector (in edge normal basis)
    rhoinf = pinf/(tin*flowvars(zr)%R)
    cinf = sqrt(flowvars(zr)%gam*(pinf/rhoinf))
    vel_inf_n(:) = vel_in_n(:)

    !If backflow then set the edge normal velocity to zero 
    if (vnorm_in .GE. 0.0d0) then !backflow
        vel_inf_n(1) = 0.0d0 
        ! vel_inf_n(1) = -vel_inf_n(1)
    end if 

    !Evaluate characteristic matrix
    ! call charmat(chm,chmi,rhoinf,cinf)
    call charmat(chm,chmi,rhoin,cin)

    !Evaluate flow deltas from freestream
    fdeltasi(1,1) = rhoin - rhoinf
    fdeltasi(2,1) = vel_in_n(1) - vel_inf_n(1) 
    fdeltasi(3,1) = vel_in_n(2) - vel_inf_n(2)  
    fdeltasi(4,1) = pin - pinf

    !Evaluate characteristics
    fchars = matmul(chm,fdeltasi)

    !Apply boundary condition
    if (vnorm_in .GE. 0.0d0) then !backflow
        fchars(1:3,1) = 0.0d0 
        ! fchars(1:3,1) = -0.1d0*fchars(1:3,1) 
    else !outflow 
        fchars(4,1) = 0.0d0
        ! fchars(4,1) = -0.1d0*fchars(4,1)
    end if 

    !Calculate boundary deltas
    fdeltasb = matmul(chmi,fchars)

    !Evaluate boundary velocity
    vel_b_n(1) = fdeltasb(2,1) + vel_inf_n(1) 
    vel_b_n(2) = fdeltasb(3,1) + vel_inf_n(2) 

    !If backflow then set the edge normal velocity to zero 
    if (vnorm_in .GE. 0.0d0) then !backflow
        vel_b_n(1) = 0.0d0 
        ! vel_b_n(2) = 0.0d0 
    end if 

    !Change boundary velocity to the global basis
    vel_b = change_basis(Ma2b,vel_b_n)
    ub = vel_b(1)
    vb = vel_b(2)

    !Evaluate boundary primatives 
    rhob = fdeltasb(1,1) + rhoinf
    if (options%force_fixed_pratio) then 
        pb = pinf
    else
        pb = fdeltasb(4,1) + pinf
    end if 
    eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
    ! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)
else !supersonic 
    ub = vel_in(1) 
    vb = vel_in(2)
    pb = pin
    rhob = rhoin 
    eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
    ! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)
end if 

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = rhob*ub*eb + pb*ub
! Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = rhob*vb*eb + pb*vb
! Gleft(4) = vb*(eb + pb)
return 
end subroutine backpressure_outflow_bc_c


end module flow_boundary_cond_mod