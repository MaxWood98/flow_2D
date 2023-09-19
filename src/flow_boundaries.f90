!Flow solver boundary conditions subroutines module
!Max Wood (mw16116@bristol.ac.uk)
!University of Bristol - Department of Aerospace Engineering 
!Version: 0.5.2
!Updated: 11/09/23

!Module
module flow_boundary_cond_mod
use flow_data_mod
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
! F(4) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%u(cellidx)*flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx)*&
! flowvars(zidx)%u(cellidx)
F(4) = flowvars(zidx)%u(cellidx)*(flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx))

G(1) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)
G(2) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)*flowvars(zidx)%u(cellidx)
G(3) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)*flowvars(zidx)%v(cellidx) + flowvars(zidx)%p(cellidx)
! G(4) = flowvars(zidx)%rho(cellidx)*flowvars(zidx)%v(cellidx)*flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx)*&
! flowvars(zidx)%v(cellidx)
G(4) = flowvars(zidx)%v(cellidx)*(flowvars(zidx)%e(cellidx) + flowvars(zidx)%p(cellidx))
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
! eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
! Fleft(4) = rhob*ub*eb + pb*ub
Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = rhob*vb*eb + pb*vb
Gleft(4) = vb*(eb + pb)
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
! eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
! Fleft(4) = rhob*ub*eb + pb*ub
Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = rhob*vb*eb + pb*vb
Gleft(4) = vb*(eb + pb)
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
! eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
! Fleft(4) = rhob*ub*eb + pb*ub
Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = rhob*vb*eb + pb*vb
Gleft(4) = vb*(eb + pb)
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
! eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
! Fleft(4) = rhob*ub*eb + pb*ub
Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = rhob*vb*eb + pb*vb
Gleft(4) = vb*(eb + pb)
return
end subroutine subsonic_outflow_bc_ri


!================================================================================
!Characteristic Farfield Conditions =============================================
!================================================================================


!Subsonic inflow/outflow boundary condition subroutine =========================
subroutine subsonic_bc_char(Fleft,Gleft,cr,zr,flowvars,if_of)
implicit none 

!Variables - Import
integer(in) :: cr,zr,if_of
real(dp), dimension(:) :: Fleft,Gleft
type(flow_var_data), dimension(:) :: flowvars

!Variables - Local
real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,cin,cinf
real(dp) :: ub,vb,rhob,pb,eb
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4)

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

!Evaluate flow deltas from freestream
fdeltasi(1,1) = rhoin - rhoinf
fdeltasi(2,1) = uin - uinf 
fdeltasi(3,1) = vin - vinf 
fdeltasi(4,1) = pin - pinf

!Evaluate characteristic matrix
call charmat(chm,chmi,rhoinf,cinf)
! call charmat(chm,chmi,0.5d0*(rhoinf+rhoin),0.5d0*(cinf+cin))

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

!Evaluate boundary velocity carrying though internal tangential component
ub = fdeltasb(2,1) + uinf
vb = fdeltasb(3,1) + vinf

!Evaluate boundary primatives 
rhob = fdeltasb(1,1) + rhoinf
pb = fdeltasb(4,1) + pinf
eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = vb*(eb + pb) 
return 
end subroutine subsonic_bc_char

! subroutine subsonic_bc_char(Fleft,Gleft,cr,zr,flowvars,mesh,edgidx,if_of)
! implicit none 

! !Variables - Import
! integer(in) :: cr,zr,if_of,edgidx
! real(dp), dimension(:) :: Fleft,Gleft
! type(mesh_data), dimension(:) :: mesh
! type(flow_var_data), dimension(:) :: flowvars

! !Variables - Local
! real(dp) :: nx,ny
! real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,vnorm_in,vnorm_inf,cin,cinf
! real(dp) :: unin,vnin,utin,vtin,uninf,vninf,utinf,vtinf
! real(dp) :: ub,vb,rhob,pb,eb
! real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
! real(dp) :: chm(4,4),chmi(4,4)

! !Cell normal vector with positive dot product with the velocity vector 
! if (if_of == 1) then !inflow
!     nx = mesh(zr)%edge_nx(edgidx)
!     ny = mesh(zr)%edge_ny(edgidx)
! elseif (if_of == -1) then !outflow
!     nx = -mesh(zr)%edge_nx(edgidx)
!     ny = -mesh(zr)%edge_ny(edgidx)
! end if 

! !Local primative vector
! rhoin = flowvars(zr)%rho(cr) 
! uin = flowvars(zr)%u(cr)
! vin = flowvars(zr)%v(cr) 
! pin = flowvars(zr)%p(cr)
! cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

! !Freestream primative vector
! rhoinf = flowvars(zr)%rhoinf
! uinf = flowvars(zr)%uinf
! vinf = flowvars(zr)%vinf
! pinf = flowvars(zr)%pinf
! cinf = flowvars(zr)%cinf

! !Surface normal and tangential velocity at internal and external 
! vnorm_in = uin*nx + vin*ny
! unin = nx*vnorm_in
! vnin = ny*vnorm_in
! utin = uin - unin
! vtin = vin - vnin

! vnorm_inf = uinf*nx + vinf*ny
! uninf = nx*vnorm_inf
! vninf = ny*vnorm_inf
! utinf = uinf - uninf
! vtinf = vinf - vninf

! ! if (if_of == 1) then !inflow

! !     vnorm_in = uin*nx + vin*ny
! !     unin = nx*vnorm_in
! !     vnin = ny*vnorm_in
! !     utin = uin - unin
! !     vtin = vin - vnin

! !     vnorm_inf = uinf*nx + vinf*ny
! !     uninf = nx*vnorm_inf
! !     vninf = ny*vnorm_inf
! !     utinf = uinf - uninf
! !     vtinf = vinf - vninf

! ! elseif (if_of == -1) then !outflow

! !     vnorm_in = uin*nx + vin*ny
! !     unin = nx*vnorm_in
! !     vnin = ny*vnorm_in
! !     utin = uin + unin
! !     vtin = vin + vnin

! !     vnorm_inf = uinf*nx + vinf*ny
! !     uninf = nx*vnorm_inf
! !     vninf = ny*vnorm_inf
! !     utinf = uinf + uninf
! !     vtinf = vinf + vninf

! ! end if 


! !Evaluate flow deltas from freestream
! ! fdeltasi(1,1) = rhoin - rhoinf
! ! fdeltasi(2,1) = unin - uninf 
! ! fdeltasi(3,1) = vnin - vninf 
! ! fdeltasi(4,1) = pin - pinf


! fdeltasi(1,1) = rhoin - rhoinf
! fdeltasi(2,1) = uin - uinf 
! fdeltasi(3,1) = vin - vinf 
! fdeltasi(4,1) = pin - pinf



! !Evaluate characteristic matrix
! call charmat(chm,chmi,rhoinf,cinf)
! ! call charmat(chm,chmi,0.5d0*(rhoinf+rhoin),0.5d0*(cinf+cin))

! !Evaluate characteristics
! fchars = matmul(chm,fdeltasi)

! !Apply boundary condition 
! if (if_of == 1) then !inflow
!     fchars(1:3,1) = 0.0d0 
! elseif (if_of == -1) then !outflow
!     fchars(4,1) = 0.0d0 
! end if 

! !Calculate boundary deltas
! fdeltasb = matmul(chmi,fchars)


! !Evaluate boundary velocity carrying though internal tangential component
! ! ub = fdeltasb(2,1) + uninf + utin
! ! vb = fdeltasb(3,1) + vninf + vtin

! ub = fdeltasb(2,1) + uinf !+ utin
! vb = fdeltasb(3,1) + vinf !+ vtin

! !Evaluate boundary primatives 
! rhob = fdeltasb(1,1) + rhoinf
! pb = fdeltasb(4,1) + pinf
! eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

! !Construct boundary fluxes
! Fleft(1) = rhob*ub
! Fleft(2) = rhob*ub*ub + pb
! Fleft(3) = rhob*ub*vb
! Fleft(4) = ub*(eb + pb)

! Gleft(1) = rhob*vb
! Gleft(2) = rhob*vb*ub
! Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = vb*(eb + pb) 
! return 
! end subroutine subsonic_bc_char




!Characteristic matrix and inverse subroutine =========================
subroutine charmat(chm,chmi,rho,c)
implicit none 

!Variables - Import 
real(dp) :: rho,c
real(dp) :: chm(4,4),chmi(4,4)

!Base matrix 
chm(:,:) = 0.0d0 
chm(1,1)= -c*c
chm(1,4)= 1.0d0
chm(2,3)= rho*c
chm(3,2)= rho*c
chm(3,4)= 1.0d0 
chm(4,2)= -rho*c
chm(4,4)= 1.0d0

!Inverse matrix 
chmi(:,:) = 0.0d0 
chmi(1,1)= -1.0d0/(c*c)
chmi(1,3)= 1.0d0/(2.0d0*c*c)
chmi(1,4)= 1.0d0/(2.0d0*c*c)
chmi(2,3)= 1.0d0/(2.0d0*rho*c)
chmi(2,4)= -1.0d0/(2.0d0*rho*c)
chmi(3,2)= 1.0d0/(rho*C)
chmi(4,3)= 0.5d0
chmi(4,4)= 0.5d0
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
real(dp) :: rhoin,uin,vin,pin,cin,H,velinN,velin,velinf
real(dp) :: Rm,a,b,c,cbp,cbm,cb,pwall
real(dp) :: Mb,pb,rhob,ub,vb,eb,uninf,vninf,tinf,rhoinf,pinf
real(dp) :: unin,vnin,utin,vtin,nx,ny
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4)

!Cell normal vector with positive dot product with the velocity vector 
nx = mesh(zr)%edge_nx(edge_idx)
ny = mesh(zr)%edge_ny(edge_idx)

!Internal primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

!Internal velocity and wall normal velocity
velin = sqrt(uin*uin + vin*vin)
velinN = uin*mesh(zr)%edge_nx(edge_idx) + vin*mesh(zr)%edge_ny(edge_idx)

!Internal boundary normal and tangential velocity
unin = mesh(zr)%edge_nx(edge_idx)*velinN
vnin = mesh(zr)%edge_ny(edge_idx)*velinN
utin = uin - unin
vtin = vin - vnin

!Enthalpy (from interior properties)
H = ((cin*cin)/(flowvars(zr)%gam - 1.0d0)) + 0.5d0*(uin*uin + vin*vin) 

!Reimann invariant - 
Rm = -velin + (2.0d0*cin/(flowvars(zr)%gam - 1.0d0))

!Solve for boundary speed of sound 
a = 1.0d0 + (2.0d0/(flowvars(zr)%gam - 1.0d0))
b = -2.0d0*Rm
c = 0.5d0*(flowvars(zr)%gam - 1.0d0)*(Rm*Rm - 2.0d0*H)
cbp = (-b + sqrt(b*b - 4.0d0*a*c))/(2.0d0*a)
cbm = (-b - sqrt(b*b - 4.0d0*a*c))/(2.0d0*a)
cb = max(cbp,cbm)

!Boundary normal velocity and mach number (from internal)
velinf = (2.0d0*cb/(flowvars(zr)%gam - 1.0d0)) - Rm
uninf = velinf*mesh(zr)%edge_nx(edge_idx)
vninf = velinf*mesh(zr)%edge_ny(edge_idx)
Mb = velinf/cb

!Boundary normal velocity and mach number (from external)
! velinf = sqrt(flowvars(zr)%uinf*flowvars(zr)%uinf + flowvars(zr)%vinf*flowvars(zr)%vinf)
! uninf = mesh(zr)%edge_nx(edge_idx)*velinf
! vninf = mesh(zr)%edge_ny(edge_idx)*velinf
! Mb = sqrt(uninf*uninf + vninf*vninf)/cb

!Cases
if (velinN .LT. 0.0d0) then !Backflow -> set as solid wall boundary 
    call wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
else !Normal inflow -> switch based on mach number
    if (Mb .GE. 1.0d0) then !Supersonic inflow
        call supersonic_inflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh) 
    else !Subsonic inflow
        !Boundary pressure temperature and density
        pinf = flowvars(zr)%p0inf/(((1.0d0 + 0.5d0*(flowvars(zr)%gam - 1.0d0)*Mb*Mb)**&
        (flowvars(zr)%gam/(flowvars(zr)%gam - 1.0d0))))
        tinf = flowvars(zr)%t0inf/(1.0d0 + 0.5d0*(flowvars(zr)%gam - 1.0d0)*Mb*Mb)
        rhoinf = pinf/(tinf*flowvars(zr)%R) 

        !Evaluate flow deltas from 'freestream'
        fdeltasi(1,1) = rhoin - rhoinf
        fdeltasi(2,1) = unin - uninf 
        fdeltasi(3,1) = vnin - vninf 
        fdeltasi(4,1) = pin - pinf

        !Evaluate characteristic matrix
        call charmat(chm,chmi,rhoinf,cb)

        !Evaluate characteristics
        fchars = matmul(chm,fdeltasi)

        !Apply inflow boundary condition 
        fchars(1:3,1) = 0.0d0 

        !Calculate boundary deltas
        fdeltasb = matmul(chmi,fchars)

        !Evaluate boundary velocity carrying though internal tangential component
        ub = fdeltasb(2,1) + uninf + utin
        vb = fdeltasb(3,1) + vninf + vtin

        !Evaluate boundary primatives 
        rhob = fdeltasb(1,1) + rhoinf
        pb = fdeltasb(4,1) + pinf
        eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

        !Accumulate inflow mass flux 
        flowvars(zr)%mflux_in = flowvars(zr)%mflux_in + rhob*(ub*nx + vb*ny)*mesh(zr)%edgelen(edge_idx)

        !Set boundary flux 
        Fleft(1) = rhob*ub
        Fleft(2) = rhob*ub*ub + pb
        Fleft(3) = rhob*ub*vb
        Fleft(4) = ub*(eb + pb)
        
        Gleft(1) = rhob*vb
        Gleft(2) = rhob*vb*ub
        Gleft(3) = rhob*vb*vb + pb
        Gleft(4) = vb*(eb + pb)
    end if
end if
return 
end subroutine stagnation_inflow_bc



! subroutine stagnation_inflow_bc(Fleft,Gleft,flowvars,mesh,edge_idx,cr,zr)
! implicit none 

! !Variables - Import
! integer(in) :: edge_idx,cr,zr
! real(dp), dimension(:) :: Fleft,Gleft 
! type(flow_var_data), dimension(:) :: flowvars
! type(mesh_data), dimension(:) :: mesh

! !Variables - Local
! real(dp) :: rhoin,uin,vin,pin,cin,H,velinN
! real(dp) :: Rm,a,b,c,cbp,cbm,cb,pwall
! real(dp) :: velb,Mb,Mnb,tb,pb,rhob,ub,vb,eb
! real(dp) :: unin,vnin,utin,vtin

! !Internal primative vector
! rhoin = flowvars(zr)%rho(cr) 
! uin = flowvars(zr)%u(cr)
! vin = flowvars(zr)%v(cr) 
! pin = flowvars(zr)%p(cr)
! cin = sqrt(flowvars(zr)%gam*(pin/rhoin))

! !Internal velocity
! velinN = uin*mesh(zr)%edge_nx(edge_idx) + vin*mesh(zr)%edge_ny(edge_idx)

! !Surface normal and tangential velocity at internal 
! unin = mesh(zr)%edge_nx(edge_idx)*velinN
! vnin = mesh(zr)%edge_ny(edge_idx)*velinN
! utin = uin - unin
! vtin = vin - vnin

! !Enthalpy (from interior properties)
! H = ((cin*cin)/(flowvars(zr)%gam - 1.0d0)) + 0.5d0*(uin*uin + vin*vin) 

! !Reimann invariant - 
! Rm = velinN - (2.0d0*cin/(flowvars(zr)%gam - 1.0d0))

! !Solve for boundary speed of sound 
! a = 1.0d0 + (2.0d0/(flowvars(zr)%gam - 1.0d0))
! b = 2.0d0*Rm
! c = 0.5d0*(flowvars(zr)%gam - 1.0d0)*(Rm*Rm - 2.0d0*H)
! cbp = (-b + sqrt(b*b - 4.0d0*a*c))/(2.0d0*a)
! cbm = (-b - sqrt(b*b - 4.0d0*a*c))/(2.0d0*a)
! cb = max(cbp,cbm)

! !Boundary velocity and mach number 
! velb = (2.0d0*cb/(flowvars(zr)%gam - 1.0d0)) + Rm
! Mnb = velb/cb
! ub = velb*mesh(zr)%edge_nx(edge_idx) + utin
! vb = velb*mesh(zr)%edge_ny(edge_idx) + vtin
! velb = sqrt(ub*ub + vb*vb)
! Mb = velb/cb

! !Cases
! if (velinN .LT. 0.0d0) then !Backflow -> set as solid wall boundary 
!     call wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
! else !Normal inflow -> switch based on mach number
!     if (Mnb .GE. 1.0d0) then !Supersonic inflow
!         call supersonic_inflow_bc_ri(Fleft,Gleft,cr,zr,edge_idx,flowvars,mesh) 
!     else !Subsonic inflow 
!         !Boundary pressure temperature and density
!         pb = flowvars(zr)%p0inf/(((1.0d0 + 0.5d0*(flowvars(zr)%gam - 1.0d0)*Mb*Mb)**(flowvars(zr)%gam/(flowvars(zr)%gam - 1.0d0))))
!         tb = flowvars(zr)%t0inf/(1.0d0 + 0.5d0*(flowvars(zr)%gam - 1.0d0)*Mb*Mb)
!         rhob = pb/(tb*flowvars(zr)%R) 

!         !Boundary energy
!         eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*velb*velb

!         !Set boundary flux 
!         Fleft(1) = rhob*ub
!         Fleft(2) = rhob*ub*ub + pb
!         Fleft(3) = rhob*ub*vb
!         ! Fleft(4) = rhob*ub*eb + pb*ub
!         Fleft(4) = ub*(eb + pb)
        
!         Gleft(1) = rhob*vb
!         Gleft(2) = rhob*vb*ub
!         Gleft(3) = rhob*vb*vb + pb
!         ! Gleft(4) = rhob*vb*eb + pb*vb
!         Gleft(4) = vb*(eb + pb)
!     end if
! end if
! return 
! end subroutine stagnation_inflow_bc





!Freestream inflow boundary condition subroutine =========================
subroutine freestream_inflow_c(Fleft,Gleft,cr,zr,flowvars,mesh,options,edgidx)
implicit none 

!Variables - Import
integer(in) :: cr,zr,edgidx
real(dp), dimension(:) :: Fleft,Gleft
type(mesh_data), dimension(:) :: mesh
type(flow_var_data), dimension(:) :: flowvars
type(options_data) :: options 

!Variables - Local
real(dp) :: nx,ny,pwall
real(dp) :: uin,vin,rhoin,pin,uinf,vinf,rhoinf,pinf,vnorm_in,vnorm_inf,cin,cinf
real(dp) :: unin,vnin,utin,vtin,uninf,vninf,utinf,vtinf
real(dp) :: ub,vb,rhob,pb,eb
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4)

!Cell normal vector with positive dot product with the velocity vector 
nx = mesh(zr)%edge_nx(edgidx)
ny = mesh(zr)%edge_ny(edgidx)

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

!Surface normal and tangential velocity at internal and external 
vnorm_in = uin*nx + vin*ny
unin = nx*vnorm_in
vnin = ny*vnorm_in
utin = uin - unin
vtin = vin - vnin
vnorm_inf = uinf*nx + vinf*ny
uninf = nx*vnorm_inf
vninf = ny*vnorm_inf
utinf = uinf - uninf
vtinf = vinf - vninf

!Switch on sign of normal velocity
if (vnorm_in .LT. 0.0d0) then !backflow so default to wall condition
    call wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
else
    !Evaluate flow deltas from freestream
    fdeltasi(1,1) = rhoin - rhoinf
    fdeltasi(2,1) = unin - uninf 
    fdeltasi(3,1) = vnin - vninf 
    fdeltasi(4,1) = pin - pinf

    !Evaluate characteristic matrix
    call charmat(chm,chmi,rhoinf,cinf)
    ! call charmat(chm,chmi,0.5d0*(rhoinf+rhoin),0.5d0*(cinf+cin))

    !Evaluate characteristics
    fchars = matmul(chm,fdeltasi)

    !Apply inflow boundary condition 
    fchars(1:3,1) = 0.0d0 

    !Calculate boundary deltas
    fdeltasb = matmul(chmi,fchars)

    !Evaluate boundary velocity carrying though internal tangential component
    ub = fdeltasb(2,1) + uninf !+ utin
    vb = fdeltasb(3,1) + vninf !+ vtin

    !Evaluate boundary primatives 
    rhob = fdeltasb(1,1) + rhoinf
    if (options%force_fixed_pratio) then 
        pb = pinf
    else
        pb = fdeltasb(4,1) + pinf
    end if 
    eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

    !Accumulate inflow mass flux 
    flowvars(zr)%mflux_in = flowvars(zr)%mflux_in + rhob*(ub*nx + vb*ny)*mesh(zr)%edgelen(edgidx)

    !Construct boundary fluxes
    Fleft(1) = rhob*ub
    Fleft(2) = rhob*ub*ub + pb
    Fleft(3) = rhob*ub*vb
    Fleft(4) = ub*(eb + pb)

    Gleft(1) = rhob*vb
    Gleft(2) = rhob*vb*ub
    Gleft(3) = rhob*vb*vb + pb
    Gleft(4) = vb*(eb + pb) 
end if 
return 
end subroutine freestream_inflow_c



!Supersonic internal inflow boundary condition subroutine =========================
! subroutine supersonic_internal_inflow(Fleft,Gleft,flowvars,mesh,edge_idx)
! implicit none 

! !Variables - Import
! integer(in) :: edge_idx
! real(dp), dimension(:) :: Fleft,Gleft
! type(flow_var_data) :: flowvars
! type(mesh_data) :: mesh

! !Variables - Local
! real(dp) :: ub,vb,rhob,pb,eb

! !Set boundary properties from the external values
! ub = flowvars%uinf*mesh%edge_nx(edge_idx) 
! vb = flowvars%vinf*mesh%edge_ny(edge_idx) 
! rhob = flowvars%rhoinf
! pb = flowvars%pinf
! eb = (pb/((flowvars%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)

! !Construct boundary fluxes
! Fleft(1) = rhob*ub
! Fleft(2) = rhob*ub*ub + pb
! Fleft(3) = rhob*ub*vb
! Fleft(4) = rhob*ub*eb + pb*ub

! Gleft(1) = rhob*vb
! Gleft(2) = rhob*vb*ub
! Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = rhob*vb*eb + pb*vb
! return 
! end subroutine supersonic_internal_inflow




!Supersonic internal inflow boundary condition subroutine =========================
! subroutine supersonic_internal_outflow(Fleft,Gleft,flowvars,cr)
! implicit none 

! !Variables - Import
! integer(in) :: cr
! real(dp), dimension(:) :: Fleft,Gleft
! type(flow_var_data) :: flowvars

! !Variables - Local
! real(dp) :: ub,vb,rhob,pb,eb

! !Set boundary properties from the internal values
! ub = flowvars%u(cr) 
! vb = flowvars%v(cr)
! rhob = flowvars%rho(cr)
! pb = flowvars%p(cr)
! eb = (pb/((flowvars%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)

! !Construct boundary fluxes
! Fleft(1) = rhob*ub
! Fleft(2) = rhob*ub*ub + pb
! Fleft(3) = rhob*ub*vb
! Fleft(4) = rhob*ub*eb + pb*ub

! Gleft(1) = rhob*vb
! Gleft(2) = rhob*vb*ub
! Gleft(3) = rhob*vb*vb + pb
! Gleft(4) = rhob*vb*eb + pb*vb
! return 
! end subroutine supersonic_internal_outflow




!Mass flux outflow boundary condition (pressure) subroutine =========================
! subroutine mass_flux_outflow_bc_P(R,mesh,flowvars,options,iter)
! implicit none 

! !Variables - Import 
! integer(in) :: iter
! real(dp), dimension(:,:) :: R
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars
! type(options_data) :: options

! !Variables - Local 
! integer(in) :: ii,cl,cr
! real(dp) :: nx,ny,rhoin,uin,vin,pin
! real(dp) :: Aout,Fm,Fp,mflux,dm
! real(dp) :: pb,ub,vb,tb,rhob,eb,pset_tgt
! real(dp) :: Fleft(4),Fright(4),Gleft(4),Gright(4)
! real(dp) :: Fedge(4),Gedge(4),Flxe(4)

! !Evaluate outflow mass flow rate and properties
! if ((mod(iter,options%Mflux_relaxiter) == 0) .OR. (iter == 1)) then 
!     Aout = 0.0d0 
!     Fm = 0.0d0 
!     Fp = 0.0d0 
!     mflux = 0.0d0  
!     do ii=1,mesh%nedge

!         !Left and right cells for this edge
!         cl = mesh%cell_left(ii)
!         cr = mesh%cell_right(ii)

!         !If outflow face
!         if (cl == -4) then 

!             !Edge unit normal vector
!             nx = -mesh%edge_nx(ii)
!             ny = -mesh%edge_ny(ii)

!             !Local primative vector
!             rhoin = flowvars%rho(cr) 
!             uin = flowvars%u(cr)
!             vin = flowvars%v(cr) 
!             pin = flowvars%p(cr)

!             !Accumulate boundary integrals 
!             Fm = Fm + rhoin*((uin*nx + vin*ny)**2)*mesh%edgelen(ii)
!             Fp = Fp + pin*mesh%edgelen(ii)
!             mflux = mflux + rhoin*(uin*nx + vin*ny)*mesh%edgelen(ii)

!             !Accumulate outflow area
!             Aout = Aout + mesh%edgelen(ii)
!         end if 
!     end do 

!     !Set required downstream static pressure
!     if (mflux .NE. 0.0d0) then
!         dm = (1.0d0 - (flowvars%mflow_in/mflux))
!     else
!         dm = 0.0d0  
!     end if
!     pset_tgt = (1.0d0/Aout)*(Fm*dm + Fp)
!     if (iter == 1) then 
!         flowvars%mflux_pset0 = pset_tgt
!     end if 
!     flowvars%mflux_pset = (pset_tgt - flowvars%mflux_pset0)*options%Mflux_relax + flowvars%mflux_pset0
! end if 

! !Set boundary pressure 
! pb = flowvars%mflux_pset

! !Set downstream properties
! do ii=1,mesh%nedge

!     !Left and right cells for this edge
!     cl = mesh%cell_left(ii)
!     cr = mesh%cell_right(ii)

!     !If outflow face
!     if (cl == -4) then 

!         !Local primative vector
!         rhoin = flowvars%rho(cr) 
!         pin = flowvars%p(cr)

!         !Set boundary properties
!         ub = flowvars%u(cr)
!         vb = flowvars%v(cr) 
!         tb = pin/(rhoin*flowvars%R)
!         rhob = pb/(tb*flowvars%R)
!         eb = (pb/((options%gamma - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)

!         !Construct boundary fluxes 
!         Fleft(1) = rhob*ub
!         Fleft(2) = rhob*ub*ub + pb
!         Fleft(3) = rhob*ub*vb
!         Fleft(4) = rhob*ub*eb + pb*ub

!         Gleft(1) = rhob*vb
!         Gleft(2) = rhob*vb*ub
!         Gleft(3) = rhob*vb*vb + pb
!         Gleft(4) = rhob*vb*eb + pb*vb

!         !Cell right flux
!         call cell_flux(Fright,Gright,flowvars,cr)

!         !Edge flux value ([dy,-dx] normal)
!         Fedge(:) = 0.5d0*(Fleft(:) + Fright(:))*mesh%edge_dy(ii)
!         Gedge(:) = 0.5d0*(Gleft(:) + Gright(:))*mesh%edge_dx(ii)
!         Flxe(:) =  Fedge(:) - Gedge(:)

!         !Accumulate to cell right 
!         R(cr,:) = R(cr,:) - Flxe(:)
!     end if
! end do 
! return 
! end subroutine mass_flux_outflow_bc_P




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
            ! eb = (pb/((options%gamma - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
            eb = (pb/(options%gamma - 1.0d0)) + 0.5d0*rhob*(ub**2 + vb**2)

            !Construct boundary fluxes 
            Fleft(1) = rhob*ub
            Fleft(2) = rhob*ub*ub + pb
            Fleft(3) = rhob*ub*vb
            ! Fleft(4) = rhob*ub*eb + pb*ub
            Fleft(4) = ub*(eb + pb)

            Gleft(1) = rhob*vb
            Gleft(2) = rhob*vb*ub
            Gleft(3) = rhob*vb*vb + pb
            ! Gleft(4) = rhob*vb*eb + pb*vb
            Gleft(4) = vb*(eb + pb)
        elseif ((M_in .LT. 1.0d0) .AND. (M_b .GE. 1.0d0)) then !internal subsonic / boundary supersonic -> set boundary to sonic to choke 

            !Force boundary normal velocity to sonic 
            velb = cb

            !Set boundary properties
            rhob = ((cb*cb)/(flowvars(zone)%gam*Sb))**(1.0d0/(flowvars(zone)%gam - 1.0d0))
            pb = (rhob*cb*cb)/flowvars(zone)%gam 
            ub = velb*nx + ut 
            vb = velb*ny + vt
            ! eb = (pb/((options%gamma - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
            eb = (pb/(options%gamma - 1.0d0)) + 0.5d0*rhob*(ub**2 + vb**2)

            !Construct boundary fluxes 
            Fleft(1) = rhob*ub
            Fleft(2) = rhob*ub*ub + pb
            Fleft(3) = rhob*ub*vb
            ! Fleft(4) = rhob*ub*eb + pb*ub
            Fleft(4) = ub*(eb + pb)

            Gleft(1) = rhob*vb
            Gleft(2) = rhob*vb*ub
            Gleft(3) = rhob*vb*vb + pb
            ! Gleft(4) = rhob*vb*eb + pb*vb
            Gleft(4) = vb*(eb + pb)
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





!Mass flux inflow boundary condition subroutine =========================
! subroutine mass_flux_inflow_bc(R,mesh,flowvars)
! implicit none 

! !Variables - Import 
! real(dp), dimension(:,:) :: R
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars
    
! !Variables - Local 
! integer(in) :: ii,cl,cr
! real(dp) :: nx,ny,rhoin
! real(dp) :: pb,ub,vb,tb,rhob,eb,Mb,velb
! real(dp) :: Fleft(4),Fright(4),Gleft(4),Gright(4)
! real(dp) :: Fedge(4),Gedge(4),Flxe(4)

! !Find inflow average density
! rhob = 0.0d0 
! do ii=1,mesh%nedge

!     !Left and right cells for this edge
!     cl = mesh%cell_left(ii)
!     cr = mesh%cell_right(ii)

!     !If inflow face
!     if (cl == -3) then 
!         rhoin = flowvars%rho(cr) 
!         rhob = rhob + rhoin*mesh%edgelen(ii)
!     end if
! end do
! rhob = rhob/flowvars%A_in

! !Boundary velocity and mach number
! velb = (flowvars%mflow_in/(rhob*flowvars%A_in))
! Mb = velb/flowvars%cinf

! !Apply conditon to inflow faces
! do ii=1,mesh%nedge

!     !Left and right cells for this edge
!     cl = mesh%cell_left(ii)
!     cr = mesh%cell_right(ii)

!     !If inflow face
!     if (cl == -3) then 

!         !Face unit normal vector
!         nx = mesh%edge_nx(ii)
!         ny = mesh%edge_ny(ii)

!         !Conditions 
!         ub = velb*nx
!         vb = velb*ny
!         if (Mb .GE. 1.0d0) then 
!             rhob = flowvars%rhoinf
!             pb = flowvars%pinf
!         else
!             tb = flowvars%t0inf - ((flowvars%gam - 1.0d0)/(2.0d0*flowvars%gam*flowvars%R))*velb*velb
!             rhob = flowvars%rho(cr) 
!             pb = rhob*flowvars%R*tb
!         end if
!         eb = (pb/((flowvars%gam - 1.0d0)*rhob)) + 0.5d0*velb*velb

!         !Construct boundary fluxes
!         Fleft(1) = rhob*ub
!         Fleft(2) = rhob*ub*ub + pb
!         Fleft(3) = rhob*ub*vb
!         Fleft(4) = rhob*ub*eb + pb*ub

!         Gleft(1) = rhob*vb
!         Gleft(2) = rhob*vb*ub
!         Gleft(3) = rhob*vb*vb + pb
!         Gleft(4) = rhob*vb*eb + pb*vb

!         !Cell right flux
!         call cell_flux(Fright,Gright,flowvars,cr)

!         !Edge flux value ([dy,-dx] normal)
!         Fedge(:) = 0.5d0*(Fleft(:) + Fright(:))*mesh%edge_dy(ii)
!         Gedge(:) = 0.5d0*(Gleft(:) + Gright(:))*mesh%edge_dx(ii)
!         Flxe(:) =  Fedge(:) - Gedge(:)

!         !Accumulate to cell right 
!         R(cr,:) = R(cr,:) - Flxe(:)
!     end if 
! end do 
! return 
! end subroutine mass_flux_inflow_bc




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




!Subroutine to find the inflow mass flow rate =========================
! subroutine find_intake_massflow(flowvars,mesh)
! implicit none 

! !Variables - Import 
! type(mesh_data) :: mesh
! type(flow_var_data) :: flowvars

! !Variables - Local 
! integer(in) :: ii,cl,cr
! real(dp) :: nx,ny,uin,vin,rhoin,A_rho_out

! !Evaluate inflow mass flow rate and properties
! A_rho_out = 0.0d0 
! flowvars%A_in = 0.0d0 
! flowvars%mflow_in = 0.0d0 
! do ii=1,mesh%nedge

!     !Left and right cells for this edge
!     cl = mesh%cell_left(ii)
!     cr = mesh%cell_right(ii)

!     !If inflow face
!     if ((cl == -3) .OR. (cl == -5)) then !If inflow face

!         !Cell unit normal vector (into domain)
!         nx = mesh%edge_nx(ii)
!         ny = mesh%edge_ny(ii)

!         !Local primative vector
!         rhoin = flowvars%rhoinf !flowvars%rho(cr) 
!         uin = flowvars%uinf !flowvars%u(cr)
!         vin = flowvars%vinf !flowvars%v(cr) 

!         !Accumulate inflow mass flow and area
!         flowvars%mflow_in  = flowvars%mflow_in  + rhoin*(uin*nx + vin*ny)*mesh%edgelen(ii)
!         flowvars%A_in = flowvars%A_in + mesh%edgelen(ii)
!     end if 

!     !If outflow face
!     if (cl == -4) then 
!         A_rho_out = A_rho_out + mesh%edgelen(ii)*flowvars%rho(cr)
!     end if 
! end do 

! !Set outflow velocity for equal outflow mass flux
! flowvars%vnset_of = flowvars%mflow_in/A_rho_out
! return 
! end subroutine find_intake_massflow




!Outflow back pressure boundary condition (Riemann) =========================
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
eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!Construct boundary fluxes
Fleft(1) = rhob*ub
Fleft(2) = rhob*ub*ub + pb
Fleft(3) = rhob*ub*vb
Fleft(4) = ub*(eb + pb)

Gleft(1) = rhob*vb
Gleft(2) = rhob*vb*ub
Gleft(3) = rhob*vb*vb + pb
Gleft(4) = vb*(eb + pb) 
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
real(dp) :: rhoin,uin,vin,pin,cin,tin,Machin,velin
real(dp) :: pb,rhob,ub,vb,eb,pwall,nx,ny

real(dp) :: uinf,vinf,rhoinf,pinf,cinf,vnorm_in,vnorm_inf
real(dp) :: unin,vnin,utin,vtin,uninf,vninf,utinf,vtinf
real(dp) :: fdeltasi(4,1),fdeltasb(4,1),fchars(4,1)
real(dp) :: chm(4,4),chmi(4,4)

!Outflow zone of this edge
zof = mesh(zr)%outflow_zone(edge_idx)

!Internal primative vector
rhoin = flowvars(zr)%rho(cr) 
uin = flowvars(zr)%u(cr)
vin = flowvars(zr)%v(cr) 
pin = flowvars(zr)%p(cr)
velin = sqrt(uin*uin + vin*vin)
cin = sqrt(flowvars(zr)%gam*(pin/rhoin))
Machin = velin/cin
tin = pin/(rhoin*flowvars(zr)%R)

!Set downstream pressure
if (Machin .GE. 1.0d0) then 
    pinf = pin
else
    pinf = flowvars(zr)%backp_pset
end if 

!Set condition 
if ((uin*mesh(zr)%edge_nx(edge_idx) + vin*mesh(zr)%edge_ny(edge_idx)) .GT. 0.0d0) then !inflow -> wall condition
    call wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
else !outflow -> pressure condition 

    !Edge normal 
    nx = -mesh(zr)%edge_nx(edge_idx)
    ny = -mesh(zr)%edge_ny(edge_idx)

    ! !Boundary primatives (average outflow)
    ! uinf = flowvars(zr)%u_outflow(zof)
    ! vinf = flowvars(zr)%v_outflow(zof)
    ! rhoinf = flowvars(zr)%rho_outflow(zof)
    ! cinf = sqrt(flowvars(zr)%gam*(pinf/rhoinf))

    !Boundary primatives (local)
    uinf = flowvars(zr)%u(cr)
    vinf = flowvars(zr)%v(cr)
    rhoinf = flowvars(zr)%rho(cr)
    cinf = sqrt(flowvars(zr)%gam*(pinf/rhoinf))

    !Internal and external boundary normal and tangential velocity
    vnorm_in = uin*nx + vin*ny
    unin = nx*vnorm_in
    vnin = ny*vnorm_in
    utin = uin - unin
    vtin = vin - vnin
    vnorm_inf = uinf*nx + vinf*ny
    uninf = nx*vnorm_inf
    vninf = ny*vnorm_inf
    utinf = uinf - uninf
    vtinf = vinf - vninf

    !Toggle subsonic/supersonic condition 
    if (Machin .LT. 1.0d0) then !subsonic 

        !Evaluate flow deltas from freestream
        fdeltasi(1,1) = rhoin - rhoinf
        fdeltasi(2,1) = unin - uninf 
        fdeltasi(3,1) = vnin - vninf 
        fdeltasi(4,1) = pin - pinf

        !Evaluate characteristic matrix
        call charmat(chm,chmi,rhoinf,cinf)

        !Evaluate characteristics
        fchars = matmul(chm,fdeltasi)

        !Apply outflow boundary condition 
        fchars(4,1) = 0.0d0 

        !Calculate boundary deltas
        fdeltasb = matmul(chmi,fchars)

        !Evaluate boundary velocity
        ub = fdeltasb(2,1) + uninf + utin
        vb = fdeltasb(3,1) + vninf + vtin

        !Evaluate boundary primatives 
        rhob = fdeltasb(1,1) + rhoin
        if (options%force_fixed_pratio) then 
            pb = pinf
        else
            pb = fdeltasb(4,1) + pinf
        end if 
        eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)
    else !supersonic 
        ub = uin 
        vb = vin 
        pb = pin
        rhob = rhoin 
        eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)
    end if 

    !Construct boundary fluxes
    Fleft(1) = rhob*ub
    Fleft(2) = rhob*ub*ub + pb
    Fleft(3) = rhob*ub*vb
    Fleft(4) = ub*(eb + pb)

    Gleft(1) = rhob*vb
    Gleft(2) = rhob*vb*ub
    Gleft(3) = rhob*vb*vb + pb
    Gleft(4) = vb*(eb + pb) 
end if 
return 
end subroutine backpressure_outflow_bc_c




! subroutine backpressure_outflow_bc(Fleft,Gleft,flowvars,mesh,edge_idx,cr,zr)
! implicit none 

! !Variables - Import
! integer(in) :: cr,zr,edge_idx
! real(dp), dimension(:) :: Fleft,Gleft 
! type(flow_var_data), dimension(:) :: flowvars
! type(mesh_data), dimension(:) :: mesh

! !Variables - Local
! real(dp) :: rhoin,uin,vin,pin,cin,tin,Machin,velin
! real(dp) :: pb,rhob,ub,vb,eb,pwall!,pshock,Machshock

! !Internal primative vector
! rhoin = flowvars(zr)%rho(cr) 
! uin = flowvars(zr)%u(cr)
! vin = flowvars(zr)%v(cr) 
! pin = flowvars(zr)%p(cr)
! velin = sqrt(uin*uin + vin*vin)
! cin = sqrt(flowvars(zr)%gam*(pin/rhoin))
! Machin = velin/cin
! tin = pin/(rhoin*flowvars(zr)%R)

! !Set downstream pressure
! if (Machin .GE. 1.0d0) then 
!     pb = pin
! else
!     pb = flowvars(zr)%backp_pset
! end if 

! !Set condition 
! if ((uin*mesh(zr)%edge_nx(edge_idx) + vin*mesh(zr)%edge_ny(edge_idx)) .GT. 0.0d0) then !inflow -> wall condition
!     call wall_flux_bc_Cextrp(pwall,Fleft,Gleft,flowvars,cr,zr)
! else !normal

!     !Check shock condition if supersonic and then set boundary properties
!     ! if (Machin .GE. 1.0d0) then 
!     !     pshock = pin*((2.0d0*flowvars%gam*Machin*Machin - (flowvars%gam - 1.0d0))/(flowvars%gam + 1.0d0))
!     !     if (pshock .LT. flowvars%backp_pset) then !assume shock exists at the boundary -> use post shock state for flux evaluation 
!     !         Machshock = (2.0d0 + (flowvars%gam - 1.0d0)*Machin*Machin)/(2.0d0*flowvars%gam*Machin*Machin - (flowvars%gam - 1.0d0))
!     !         pb = flowvars%backp_pset!pshock 
!     !         rhob = rhoin*((flowvars%gam + 1.0d0)*Machin*Machin)/(2.0d0 + (flowvars%gam - 1.0d0)*Machin*Machin)
!     !         ub = Machshock*cin
!     !         vb = 0.0d0 !Machshock*cin
!     !         eb = (pb/((flowvars%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
!     !     else
!     !         ub = flowvars%u(cr)
!     !         vb = flowvars%v(cr)
!     !         rhob = pb/(tin*flowvars%R)
!     !         eb = (pb/((flowvars%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
!     !     end if 
!     ! else
!         ! ub = flowvars%u(cr)
!         ! vb = flowvars%v(cr)
!         ! rhob = pb/(tin*flowvars%R)
!         ! eb = (pb/((flowvars%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
!     ! end if 

!     !Set boundary properties
!     ub = flowvars(zr)%u(cr)
!     vb = flowvars(zr)%v(cr)
!     rhob = pb/(tin*flowvars(zr)%R)
!     ! eb = (pb/((flowvars(zr)%gam - 1.0d0)*rhob)) + 0.5d0*(ub**2 + vb**2)
!     eb = (pb/((flowvars(zr)%gam - 1.0d0))) + 0.5d0*rhob*(ub**2 + vb**2)

!     !Construct boundary fluxes 
!     Fleft(1) = rhob*ub
!     Fleft(2) = rhob*ub*ub + pb
!     Fleft(3) = rhob*ub*vb
!     ! Fleft(4) = rhob*ub*eb + pb*ub
!     Fleft(4) = ub*(eb + pb)

!     Gleft(1) = rhob*vb
!     Gleft(2) = rhob*vb*ub
!     Gleft(3) = rhob*vb*vb + pb
!     ! Gleft(4) = rhob*vb*eb + pb*vb
!     Gleft(4) = vb*(eb + pb)
! end if 
! return 
! end subroutine backpressure_outflow_bc


end module flow_boundary_cond_mod