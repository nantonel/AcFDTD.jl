
using AcFDTD
using Base.Test

srand(123)
geo = CuboidRoom(7,8,10,abs(10*randn(6)), FDTDEnv(0.1,SLF()))
Nxyz = geo.Nx*geo.Ny*geo.Nz

Nt = 20
xs = [6 7 9; 3 3 3]'
xr = [9 2 4; 8 7 8; 2 5 2]'

B = get_B(geo,Nt)

E = get_Exp(geo,Nt,xs)
S = get_Sel(geo,Nt,xr)

s = randn(Nt,size(xs,2))
p = randn(Nt,size(xr,2))

SBE = S*(B\full(E))
p1 =  SBE*s[:]
s1 =  SBE'*p[:]

p2 = fdtd(s,xs,xr,Nt,geo)
s2 = fdtdAdj(p,xs,xr,Nt,geo)

@test abs(vecdot(p[:],p1)-vecdot(s1,s[:]))<1e-12
@test abs(vecdot(p,p2)-vecdot(s2,s))<1e-12

Eic = [sparse(1:Nxyz,1:Nxyz,ones(Nxyz),(2+Nt)*Nxyz,2*Nxyz) E]
SBE = S*(B\full(Eic))

pic0,pic1 = randn(Nxyz),randn(Nxyz)
sic = [pic0;pic1;s[:]]
p1 =  SBE*sic
s1 =  SBE'*p[:]

@test abs(vecdot(p,p1)-vecdot(s1,sic))<1e-12
#TODO make adjoint with ICs
