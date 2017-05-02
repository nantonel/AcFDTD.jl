
using AcFDTD
using Base.Test

srand(123)
geo = CuboidRoom(7,8,10,abs(10*randn(6)), FDTDEnv(0.1,IWB()))
f = FDTD(geo) 
Nxyz = geo.Nx*geo.Ny*geo.Nz

Nt = 40
xs = [(6, 7, 9), (3, 3, 3), (2, 5, 2)]
xr = [(2, 2, 4), (2, 7, 8)]

B = get_B(geo,Nt)

Exp = get_Exp(geo,Nt,xs)
Sel = get_Sel(geo,Nt,xr)

s = randn(Nt,length(xs))
p = randn(Nt,length(xr))

S = randn(size(geo)...,Nt)
P = randn(size(geo)...,Nt)

println("testing Adjoint with receiver no IC")
SBE = Sel*(B\full(Exp))
p1 =  SBE*s[:]
s1 =  SBE'*p[:]

p2 = fdtd(xr,xs,Nt,f,s)
s2 = fdtdAdj(xr,xs,Nt,f,p)

println("testing calls")
@time s2 = fdtdAdj(xr,xs,Nt, f,p)

s3 = 0.*copy(s2)
p3 = 0.*copy(p2)

fdtd!(p3,xr,xs,Nt,f,s)
fdtdAdj!(s3,xr,xs,Nt,f,p)
fdtd!(p3,xr,xs,Nt,f,s)
println("testing in-place call")
@time fdtdAdj!(s3,xr,xs,Nt,f,p)

@test abs(vecdot(p[:],p1)-vecdot(s1,s[:]))<1e-10
@test abs(vecdot(p,p2)-vecdot(s2,s))<1e-10
@test abs(vecdot(p,p3)-vecdot(s3,s))<1e-10

println("testing Adjoint full no IC")
BE = (B\full(Exp))
P1 =  (BE*s[:])[2*prod(size(geo))+1:end]
s1 =  BE'*[zeros(2*prod(size(geo)));P[:]]

P2 = fdtd(xs,Nt,f,s)
s2 = fdtdAdj(xs,Nt,f,P)

println("testing calls")
@time s2 = fdtdAdj(xs,Nt,f,P)

s3 = 0.*copy(s2)
P3 = 0.*copy(P2)

fdtd!(P3,xs,Nt,f,s)
fdtdAdj!(s3,xs,Nt,f,P)
println("testing in-place call")
@time fdtdAdj!(s3,xs,Nt,f,P)

@test abs(vecdot(P[:],P1)-vecdot(s1,s[:]))<1e-10
@test abs(vecdot(P,P2)-vecdot(s2,s))<1e-10
@test abs(vecdot(P,P3)-vecdot(s3,s))<1e-10

println("testing Adjoint full, source everywhere no IC")

P2 = fdtd(Nt,f,S)
S2 = fdtdAdj(Nt,f,P)

println("testing calls")
@time S2 = fdtdAdj(Nt,f,P)

S3 = 0.*copy(S2)
P3 = 0.*copy(P2)

fdtd!(P3,Nt,f,S)
fdtdAdj!(S3,Nt,f,P)
println("testing in-place call")
@time fdtdAdj!(S3,Nt,f,P)

#@test abs(vecdot(P[:],P1)-vecdot(S1,S[:]))<1e-12
@test abs(vecdot(P,P2)-vecdot(S2,S))<1e-10
@test abs(vecdot(P,P3)-vecdot(S3,S))<1e-10

##TODO make adjoint with ICs
