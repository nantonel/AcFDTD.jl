using AcFDTD
using Base.Test
using Base.Profile
srand(123)

geo = CuboidRoom(10,15,16,[10.;20.;30.;40.;50.;60.], FDTDEnv(0.1,IWB()))
Nt = 200
s_t = [randn(Nt,1) randn(Nt,1)]
S = randn(size(geo)...,Nt)
xs = [(5, 5, 5),(5,2,4)]
xr = [(5, 8, 7),(2,2,2)]
f = FDTD(geo) 

p0ic = randn(10,15,16)
p1ic = randn(10,15,16)

println("test calls: \n")
println("----------      no IC")
p    = fdtd(xr, xs, Nt,geo, s_t)
p    = fdtd(xr, xs, Nt,  f, s_t)
p2 = copy(p)
fdtd!(p2, xr, xs, Nt,  f, s_t)
@time p    = fdtd(xr, xs, Nt,  f, s_t)
println("                       in-place call")
@time        fdtd!(p2,xr, xs, Nt,  f, s_t)

@test norm(p-p2) < 1e-8

println("----------         IC")
pic  = fdtd(xr, xs, Nt,geo, s_t, p0ic, p1ic)
ppic = fdtd(xr, xs, Nt,  f, s_t, p0ic, p1ic)
pic2 = copy(pic)
fdtd!(pic2, xr, xs, Nt,  f, s_t, p0ic, p1ic)
@time ppic = fdtd(xr, xs, Nt,  f, s_t, p0ic, p1ic)
println("                       in-place call")
@time        fdtd!(pic2, xr, xs, Nt,  f, s_t, p0ic, p1ic)

@test norm(pic-pic2) < 1e-8
@test norm(pic-ppic) < 1e-8

println("---------- full no IC")
P    = fdtd(    xs, Nt,geo, s_t)
P    = fdtd(    xs, Nt,  f, s_t)
P2 = copy(P)
fdtd!(P2,    xs, Nt,  f, s_t)
@time P    = fdtd(    xs, Nt,  f, s_t)
println("                       in-place call")
@time        fdtd!(P2, xs, Nt,  f, s_t)

@test vecnorm(P-P2) < 1e-8

println("---------- full    IC")
PIC  = fdtd(    xs, Nt,geo, s_t, p0ic, p1ic)
PIC  = fdtd(    xs, Nt,  f, s_t, p0ic, p1ic)
PIC2 = copy(PIC)
fdtd!(PIC2,    xs, Nt,  f, s_t, p0ic, p1ic)
@time PIC  = fdtd(    xs, Nt,  f, s_t, p0ic, p1ic)
println("                       in-place call")
@time      fdtd!(PIC2,xs, Nt,  f, s_t, p0ic, p1ic)

@test vecnorm(PIC-PIC2) < 1e-8

println("---------- source everywhere, no IC")
ps    = fdtd(    xr, Nt,geo, S)
ps    = fdtd(    xr, Nt,  f, S)
ps2 = copy(ps)
fdtd!(ps2, xr, Nt,  f, S)
@time ps    = fdtd(      xr, Nt,  f, S)
println("                       in-place call")
@time        fdtd!(ps2,  xr, Nt,  f, S)

println("---------- source everywhere, IC")
psic    = fdtd(    xr, Nt,geo, S, p0ic, p1ic)
psic    = fdtd(    xr, Nt,  f, S, p0ic, p1ic)
psic2 = copy(ps)
fdtd!(psic2, xr, Nt,  f, S, p0ic, p1ic)
@time psic    = fdtd(      xr, Nt,  f, S, p0ic, p1ic)
println("                       in-place call")
@time        fdtd!(psic2,  xr, Nt,  f, S, p0ic, p1ic)

@test vecnorm(psic-psic2) < 1e-8

println("---------- full no IC, source everywhere")
Ps    = fdtd(        Nt,geo, S)
Ps    = fdtd(        Nt,  f, S)
Ps2 = copy(Ps)
fdtd!(Ps2,    Nt,  f, S)
@time Ps    = fdtd(        Nt,  f, S)
println("                       in-place call")
@time        fdtd!(Ps2,    Nt,  f, S)

@test vecnorm(Ps-Ps2) < 1e-8

println("---------- full IC, source everywhere")
Psic    = fdtd(        Nt,geo, S, p0ic, p1ic)
Psic    = fdtd(        Nt,  f, S, p0ic, p1ic)
Psic2 = copy(Ps)
fdtd!(Psic2,    Nt,  f, S)
@time Psic    = fdtd(        Nt,  f, S, p0ic, p1ic)
println("                       in-place call")
@time        fdtd!(Psic2,    Nt,  f, S, p0ic, p1ic)

@test vecnorm(Psic-Psic2) < 1e-8

println("\ntest equivalence with matrix formulation. \n")

B = get_B(geo,Nt)
E = get_Exp(geo,Nt,xs)
F = get_Sel(geo,Nt,xr)
s = E*s_t[:]
s =[p0ic[:];p1ic[:];s[2*geo.Nx*geo.Ny*geo.Nz+1:end]]
S =[p0ic[:];p1ic[:];S[:]]
pm = F*(B\s)
println("sparse matrix      call")
@time pm = F*(B\s)
println("sparse matrix full call")
Pm = B\s
@time Pm = B\s

Pms = B\S

@test norm(PIC[:]-Pm[2*geo.Nx*geo.Ny*geo.Nz+1:end]) <1e-8
@test norm(pic[:]-pm[:]) <1e-8

@test norm(Psic[:]-Pms[2*geo.Nx*geo.Ny*geo.Nz+1:end]) <1e-8


println("\ntest resonances on cuboid room are correct \n")
Nt = 3000
X = 0.4
geo = CuboidRoom(10, 11, 12, 1e8*ones(6), FDTDEnv(X,IISO()))

using DSP
s_t = zeros(Nt)   #source signal
s_t[3] = 1
f2 = geo.env.Fs/2*0.14#cut-off frequency of source
filt!(s_t,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(15)),s_t)

xr = (2, 2, 2)
xs = (geo.Nx-1, geo.Ny-1, geo.Nz-1)
p_t2 = fdtd(xr, xs, Nt, geo, s_t)

f = linspace(0,geo.env.Fs/2,div(Nt,2)+1)
fres  = zeros(5,5,5)
for l = 0:4,m=0:4,i =0:4

	fres[l+1,m+1,i+1] = 343/(2*π)*sqrt(  (l*π/(geo.Lx))^2+(m*π/(geo.Ly))^2+(i*π/(geo.Lz))^2 )
end
fres = sort(fres[:])[2:7]

P = 10*log10(abs(rfft(p_t2)))

##using PyPlot
##figure()
##plot(f,P)
##plot(fres[:],5.* ones(length(fres[:])), "r*" )
##xlim([0,400])

@test norm(f[P.>5.5][1:6]-fres)<0.3 #test on resonance frequencies


println("\ntest acoustic delay is correct \n")
X = 0.4
Nt = 500
geo = CuboidRoom(20, 21, 24, 4*ones(6), FDTDEnv(X,IWB()))
s_t = zeros(Nt)   #source signal
s_t[3] = 1
f2 = geo.env.Fs/2*0.1#cut-off frequency of source
filt!(s_t,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(2)),s_t)

xs = (div(geo.Nx,2), div(geo.Ny,2)  , div(geo.Nz,2))
xr = (div(geo.Nx,2), div(geo.Ny,2)+4, div(geo.Nz,2))
d  = norm(X*[ xs[1]-xr[1],xs[2]-xr[2],xs[3]-xr[3]])

p_t2 = fdtd(xr, xs, Nt, geo, s_t)
t = 1/geo.env.Fs:1/geo.env.Fs:Nt/geo.env.Fs

##using PyPlot
##plot(t,p_t2)
##plot(d/343.+2*X/343, 0., "r*")#+2 due to initial conditions

#test delay
@test norm((indmax(xcorr(p_t2[:],s_t)[end-length(p_t2):end])-1)/geo.env.Fs-d/343.)<1e-3










