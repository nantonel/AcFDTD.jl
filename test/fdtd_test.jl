using AcFDTD
using Base.Test

geo = CuboidRoom(7,8,10,10*ones(6), FDTDEnv(0.1,SLF()))

srand(123)
Nt = 500
#@time L = get_L(geo)
#@time Q = get_Q(geo)
B = get_B(geo,Nt)
@time B = get_B(geo,Nt)

xs = [6 7 9; 3 3 3]'
E = get_Exp(geo,Nt,xs)
s_t =[exp(-((1:Nt)-30).^2./20) 1e-2*randn(Nt)] 

s = E*s_t[:]
p = B\s
@time p = B\s

xr = [2 3 2; 4 4 4]'
F = get_Sel(geo,Nt,xr)
p_t = F*p

p_t2 = fdtd(s_t,xs,xr,Nt,geo)
@time p_t2 = fdtd(s_t,xs,xr,Nt,geo)

Qm,A,Qp = get_QmAQp(geo)

p_t2 = fdtd(s_t,xs,xr,Nt,geo,Qm,A,Qp)
@time p_t2 = fdtd(s_t,xs,xr,Nt,geo,Qm,A,Qp)
#test matrix inv and fdtd return same stuff
@test norm(p_t2[:]-p_t)<=1e-11

@time p2 = fdtd(s_t,xs,Nt,geo,Qm,A,Qp)
@test norm(p2'[:]-p[2*geo.Nx*geo.Ny*geo.Nz+1:end])<=1e-10

p0,p1 = randn(geo.Nx*geo.Ny*geo.Nz),randn(geo.Nx*geo.Ny*geo.Nz)
s = zeros(s)
s[1:2*geo.Nx*geo.Ny*geo.Nz]=[p0;p1]
p = B\s
p_t = F*p

p_t2 = fdtd(zeros(s_t),xs,xr,Nt,geo,Qm,A,Qp; pic0 = p0,pic1 = p1)
@test norm(p_t2[:]-p_t)<=1e-11


Nt = 3000
X = 0.4
geo = CuboidRoom(10, 11, 12, 1e8*ones(6), FDTDEnv(X,IISO()))
using DSP
s_t = zeros(Nt)   #source signal
s_t[3] = 1
f2 = geo.env.Fs/2*0.14#cut-off frequency of source
filt!(s_t,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(15)),s_t)

xr = [2 2 2]'
xs = [geo.Nx-1 geo.Ny-1 geo.Nz-1]'
p_t2 = fdtd(s_t,xs,xr,Nt,geo)

f = linspace(0,geo.env.Fs/2,round(Int64,Nt/2)+1)
fres  = zeros(5,5,5)
for l = 0:4,m=0:4,i =0:4

	fres[l+1,m+1,i+1] = 343/(2*π)*sqrt(  (l*π/(geo.Lx))^2+(m*π/(geo.Ly))^2+(i*π/(geo.Lz))^2 )
end
fres = sort(fres[:])[2:7]

P = 10*log10(abs(rfft(p_t2)))

#using PyPlot
#figure()
#plot(f,P)
#plot(fres[:],5.* ones(length(fres[:])), "r*" )
#xlim([0,400])

@test norm(f[P.>5.5][1:6]-fres)<0.3 #test on resonance frequencies

X = 0.4
Nt = 500
geo = CuboidRoom(20, 21, 24, 4*ones(6), FDTDEnv(X,IWB()))
s_t = zeros(Nt)   #source signal
s_t[3] = 1
f2 = geo.env.Fs/2*0.1#cut-off frequency of source
filt!(s_t,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(2)),s_t)

xs = round(Int64,[geo.Nx/2 geo.Ny/2 geo.Nz/2]')
xr = round(Int64,[geo.Nx/2 geo.Ny/2+4 geo.Nz/2]')

p_t2 = fdtd(s_t,xs,xr,Nt,geo)
t = 1/geo.env.Fs:1/geo.env.Fs:Nt/geo.env.Fs

d = norm(X*(xr-xs))

#using PyPlot
#plot(t,p_t2)
#plot(d/343.+2*X/343, 0., "r*")#+2 due to initial conditions

#test delay
@test norm((indmax(xcorr(p_t2[:],s_t)[end-length(p_t2):end])-1)/geo.env.Fs
	   -norm(X*(xr-xs))/343.)<1e-3










