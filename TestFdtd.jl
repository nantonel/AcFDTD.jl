using PyPlot
using DSP

include("src/Fdtd.jl")
include("src/FdtdGeo.jl")
include("src/FdtdLossMatrices.jl")

Nx, Ny, Nz = 10, 11, 12           # Size of grid      
λ, a, b = sqrt(3/4), 1/6, 0       # IISO scheme
#λ, a, b = sqrt(1/3), 0, 0        # SLF  scheme
X = 4/(Nx)                        # Spatial step
Fs = 343/(X*λ)                    # Sampling Frequency
Nt = iround(7*Fs)                 # Time samples (7 secs) 

s = zeros(Nt)
s[5] = 1       #impulse
fl,fh = 10,400
myHP = digitalfilter(Bandpass(fl,fh,fs = Fs),Butterworth(15))
s = filt(myHP,s) #bandpass impulse

ξ = [1000, 1000, 1000, 1000, 1000, 1000] # faces main impedance
υ = 1./ξ                                 # convert to admittances

pos = [2 2 2]'         # source position
posm = [Nx-1 2 Nz-1]'  # microphone position

G = CreateGeometry(Nx,Ny,Nz)
A, L, Gx,Gy,Gz = GetGeometryMatrix(G, a, b, λ^2)
Qp, Qm = GetLossMatrixFaces(Gx,Gy,Gz, υ, λ)

p_tild = Fdtd(Nx,Ny,Nz,Nt,A,Qp,Qm,s,pos,posm)

G2 = CreateGeometryLshaped(Nx,Ny,Nz)
A2, L2, Gx2,Gy2,Gz2 = GetGeometryMatrix(G2, a, b, λ^2)
Qp2, Qm2 = GetLossMatrixFaces(Gx2,Gy2,Gz2, υ, λ)

p_tild2 = Fdtd(Nx,Ny,Nz,Nt,A2,Qp2,Qm2,s,pos,posm)
Lx,Ly,Lz = X*Nx,X*Ny,X*Nz
fres  = zeros(5,5,5)

for l = 0:4,m=0:4,i =0:4

	fres[l+1,m+1,i+1] = 343/(2*π)*sqrt(  (l*π/Lx)^2+(m*π/Ly)^2+(i*π/Lz)^2 )
end

t = linspace(0,1/Fs*length(s),length(s))
ff = linspace(0,Fs, Nt)
figure()
plot(ff,10*log10(abs(fft(p_tild))),  label = "cuboid")
plot(ff,10*log10(abs(fft(p_tild2))), label = "L-shaped")
plot(fres[:],10.* ones(length(fres[:])), "r*",  label = "cuboid analytical resonances" )
xlim([20,170])
legend(loc = 3 )

figure()
plot(t, p_tild)
plot(t, p_tild2)
