# AcFDTD (this library is unmantained) 

Finite Difference Time Domain (FDTD) method for room acoustic simulation

## Installation

From the Julia command line hit:

```julia
Pkg.clone("https://github.com/nantonel/AcFDTD.jl.git")
```

Once the package is installed you can update it along with the others issuing `Pkg.update()` in the command line.

## Usage 

Import the package by typing `using AcFDTD`. 
First you need to specify an acoustic environment 
and FDTD scheme: 
```julia
using AcFDTD
X = 0.1                   # spatial sampling
env = FDTDEnv(X,IISO())   # create new acoustic env with default values
```
where `IISO()` returns the Interpolated Isotropic scheme.
Alternatively one can choose a samplng frequency 
instead of a spatial sampling:
```julia
Fs = 2000.                          # sampling frequency in Hz
env = FDTDEnv(IISO(),Fs; c = 340)   # create new acoustic env with default values
```
notice that in the latter line the speed of sound was 
chosen to be `340` m/s. 
By default this is set to `343` m/s.
Set the acoustic impedance `ξ` and room geometry room geometry:
```julia
ξ = [50.;50.;100.;30.;50.;50.]; # [   ξx1    ;    ξx2   ;    ξy1   ;    ξy2    ;  ξz1 ;   ξz2  ]
                                # [front wall; rear wall; left wall; right wall; floor; ceiling]
geo = CuboidRoom(10, 11, 12, ξ, env)
```
The first three parameters 
indicate the number of spatial samples 
of the `x`, `y` and `z` directions.  
Alternatively one can specify 
the room dimensions in meters:
```julia
geo = CuboidRoom(4., 5., 3., ξ, env)
```
which are then approximated on the grid.
Set the number of time steps `Nt`
Create a band-limited sound source
with e.g. the `DSP` package:
```julia
using DSP
Nt = round(Int,env.Fs)         # number of time steps (1 sec)
s = zeros(Nt)                    # source signal
s[3] = 1
f2 = geo.env.Fs/2*0.175          # cut-off frequency of source
filt!(s,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(5)),s)
```
Define the position of microphone 
and sound sources:
```julia
xr = [(2, 2, 2), (2, geo.Ny-1, geo.Nz-1)] # mic positions
xs = (geo.Nx-1, geo.Ny-1, geo.Nz-1) # sound source position
#positions must be Tuple{Int,Int,Int} or Array{Tuple{Int,Int,Int},1} 
```
Now type:
```julia
p = fdtd(xr,xs,Nt,geo,s)
```
to obtain the sound pressure of the microphones.

For more details on the methods type:
```julia
?fdtd
```

## Credits

AcFDTD.jl is developed by [Niccolò Antonello](http://homes.esat.kuleuven.be/~nantonel/) at [KU Leuven, ESAT/Stadius](https://www.esat.kuleuven.be/stadius/).











