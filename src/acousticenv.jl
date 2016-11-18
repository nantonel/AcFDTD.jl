export FDTDEnv

abstract Acoustics
abstract AcousticEnvironment <:Acoustics
abstract FDTDAcousticEnvironment <: AcousticEnvironment

immutable FDTDEnv <: FDTDAcousticEnvironment
	Fs::Float64        # sampling frequency
	c::Float64         # speed of sound
	X::Float64         # spatial step
	T::Float64         # temp. step
        scheme::FDTDScheme #FDTD scheme
end

"""
Create acoustic environment for FDTD simulation 

## Usage

* 'FDTD(X::Float64, scheme::FDTDScheme)' 
   * Create an acoustic env. for a given spatial step `X` 
* 'FDTD(scheme::FDTDScheme, Fs::Float64)': 
   * Create an acoustic env. for a given sampling frequency `Fs`
"""
function FDTDEnv(X::Float64, scheme::FDTDScheme; c = 343.) 

	Fs = 343/(X*scheme.λ)   # Sampling Frequency
	return FDTDEnv(Fs,c,X,1/Fs,scheme)

end

function FDTDEnv(scheme::FDTDScheme, Fs::Float64; c = 343.) 

	X = 343/(Fs*scheme.λ)   # Sampling Frequency
	return FDTDEnv(Fs,c,X,1/Fs,scheme)

end

function Base.show(io::IO, env::FDTDEnv)
	println(io, "FDTD Acoustic Env.")
	println(io, "Sampling freq. : $(round(env.Fs/1000,1))   kHz")
	println(io, "Spatial step   : $(round(env.X*100,2)) cm")
	println(io, "Speed of sound : $(round(env.c,1)) m/s")
	println(io, "FDTD Scheme    : $(env.scheme.name)")
end
