export FDTD, fdtd, fdtd!

immutable FDTD{T,N}
	p0::Array{T,N}
	p1::Array{T,N}
	p2::Array{T,N}
	d::Tuple{T,T,T,T}
	q::Vector{T}
	G::Array{UInt8,3}     
	env::FDTDEnv

	function FDTD(p0,p1,p2,geo)
		a,b,c,λ = geo.env.scheme.a, geo.env.scheme.b, geo.env.scheme.c, geo.env.scheme.λ
		d = (λ^2*(1-4*a+4*b), 
      		     λ^2*(a-2*b),
                     λ^2*b,
                     2*(1+λ^2*(-3+6*a-4*b)) )
		q = λ./(geo.ξ*2.0)
		new(p0,p1,p2,d,q,geo.G,geo.env)
	end


	FDTD(p0,p1,p2,d,q,G,env) = new(p0,p1,p2,d,q,G,env) 
	
end

"""
`FDTD(S::Type,geo::AbstractGeometry)`

Returns a `FDTD` object that can be used in `fdtd` and `fdtd!`. 
`S` must be an `AbstractFloat` and determinates the numerical precision of the FDTD simulation.
"""
function FDTD(S::Type,geo::AbstractGeometry)

	p0 = zeros(S,geo.Nx,geo.Ny,geo.Nz) 
	p1 = zeros(S,geo.Nx,geo.Ny,geo.Nz) 
	p2 = zeros(S,geo.Nx,geo.Ny,geo.Nz) 
	FDTD{S,3}(p0,p1,p2,geo)
end

FDTD(geo::AbstractGeometry) = FDTD(Float64, geo)

function setIC!{T,N}(f::FDTD{T,N},p0::Array{T,N},p1::Array{T,N})
	f.p0 .= p0
	f.p1 .= p1
end

function resetIC!{T,N}(f::FDTD{T,N})
	f.p0 .= 0.0
	f.p1 .= 0.0
end

fdtd(xr::Tuple{Int,Int,Int}, xs::Tuple{Int,Int,Int}, args...) = 
fdtd([xr], [xs], args...)

fdtd(xr::Vector{Tuple{Int,Int,Int}}, xs::Tuple{Int,Int,Int}, args...) = 
fdtd(xr, [xs], args...)

fdtd(xs::Tuple{Int,Int,Int}, args...) = 
fdtd([xs], args...)


"""
* `p = fdtd([xr, xs], Nt, geo, s, [p0ic, p1ic])`
  * `xs`  : microphone positions 
  * `xr`  : source positions (must be a `Tuple{Int,Int,Int}` or and `Vector{Tuple}` like `xs`)
  * `Nt`  : number of time samples `Int` 
  * `geo` : geometry specification, must be a `AbstractGeometry` or a `FDTD` object (the latter saves some memory allocation)
  * `s`   : 
    * if `xs` is given must be an `Array` of `size(s) = (Nt, length(xs))` where `s[:,k]` is the `k`th source signal positioned at `xs[k]`
    * otherwise must be a an `Array` of size `(size(geo)...,Nt)`
  * `p0ic` and `p1ic`: `Array` of `size(geo)` consisting of the initial conditions  		
	
  * reurns `p`: 
    * if `xr` is given `p` is an `Array` of size `(Nt, length(xr))` where the `p[:,k]` is the `k`th mic signal positioned at `xr[k]`
    * otherwise `p` will be an `Array` os size `(size(geo)...,Nt)`

"""
function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(xr,xs,Nt,f,s)
end

function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}}, Nt::Int, 
		 f::FDTD{T,3}, s::Array{T}) 
	size(s,2) != length(xs) && throw(ArgumentError("size(s,2) must be equal to length(xs)")) 

	#initialize output
	p_out = Array{T}(Nt,length(xr))
	fdtd!(p_out, xr, xs, Nt, f, s)
	return p_out
end

"""
`p = fdtd!(p::Array, xr, xs, Nt, f::FDTD, s, [p0ic, p1ic])`


In-place version of `fdtd` overwrites the `Array` `p`. See `fdtd` for more details.	
"""
function fdtd!{T}(p_out::Array{T},
		  xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}},Nt::Int,
		  f::FDTD{T,3}, s::Array{T}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,xr,xs,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end

## full with xr, null IC 
function fdtd{T}(xs::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(xs,Nt,f,s)
end

function fdtd{T}(xs::Vector{Tuple{Int,Int,Int}}, Nt::Int, 
		 f::FDTD{T,3}, s::Array{T}) 
	size(s,2) != length(xs) && throw(ArgumentError("size(s,2) must be equal to length(xs)")) 

	#initialize output
	p_out = Array{T}(size(f.p0)...,Nt)
	fdtd!(p_out, xs, Nt, f, s)
	return p_out
end

function fdtd!{T}(p_out::Array{T,4},
		  xs::Vector{Tuple{Int,Int,Int}},Nt::Int,
		  f::FDTD{T,3}, s::Array{T}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,xs,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end

## xr with source everywhere, null IC 
function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}},Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T,4}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(xr,Nt,f,s)
end

function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}},Nt::Int, 
		 f::FDTD{T,3}, s::Array{T,4}) 

	#initialize output
	p_out = Array{T}(Nt,length(xr))
	fdtd!(p_out, xr, Nt, f, s)
	return p_out
end

function fdtd!{T}(p_out::Array{T},
		  xr::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		  f::FDTD{T,3}, s::Array{T,4}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,xr,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end


## full with source everywhere, null IC 
function fdtd{T}(Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T,4}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(Nt,f,s)
end

function fdtd{T}(Nt::Int, 
		 f::FDTD{T,3}, s::Array{T,4}) 

	#initialize output
	p_out = Array{T}(size(f.p0)...,Nt)
	fdtd!(p_out, Nt, f, s)
	return p_out
end

function fdtd!{T}(p_out::Array{T,4},
		  Nt::Int,
		  f::FDTD{T,3}, s::Array{T,4}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end


## xr and xs with IC 

function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T},
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(xr,xs,Nt,f,s,p0ic,p1ic)
end

function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}}, Nt::Int, 
		 f::FDTD{T,3}, s::Array{T}, 
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 
	size(s,2) != length(xs) && throw(ArgumentError("size(s,2) must be equal to length(xs)")) 

	#initialize output
	p_out = Array{T}(Nt,length(xr))
	fdtd!(p_out, xr, xs, Nt, f, s, p0ic, p1ic)
	return p_out
end

function fdtd!{T}(p_out::Array{T},
		  xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}},Nt::Int,
		  f::FDTD{T,3}, s::Array{T}, 
		  p0ic::Array{T,3}, p1ic::Array{T,3}) 
	setIC!(f,p0ic,p1ic)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,xr,xs,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end


#full soundfield & IC
function fdtd{T}(xs::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T},
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(xs,Nt,f,s,p0ic,p1ic)
end

function fdtd{T}(xs::Vector{Tuple{Int,Int,Int}}, Nt::Int, 
		 f::FDTD{T,3}, s::Array{T}, 
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 
	size(s,2) != length(xs) && throw(ArgumentError("size(s,2) must be equal to length(xs)")) 

	#initialize output
	p_out = Array{T}(size(f.p0)...,Nt)
	fdtd!(p_out, xs, Nt, f, s, p0ic, p1ic)
	return p_out
end

function fdtd!{T}(p_out::Array{T,4},
		  xs::Vector{Tuple{Int,Int,Int}},Nt::Int,
		  f::FDTD{T,3}, s::Array{T}, 
		  p0ic::Array{T,3}, p1ic::Array{T,3}) 
	setIC!(f,p0ic,p1ic)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,xs,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end

## xr with source everywhere, & IC 
function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}},Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T,4}, 
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(xr,Nt,f,s,p0ic,p1ic)
end

function fdtd{T}(xr::Vector{Tuple{Int,Int,Int}},Nt::Int, 
		 f::FDTD{T,3}, s::Array{T,4},
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 

	#initialize output
	p_out = Array{T}(Nt,length(xr))
	fdtd!(p_out, xr, Nt, f, s,p0ic,p1ic)
	return p_out
end

function fdtd!{T}(p_out::Array{T},
		  xr::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		  f::FDTD{T,3}, s::Array{T,4}, 
		  p0ic::Array{T,3}, p1ic::Array{T,3}) 
	setIC!(f,p0ic,p1ic)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,xr,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end


## full with source everywhere, & IC 
function fdtd{T}(Nt::Int,
		 geo::AbstractGeometry,
		 s::Array{T,4},
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 
	f = FDTD(eltype(s),geo) 
	return fdtd(Nt,f,s,p0ic,p1ic)
end

function fdtd{T}(Nt::Int, 
		 f::FDTD{T,3}, s::Array{T,4},
		 p0ic::Array{T,3}, p1ic::Array{T,3}) 

	#initialize output
	p_out = Array{T}(size(f.p0)...,Nt)
	fdtd!(p_out, Nt, f, s,p0ic,p1ic)
	return p_out
end

function fdtd!{T}(p_out::Array{T,4},
		  Nt::Int,
		  f::FDTD{T,3}, s::Array{T,4}, 
		  p0ic::Array{T,3}, p1ic::Array{T,3}) 
	setIC!(f,p0ic,p1ic)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:Nt
		fdtdKernel!(n,p_out,p0,p1,p2,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return p_out
end




##walls
include("kernels.jl")

function air!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N})
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l+1,m,i-1] 
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l,m-1,i-1]
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l+1,m-1,i-1]
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m-1,i+1]
		sum3 += p1[l-1,m+1,i-1]
		sum3 += p1[l-1,m-1,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4-p0[l,m,i]
end

##adding source

function air!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},s::Array{T,4})
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l+1,m,i-1] 
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l,m-1,i-1]
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l+1,m-1,i-1]
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m-1,i+1]
		sum3 += p1[l-1,m+1,i-1]
		sum3 += p1[l-1,m-1,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4-p0[l,m,i]+s[l,m,i,n]
end

##walls
include("walls.jl")
##edges
include("edges.jl")
##corners
include("corners.jl")

