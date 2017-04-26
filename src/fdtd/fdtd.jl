export FDTD, fdtd, fdtd!, fdtdAdj

immutable FDTD{T<:AbstractFloat,N}
	p0::AbstractArray{T,N}
	p1::AbstractArray{T,N}
	p2::AbstractArray{T,N}
	p0ic::Union{AbstractArray{T,N},T}
	p1ic::Union{AbstractArray{T,N},T}
	xr::Vector{Tuple{Int,Int,Int}}
	xs::Vector{Tuple{Int,Int,Int}}
	d::Tuple{T,T,T,T}
	q::AbstractArray{T}
	Nt::Int
	G::Array{UInt8,3}     
	env::FDTDEnv

	function FDTD(p0,p1,p2,p0ic,p1ic,xr,xs,Nt,geo)
		a,b,c,λ = geo.env.scheme.a, geo.env.scheme.b, geo.env.scheme.c, geo.env.scheme.λ
		d = (λ^2*(1-4*a+4*b), 
      		     λ^2*(a-2*b),
                     λ^2*b,
                     2*(1+λ^2*(-3+6*a-4*b)) )
		q = λ./(geo.ξ*2.0)
		new(p0,p1,p2,p0ic,p1ic,xr,xs,d,q,Nt,geo.G,geo.env)
	end


	FDTD(p0,p1,p2,xr,xs,d,q,Nt,G,env) = new(p0,p1,p2,xr,xs,d,q,Nt,G,env) 
	
end
	
function FDTD(S::Type, 
	      xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}}, 
	      Nt::Int,geo::AbstractGeometry)

	p0 = zeros(S,geo.Nx,geo.Ny,geo.Nz) 
	p1 = zeros(S,geo.Nx,geo.Ny,geo.Nz) 
	p2 = zeros(S,geo.Nx,geo.Ny,geo.Nz) 
	FDTD{S,3}(p0,p1,p2,zero(S),zero(S),xr,xs,Nt,geo)
end

FDTD(S::Type, xs::Vector{Tuple{Int,Int,Int}}, Nt::Int, args...) = FDTD(S,[(0,0,0)],xs,Nt,args...)

function resetIC!(f::FDTD)
	f.p0 .= f.p0ic
	f.p1 .= f.p1ic
	f.p2 .= 0.
end


"""
`finite difference time domain method for room acoustic simulations`

# Usage 

* `p = fdtd(s::Array{Float64}, xr::Array{Int64}, xs::Array{Int64}, Nt::Int64, geo::CuboidRoom)`
  * `s` : must be a `Nt` × `size(xs,1)` matrix where the `k`th column contains the `k`th source signal
  * `xr`: source positions
  * `xs`: microphone positions (if this argument is not given sound pressure is return everywhere)
  * `Nt`: number of time samples
  * `geo`: geometry specification
	
* reurns `p`, a `Nt+2` × `size(xr,1)`  matrix where the `k`th column contains the `k`th mic signal

# Keyword arguments:
	
	
* `p0`: `Array{Float64,1}`, initial sound pressure at time `n-2`
* `p1`: `Array{Float64,1}`, initial sound pressure at time `n-1`
  * this must be `geo.Nx*geo.Ny*geo.Nz` long vectors

# Note

Everytime `fdtd` is run with the commands specified above some matrices some matrices are created.
If one has to call the function on the same geometry multiple times, one can store these matrices by running the command:

* `Qm,A,Qp = get_QmAQp(geo)`

and calling 'fdtd' with:

* `p = fdtd(s, xr, xs, Nt, geo, Qm, A, Qp)`.

"""
function fdtd{T<:AbstractFloat}(xr::Vector{Tuple{Int,Int,Int}},
				xs::Vector{Tuple{Int,Int,Int}},
				Nt::Int64,
				geo::AbstractGeometry,s::AbstractArray{T}) 

	f = FDTD(eltype(s),xr,xs,Nt,geo) 
	return fdtd(f,s)
end

fdtd(xr::Tuple{Int,Int,Int}, xs::Tuple{Int,Int,Int}, Nt::Int, args...) = fdtd([xr], [xs], Nt, args...)
fdtd(xs::Vector{Tuple{Int,Int,Int}}, Nt::Int, args...) = fdtd([(0,0,0)], xs, Nt, args...)
fdtd(xs::Tuple{Int,Int,Int}, Nt::Int, args...) = fdtd([(0,0,0)], [xs], Nt, args...)

function fdtd{T<:AbstractFloat}(f::FDTD{T,3},s::AbstractArray{T}) 
	size(s,2) != length(f.xs) && throw(ArgumentError("size(s,2) must be equal to length(xs)")) 
	#initialize output
	f.xr == [(0,0,0)] ? p_out = zeros(T,size(f.p0)...,f.Nt) : 
	p_out = zeros(T,f.Nt,length(f.xr))
	fdtd!(p_out,f,s)
	return p_out
end

function fdtd!{T<:AbstractFloat}(p_out::AbstractArray{T},f::FDTD{T,3},s::AbstractArray{T}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = 1:f.Nt
		fdtdKernel!(n,p_out,p0,p1,p2,f.xr,f.xs,s,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
end

#function fdtdAdj(p::Array{Float64},xs::Array{Int64}, xr::Array{Int64}, args...)
#	s = fdtd(flipdim(p,1),xr,xs,args...)
#	return flipdim(s,1)
#end

function fdtdKernel!{T<:AbstractFloat}(n::Int,
		     p_out::AbstractArray{T,2}, 
		     p0::AbstractArray{T,3}, 
		     p1::AbstractArray{T,3}, 
		     p2::AbstractArray{T,3}, 
		     xr::Vector{Tuple{Int,Int,Int}}, 
		     xs::Vector{Tuple{Int,Int,Int}}, 
		     s::AbstractArray{T,2}, 
		     d1::T, 
		     d2::T, 
		     d3::T, 
		     d4::T, 
		     q::Vector{T}, 
		     G::Array{UInt8,3} )
	for i = 1:size(p0,3), m = 1:size(p0,2), l = 1:size(p0,1) 
		if G[l,m,i] == 0x01
			air!(l,m,i,d1,d2,d3,d4,p0,p1,p2)
			
			#walls
		elseif G[l,m,i] == 0x02
			l1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x03
			le!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x04
			m1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x05
			me!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x06
			i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x07
			ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
			
			#edges
		elseif G[l,m,i] == 0x08
			l1m1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x09
			l1me!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0a
			l1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0b
			l1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)

		elseif G[l,m,i] == 0x0c
			lem1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0d
			leme!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0e
			lei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0f
			leie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)

		elseif G[l,m,i] == 0x10
			m1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x11
			m1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x12
			mei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x13
			meie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
			
			#corners
		elseif G[l,m,i] == 0x14
			l1m1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x15
			l1m1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x16
			l1mei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x17
			l1meie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)

		elseif G[l,m,i] == 0x18
			lem1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x19
			lem1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x1a
			lemei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x1b
			lemeie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		end
		## add source
		for k in eachindex(xs)
			if xs[k][1] == l && xs[k][2] == m && xs[k][3] == i 
				p2[l,m,i] += s[n,k] 
			end
		end
		## get pressure
		for k in eachindex(xr)
			if xr[k][1] == l && xr[k][2] == m && xr[k][3] == i 
				p_out[n,k] = p2[l,m,i]   
			end
		end
	end
end
		
function fdtdKernel!{T<:AbstractFloat}(n::Int,
		     p_out::AbstractArray{T,4}, #full output 
		     p0::AbstractArray{T,3}, 
		     p1::AbstractArray{T,3}, 
		     p2::AbstractArray{T,3}, 
		     xr::Vector{Tuple{Int,Int,Int}}, 
		     xs::Vector{Tuple{Int,Int,Int}}, 
		     s::AbstractArray{T,2}, 
		     d1::T, 
		     d2::T, 
		     d3::T, 
		     d4::T, 
		     q::Vector{T}, 
		     G::Array{UInt8,3} )
	for i = 1:size(p0,3), m = 1:size(p0,2), l = 1:size(p0,1) 
		if G[l,m,i] == 0x01
			air!(l,m,i,d1,d2,d3,d4,p0,p1,p2)
			
			#walls
		elseif G[l,m,i] == 0x02
			l1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x03
			le!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x04
			m1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x05
			me!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x06
			i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x07
			ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
			
			#edges
		elseif G[l,m,i] == 0x08
			l1m1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x09
			l1me!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0a
			l1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0b
			l1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)

		elseif G[l,m,i] == 0x0c
			lem1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0d
			leme!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0e
			lei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x0f
			leie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)

		elseif G[l,m,i] == 0x10
			m1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x11
			m1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x12
			mei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x13
			meie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
			
			#corners
		elseif G[l,m,i] == 0x14
			l1m1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x15
			l1m1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x16
			l1mei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x17
			l1meie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)

		elseif G[l,m,i] == 0x18
			lem1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x19
			lem1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x1a
			lemei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		elseif G[l,m,i] == 0x1b
			lemeie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
		end
		## add source
		for k in eachindex(xs)
			if xs[k][1] == l && xs[k][2] == m && xs[k][3] == i 
				p2[l,m,i] += s[n,k] 
			end
		end
		p_out[l,m,i,n] = p2[l,m,i]   
	end
end
	


function air!(l,m,i,d1,d2,d3,d4,p0,p1,p2)
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

function l1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
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
		#ghost
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l+1,m-1,i-1]
		#ghost
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l,m-1,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+1.0)
end

function l1m1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l+1,m,i-1] 
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m+1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1] 
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m+1,i-1]
		#ghost
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[3]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[3]+1.0)
end

function lem1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m+1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1] 
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m+1,i-1]
		#ghost
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l-1,m,i-1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[3]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[3]+1.0)
end

function l1me!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l+1,m,i-1] 
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1] 
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l+1,m-1,i-1]
		#ghost
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l,m-1,i-1]
		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[4]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[4]+1.0)
end

function leme!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1] 
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m-1,i+1]
		sum3 += p1[l-1,m-1,i-1]
		#ghost
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l,m-1,i-1]
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l-1,m,i-1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[4]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[4]+1.0)
end

function l1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m-1,i+1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i] 
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m-1,i+1]
		#ghost
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m-1,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[5]+1.0)
end

function lei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m-1,i+1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i] 
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m-1,i+1]
		#ghost
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l-1,m-1,i]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m-1,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[5]+1.0)
end

function leie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i-1]
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i-1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i] 
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m+1,i-1]
		sum3 += p1[l-1,m-1,i-1]
		#ghost
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l,m-1,i-1]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l-1,m-1,i]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m-1,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[6]+1.0)
end

function l1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i-1]
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i-1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i] 
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l+1,m-1,i-1]
		#ghost
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l,m-1,i-1]
		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m-1,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[6]+1.0)
end


function le!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m-1,i+1]
		sum3 += p1[l-1,m+1,i-1]
		sum3 += p1[l-1,m-1,i-1]
		#ghost
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l,m-1,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+1.0)
end

function m1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i] 
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 

		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l+1,m,i-1] 
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 

		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m+1,i-1]

		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l-1,m,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 = d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[3]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[3]+1.0)
end

function me!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m,i-1] 

		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l,m-1,i-1]
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l+1,m,i-1] 
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l-1,m,i-1] 

		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l+1,m-1,i-1]
		sum3 += p1[l-1,m-1,i+1]
		sum3 += p1[l-1,m-1,i-1]

		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l-1,m,i-1]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 = d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[4]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[4]+1.0)
end

function i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i+1]

		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]


		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l-1,m+1,i+1]
		sum3 += p1[l-1,m-1,i+1]

		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l-1,m-1,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[5]+1.0)
end


function ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m-1,i] 
		sum1 += p1[l,m,i-1]

		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m,i-1]
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l,m-1,i-1]
		sum2 += p1[l-1,m,i-1]
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m-1,i]


		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l+1,m-1,i-1]
		sum3 += p1[l-1,m+1,i-1]
		sum3 += p1[l-1,m-1,i-1]

		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l-1,m-1,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[6]+1.0)
end

function m1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l,m+1,i+1]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l-1,m,i+1] 
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l-1,m+1,i]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i] 
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		sum3 += p1[l-1,m+1,i+1]
		#ghost
		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l-1,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[3]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[3]+q[5]+1.0)
end

function m1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l,m,i-1]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l,m+1,i-1]
		sum2 += p1[l+1,m,i-1]
		sum2 += p1[l-1,m,i-1] 
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l-1,m+1,i]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i] 
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i-1]
		sum3 += p1[l-1,m+1,i-1]
		#ghost
		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l-1,m,i-1]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l-1,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[3]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[3]+q[6]+1.0)
end


function mei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l,m,i+1]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l,m-1,i+1]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l-1,m,i+1] 
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l-1,m-1,i]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i] 
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m-1,i+1]
		sum3 += p1[l-1,m-1,i+1]
		#ghost
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l-1,m-1,i]
		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l-1,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[4]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[4]+q[5]+1.0)
end

function meie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l,m,i-1]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l+1,m,i]
		sum1 += p1[l-1,m,i] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l,m-1,i-1]
		sum2 += p1[l+1,m,i-1]
		sum2 += p1[l-1,m,i-1] 
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l-1,m-1,i]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i] 
		sum2 += p1[l+1,m,i]
		sum2 += p1[l-1,m,i]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m-1,i-1]
		sum3 += p1[l-1,m-1,i-1]
		#ghost
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l-1,m-1,i]
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l-1,m,i-1]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l-1,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[4]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[4]+q[6]+1.0)
end

##corners

function l1m1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m,i+1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l,m+1,i+1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i+1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i+1]
		#ghost
		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[3]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[3]+q[5]+1.0)
end

function l1mei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l,m,i+1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i+1]
		sum2 += p1[l,m-1,i+1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i+1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m-1,i+1]
		#ghost
		sum3 += p1[l+1,m,i+1]
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m-1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[4]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[4]+q[5]+1.0)
end

function l1m1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m+1,i]
		sum2 += p1[l+1,m,i-1]
		sum2 += p1[l,m+1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m+1,i-1]
		#ghost
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l+1,m+1,i]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l,m,i-1]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[3]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[3]+q[6]+1.0)
end

function l1meie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l+1,m,i]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l+1,m-1,i]
		sum2 += p1[l+1,m,i-1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l+1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l+1,m-1,i-1]
		#ghost
		sum3 += p1[l+1,m,i-1]
		sum3 += p1[l+1,m-1,i]
		sum3 += p1[l+1,m,i]
		sum3 += p1[l,m-1,i-1]
		sum3 += p1[l,m,i-1]
		sum3 += p1[l,m-1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[4]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[1]+q[4]+q[6]+1.0)
end

function lem1i1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m,i+1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l,m+1,i+1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i+1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m+1,i+1]
		#ghost
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l-1,m,i]
		sum3 += p1[l,m+1,i+1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[3]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[3]+q[5]+1.0)
end

function lemei1!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l,m,i+1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i+1]
		sum2 += p1[l,m-1,i+1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i+1]
		sum2 += p1[l,m,i+1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m-1,i+1]
		#ghost
		sum3 += p1[l-1,m,i+1]
		sum3 += p1[l-1,m-1,i]
		sum3 += p1[l-1,m,i]
		sum3 += p1[l,m-1,i+1]
		sum3 += p1[l,m,i+1]
		sum3 += p1[l,m-1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[4]+q[5]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[4]+q[5]+1.0)
end

function lem1ie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m+1,i]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m+1,i]
		sum2 += p1[l-1,m,i-1]
		sum2 += p1[l,m+1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m+1,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m+1,i-1]
		#ghost
		sum3 += p1[l-1,m,i-1]
		sum3 += p1[l-1,m+1,i]
		sum3 += p1[l-1,m,i]
		sum3 += p1[l,m+1,i-1]
		sum3 += p1[l,m,i-1]
		sum3 += p1[l,m+1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[3]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[3]+q[6]+1.0)
end

function lemeie!(l,m,i,d1,d2,d3,d4,p0,p1,p2,q)
	sum1,sum2,sum3,sum4 = 0.,0.,0.,0.
	if d1 != 0.
		sum1 += p1[l-1,m,i]
		sum1 += p1[l,m-1,i]
		sum1 += p1[l,m,i-1] 
		#ghost
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 += p1[l,m,i]
		sum1 *= d1
	end
	if d2 != 0.
		sum2 += p1[l-1,m-1,i]
		sum2 += p1[l-1,m,i-1]
		sum2 += p1[l,m-1,i-1]
		#ghost
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l,m,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l-1,m,i]
		sum2 += p1[l,m-1,i]
		sum2 += p1[l,m,i-1]
		sum2 += p1[l,m,i-1]
		sum2 *= d2
	end
	if d3 != 0.
		sum3 += p1[l-1,m-1,i-1]
		#ghost
		sum3 += p1[l-1,m,i-1]
		sum3 += p1[l-1,m-1,i]
		sum3 += p1[l-1,m,i]
		sum3 += p1[l,m-1,i-1]
		sum3 += p1[l,m,i-1]
		sum3 += p1[l,m-1,i]
		sum3 += p1[l,m,i]
		sum3 *= d3
	end
	if d4 != 0.
		sum4 += d4*p1[l,m,i]
	end
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[4]+q[6]-1.)*p0[l,m,i]
	p2[l,m,i] /= (q[2]+q[4]+q[6]+1.0)
end
