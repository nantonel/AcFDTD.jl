export fdtdAdj, fdtdAdj!

#xr and xs, no IC
function fdtdAdj{T,N}(xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		      f::FDTD{T,N},
		      p::Array{T})
	s = fdtd(xs,xr,Nt,f,flipdim(p,1))
	return flipdim(s,1)
end

function fdtdAdj!{T}(s_out::Array{T},
		     xr::Vector{Tuple{Int,Int,Int}}, xs::Vector{Tuple{Int,Int,Int}},Nt::Int,
		     f::FDTD{T,3}, p::Array{T}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = Nt:-1:1
		fdtdKernel!(n,s_out,p0,p1,p2,xs,xr,p,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return s_out
end

#full and xs, no IC
function fdtdAdj{T,N}(xs::Vector{Tuple{Int,Int,Int}}, Nt::Int,
		      f::FDTD{T,N},
		      p::Array{T,4})
	s = fdtd(xs,Nt,f,flipdim(p,4))
	return flipdim(s,1)
end

function fdtdAdj!{T}(s_out::Array{T},
		     xs::Vector{Tuple{Int,Int,Int}},Nt::Int,
		     f::FDTD{T,3}, p::Array{T,4}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = Nt:-1:1
		fdtdKernel!(n,s_out,p0,p1,p2,xs,p,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return s_out
end

#full and source everywhere, no IC
function fdtdAdj{T,N}(Nt::Int,
		      f::FDTD{T,N},
		      p::Array{T,4})
	s = fdtd(Nt,f,flipdim(p,4))
	return flipdim(s,4)
end

function fdtdAdj!{T}(s_out::Array{T,4},
		     Nt::Int,
		     f::FDTD{T,3}, p::Array{T,4}) 
	resetIC!(f)
	p2 = f.p2
	p1 = f.p1
	p0 = f.p0
	for n = Nt:-1:1
		fdtdKernel!(n,s_out,p0,p1,p2,p,f.d[1],f.d[2],f.d[3],f.d[4],f.q,f.G)
		p2, p1, p0 = p0, p2, p1
	end
	return s_out
end

