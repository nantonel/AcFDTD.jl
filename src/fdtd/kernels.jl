# with xr and xs
function fdtdKernel!{T}(n::Int,
		     p_out::Array{T}, 
		     p0::Array{T,3}, 
		     p1::Array{T,3}, 
		     p2::Array{T,3}, 
		     xr::Vector{Tuple{Int,Int,Int}}, 
		     xs::Vector{Tuple{Int,Int,Int}}, 
		     s::Array{T}, 
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
		
# full and xs
function fdtdKernel!{T}(n::Int,
		     p_out::Array{T,4}, #full output 
		     p0::Array{T,3}, 
		     p1::Array{T,3}, 
		     p2::Array{T,3}, 
		     xs::Vector{Tuple{Int,Int,Int}}, 
		     s::Array{T}, 
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

# xr and full s
function fdtdKernel!{T}(n::Int,
		     p_out::Array{T}, 
		     p0::Array{T,3}, 
		     p1::Array{T,3}, 
		     p2::Array{T,3}, 
		     xr::Vector{Tuple{Int,Int,Int}}, 
		     s::Array{T,4}, 
		     d1::T, 
		     d2::T, 
		     d3::T, 
		     d4::T, 
		     q::Vector{T}, 
		     G::Array{UInt8,3} )
	for i = 1:size(p0,3), m = 1:size(p0,2), l = 1:size(p0,1)
		if G[l,m,i] == 0x01
			air!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,s)
			#walls
		elseif G[l,m,i] == 0x02
			l1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x03
			le!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x04
			m1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x05
			me!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x06
			i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x07
			ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
			
			#edges
		elseif G[l,m,i] == 0x08
			l1m1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x09
			l1me!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0a
			l1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0b
			l1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)

		elseif G[l,m,i] == 0x0c
			lem1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0d
			leme!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0e
			lei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0f
			leie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)

		elseif G[l,m,i] == 0x10
			m1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x11
			m1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x12
			mei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x13
			meie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
			
			#corners
		elseif G[l,m,i] == 0x14
			l1m1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x15
			l1m1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x16
			l1mei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x17
			l1meie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)

		elseif G[l,m,i] == 0x18
			lem1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x19
			lem1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x1a
			lemei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x1b
			lemeie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		end
		## get pressure
		for k in eachindex(xr)
			if xr[k][1] == l && xr[k][2] == m && xr[k][3] == i 
				p_out[n,k] = p2[l,m,i]   
			end
		end
	end
end

# full p and full s
function fdtdKernel!{T}(n::Int,
		     p_out::Array{T,4}, #full output 
		     p0::Array{T,3}, 
		     p1::Array{T,3}, 
		     p2::Array{T,3}, 
		     s::Array{T,4}, 
		     d1::T, 
		     d2::T, 
		     d3::T, 
		     d4::T, 
		     q::Vector{T}, 
		     G::Array{UInt8,3} )
	for i = 1:size(p0,3), m = 1:size(p0,2), l = 1:size(p0,1) 
		if G[l,m,i] == 0x01
			air!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,s)
			#walls
		elseif G[l,m,i] == 0x02
			l1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x03
			le!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x04
			m1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x05
			me!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x06
			i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x07
			ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
			
			#edges
		elseif G[l,m,i] == 0x08
			l1m1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x09
			l1me!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0a
			l1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0b
			l1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)

		elseif G[l,m,i] == 0x0c
			lem1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0d
			leme!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0e
			lei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x0f
			leie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)

		elseif G[l,m,i] == 0x10
			m1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x11
			m1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x12
			mei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x13
			meie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
			
			#corners
		elseif G[l,m,i] == 0x14
			l1m1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x15
			l1m1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x16
			l1mei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x17
			l1meie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)

		elseif G[l,m,i] == 0x18
			lem1i1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x19
			lem1ie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x1a
			lemei1!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		elseif G[l,m,i] == 0x1b
			lemeie!(n,l,m,i,d1,d2,d3,d4,p0,p1,p2,q,s)
		end
		p_out[l,m,i,n] = p2[l,m,i]   
	end
end
