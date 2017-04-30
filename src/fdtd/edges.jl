function l1m1!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function lem1!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function l1me!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function leme!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function l1i1!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function lei1!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function leie!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function l1ie!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function m1i1!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function m1ie!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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


function mei1!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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

function meie!{T,N}(l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T})
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



















#adding source
function l1m1!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[3]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[1]+q[3]+1.0)
end

function lem1!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[3]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[2]+q[3]+1.0)
end

function l1me!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[4]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[1]+q[4]+1.0)
end

function leme!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[4]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[2]+q[4]+1.0)
end

function l1i1!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[5]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[1]+q[5]+1.0)
end

function lei1!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[5]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[2]+q[5]+1.0)
end

function leie!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[2]+q[6]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[2]+q[6]+1.0)
end

function l1ie!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[1]+q[6]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[1]+q[6]+1.0)
end

function m1i1!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[3]+q[5]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[3]+q[5]+1.0)
end

function m1ie!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[3]+q[6]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[3]+q[6]+1.0)
end


function mei1!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[4]+q[5]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[4]+q[5]+1.0)
end

function meie!{T,N}(n::Int,l::Int,m::Int,i::Int,
		   d1::T,d2::T,d3::T,d4::T,
		   p0::Array{T,N},p1::Array{T,N},p2::Array{T,N},q::Vector{T},s::Array{T,4})
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
	p2[l,m,i]  = sum1+sum2+sum3+sum4+(q[4]+q[6]-1.)*p0[l,m,i]+s[l,m,i,n]
	p2[l,m,i] /= (q[4]+q[6]+1.0)
end
