export fdtd

"""
`finite difference time domain method for room acoustic simulations`

# Usage 

* `p = fdtd(s::Array{Float64}, xs::Array{Int64}, xr::Array{Int64}, Nt::Int64, geo::CuboidRoom)`
  * `s` : must be a `Nt` × `size(xs,1)` matrix where the `k`th column contains the `k`th source signal
  * `xs`: source positions
  * `xr`: microphone positions (if this argument is not given sound pressure is return everywhere)
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

* `p = fdtd(s, xs, xr, Nt, geo, Qm, A, Qp)`.

"""

function fdtd(s::Array{Float64}, xs::Array{Int64}, xr::Array{Int64}, Nt::Int64, geo::CuboidRoom;
	     kwargs...) 

	Qm,A,Qp = get_QmAQp(geo)

	return fdtd(s,xs,xr,Nt,geo,Qm,A,Qp;kwargs...)

end

function fdtd(s::Array{Float64}, xs::Array{Int64}, Nt::Int64, geo::CuboidRoom, args...; kwargs...) 
	xr = zeros(Int64,3,geo.Nx*geo.Ny*geo.Nz)

	for i = 1:geo.Nx*geo.Ny*geo.Nz
		xr[1,i],xr[2,i],xr[3,i] =  ind2sub((geo.Nx,geo.Ny,geo.Nz),i)
	end

	return fdtd(s,xs,xr,Nt,geo,args...;kwargs...)
end

function fdtd(s::Array{Float64}, xs::Array{Int64}, xr::Array{Int64}, Nt::Int64, geo::CuboidRoom, 
	      Qm::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, Qp::Array{Float64,1}; 
	      p1::Array{Float64,1} = zeros(Float64, geo.Nx*geo.Ny*geo.Nz ), #initial conditions
	      p0::Array{Float64,1} = zeros(Float64, geo.Nx*geo.Ny*geo.Nz ) ) 

	if(size(xs,1)!=3) error("size(xs,1) must be equal to 3") end
	if(size(xs,2)!= size(s,2)) error("size(xs,1) must be equal to size(s,2)") end
	if(size(xr,1)!=3) error("size(xr,1) must be equal to 3") end
	if(size(s,1) != Nt) s = [s; zeros(Nt-size(s,1),size(s,2))]  end

	Nxyz = geo.Nx*geo.Ny*geo.Nz
	#initialize arrays
	p2 = zeros(Float64, Nxyz )

	#initialize output
	p_out = zeros(Float64,Nt+2,size(xr,2))
		
	indxr = sub2ind((geo.Nx,geo.Ny,geo.Nz),xr[1,:],xr[2,:],xr[3,:])
	indxs = sub2ind((geo.Nx,geo.Ny,geo.Nz),xs[1,:],xs[2,:],xs[3,:])

	p_out[1,:] = copy(p0[indxr]) 
	p_out[2,:] = copy(p1[indxr]) 

	for n = 3:Nt+2

		p2 = A*p1 + Qm.*p0
		p2[indxs] += s[n-2,:]
		p2 = p2./Qp

		p_out[n,:] = copy(p2[indxr]) 

		copy!(p0,p1)
		copy!(p1,p2)

	end

	return p_out
end

