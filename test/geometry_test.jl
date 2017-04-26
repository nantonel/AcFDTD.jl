
Nx,Ny,Nz = 10,20,11
G = ones(Bool,Nx+2,Ny+2,Ny+2) 
G[1,:,:] = G[end,:,:] = false
G[:,1,:] = G[:,end,:] = false
G[:,:,1] = G[:,:,end] = false

I = AcFDTD.getNodeProperty(G)
@test all(I[2:end-1,2:end-1,2:end-1] .== 0x01)
#walls
@test all(I[1  ,2:end-1,2:end-1].==0x02)
@test all(I[end,2:end-1,2:end-1].==0x03)
@test all(I[2:end-1,  1,2:end-1].==0x04)
@test all(I[2:end-1,end,2:end-1].==0x05)
@test all(I[2:end-1,2:end-1,  1].==0x06)
@test all(I[2:end-1,2:end-1,end].==0x07)
#edges
@test all(I[  1,  1,2:end-1].==0x08)
@test all(I[  1,end,2:end-1].==0x09)
@test all(I[  1,2:end-1,  1].==0x0a)
@test all(I[  1,2:end-1,end].==0x0b)

@test all(I[end,  1,2:end-1].==0x0c)
@test all(I[end,end,2:end-1].==0x0d)
@test all(I[end,2:end-1,  1].==0x0e)
@test all(I[end,2:end-1,end].==0x0f)

@test all(I[2:end-1,   1,  1].==0x10)
@test all(I[2:end-1,   1,end].==0x11)
@test all(I[2:end-1, end,  1].==0x12)
@test all(I[2:end-1, end,end].==0x13)

#corners
@test I[  1,  1,  1] == 0x14
@test I[  1,  1,end] == 0x15
@test I[  1,end,  1] == 0x16
@test I[  1,end,end] == 0x17

@test I[end,  1,  1] == 0x18
@test I[end,  1,end] == 0x19
@test I[end,end,  1] == 0x1a
@test I[end,end,end] == 0x1b

