using AcFDTD
using Base.Test

srand(123)

include("schemes_test.jl")
include("geometry_test.jl")
include("fdtd_test.jl")
include("adjoint_test.jl")

