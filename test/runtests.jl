
import SpatialGSP as GSP

import Random
Random.seed!(25)

using LinearAlgebra
using Test

@testset "kernel params solve" begin

    N_tests = 500
    T = Float64
    
    zero_tol = eps(T)*100
    kernel_ref = GSP.SqExp(T)

    for _ = 1:N_tests
        
        x = convert(T, 50.0*rand(T))
        lb = rand(T)

        a = GSP.solveparam(kernel_ref, x, lb)
        kernel = GSP.createkernel(kernel_ref, a)
        @test abs( GSP.evalkernel(kernel, x) - lb) < zero_tol
    end

end