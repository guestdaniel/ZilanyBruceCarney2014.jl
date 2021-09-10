using Test

const libihc = "/home/daniel/AuditoryNerveFiber.jl/external/libihc.so"

# First, we will test the super simplified C in test.c 
@testset "test.c" begin
    # Start by trying to call test_change_input_vector and checking it works
    @test begin
        x = zeros(10)
        ccall((:test_change_input_vector, libihc), Cvoid, (Ptr{Cdouble}, ), x)
        x[1] == 99.99
    end
    # Next, we write a function we pass to test_manipulate_cvector_with_julia_function
    @test begin
        # Create output vector
        output = zeros(10)
        # Write function in Julia 
        function testfunc(input::Ptr{Cdouble})
            temp = unsafe_wrap(Array, input, 10)
            temp[1] = 99.99
        end
        # Call 
        temp = @cfunction(x -> unsafe_wrap(Array, x, 10)[1] = 99.99, Vector{Cdouble}, (Vector{Cdouble}, ))
        ccall((:test_manipulate_cvector_with_julia_function, libihc), Cvoid, (Ptr{Cdouble}, Ptr{Cvoid}), output, temp)
        print(output)
        output[1] == 99.99
    end
end

