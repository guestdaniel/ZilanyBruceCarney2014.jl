using Test

const libihc = "/home/daniel/ZilanyBruceCarney2014.jl/external/libihc.so"

"""
    edit_one_element_of_c_array(a)

Attempts to take a pointer and edit the value of the first element
"""
function edit_one_element_of_c_array(a::Ptr{Cdouble})
    a_wrapped = unsafe_wrap(Array, a, 10)
    a_wrapped[1] = 99.99
    return nothing
end

# First, we will test the super simplified C in test.c 
@testset "test.c" begin
    # Start by trying to call test_change_input_vector()
    # test_change_input_vector() is a C function that takes a double * as input
    # The function then edits the first element of that array to be 99.99.
    # If we can read out that value correctly on the Julia side, we can be sure that we can edit Julia arrays in C.
    @test begin
        x = zeros(10)
        ccall((:test_change_input_vector, libihc), Cvoid, (Ptr{Cdouble}, ), x)
        x[1] == 99.99
    end
    # Next, we try to test test_manipulate_cvector_with_julia_function
    # This function accepts a void * as an input which should be a pointer to a function
    # In this case, we test the edit_one_element_of_c_array() function, which is a Julia function that accepts 
    # a double * as input and attempts to edit one element of it in-place. This merely confirms that editing the
    # value does not immediately throw an error or cause a seg fault.
    @test begin
        # Call 
        temp = @cfunction(edit_one_element_of_c_array, Cvoid, (Ptr{Cdouble}, ))
        ccall((:test_manipulate_cvector_with_julia_function, libihc), Cvoid, (Ptr{Cvoid}, ), temp)
        true
    end
    # Next we try to test test_manipulate_julia_vector_with_julia_function
    # This is the same as above, except instead of generating a vector and handling memory in C, we
    # do it in Julia and then pass it through to C.
    @test begin
        # Create output vector
        output = zeros(10)
        # Call 
        temp = @cfunction(edit_one_element_of_c_array, Cvoid, (Ptr{Cdouble}, ))
        ccall((:test_manipulate_julia_vector_with_julia_function, libihc), Cvoid, (Ptr{Cdouble}, Ptr{Cvoid}), output, temp)
        output[1] == 99.99
    end
end