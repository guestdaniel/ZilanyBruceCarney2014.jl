using DSP

function ffGn(N::Int32)
    return zeros((N, ))
end

function decimate(original_signal::Ptr{Cdouble}, k::Int32, resamp::Int32)
    temp_orig = unsafe_wrap(Array, original_signal, k)
    _resampled = resample(temp_orig, 1/resamp)
    return pointer(_resampled)
end

function test(signal::Ptr{Cdouble})
    temp = unsafe_wrap(Array, signal, 100000)
    #print(temp)
    return temp
end

function test_float(x::Float64)
    print(x)
    print(' ')
    return nothing
end

function test_array(x::Ptr{Cdouble})
    temp = unsafe_wrap(Array, x, 100000)
    for ii = 1:100000
        if ii % 10000 == 0
            print(temp[ii])
            print(' ')
        end
    end
    print('\n')
    return nothing
end
