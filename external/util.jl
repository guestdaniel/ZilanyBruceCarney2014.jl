using DSP

function ffGn(N::Cint)::Vector{Cdouble}
    return zeros((N, ))
end

function decimate(signal::Vector{Cdouble}, q::Cint)
    return signal[1:20000]
    #return resample(signal, 1/q)
end
