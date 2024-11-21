export random_numbers, decimate, upsample, ffGn_native

"""
    random_numbers(n)

Generates random n numbers like MATLAB's rand or Python's numpy.random.rand.

Note that this function returns a pointer instead of a Julia array as we expect that it will
be called from C and not Julia.
"""
function random_numbers(n::Int32)
    return pointer(rand(Float64, n))
end

"""
    ffGn_old(N)

Synthesize a sample of fractional Gaussian noise.

Adapted from the original MATLAB code avaialble at
https://www.urmc.rochester.edu/labs/carney/publications-code/auditory-models.aspx.

Note that this function returns a pointer instead of a Julia array, as it is expected to be
called from within C.

# Arguments
- `N::Int32`: length of the output sequence
- `tdres::Float64`: reciprocal of the sampling rate (1/Hz)
- `Hinput::Float64`: Hurst index. For 0 < H <= 1, we synthesize fractional Gaussian noise with Hurst index H. For 1 < H <= 2, we synthesize fractional Brownian motion with Hurst index H-1.
- `noiseType::Float64`: If noiseType == 0, we return zeros, else we return fractional Gaussian noise
- `mu::Float64`: Mean of the noise

# Warnings
- This code was adapted from the MATLAB original, but has not been tested against the original. Use caution and test carefully! Report any bugs on GitHub.
"""
function ffGn_old(
    N::Int32,
    tdres::Float64,
    Hinput::Float64,
    noiseType::Float64,
    mu::Float64;
    safety::Int64=4
)
    # Start by handling noiseType
    if noiseType == 0.0
        y = zeros(N)
        return pointer(y)
    end

    # If noiseType != 0, we're synthesizing fractional Gaussian noise
    # First, we downsample the number of points
    resamp = Int(ceil(1e-1 / tdres))
    nop = N
    N = Int(ceil(N / resamp) + 1)
    if N < 10
        N = 10
    end

    # Next, determine if fGn or fBm should be produced
    if Hinput < 1.0
        H = Hinput
        fBn = false
    else
        H = Hinput - 1
        fBn = true
    end

    # Calculate fGn
    if H == 0.5
        y = randn(N)
    else
        Nfft = Int(2 ^ ceil(log2(2*(N-1))))
        NfftHalf = Int(round(Nfft / 2))

        k = [0:(NfftHalf-1); NfftHalf:-1:1]
        Zmag = 0.5 * ( (k.+1) .^ (2*H) -2*k .^ (2*H) + abs.(k .- 1) .^ (2*H))

        Zmag = real.(fft(Zmag))
        Zmag = sqrt.(Zmag)

        Z = Zmag .* (randn(Nfft) + randn(Nfft) .* 1im)

        y = real.(ifft(Z)) * sqrt(Nfft)

        y = y[1:(N+safety)]
    end

    # Convert fGn to fBn if needed
    if fBn == 1.0
        y = cumsum(y)
    end

    # Resample to match AN model
    y = upsample(y, resamp)

    # Handle mu and sigma
    if mu < 0.5
        sigma = 3.0
    elseif mu < 18.0
        sigma = 30.0
    else
        sigma = 200.0
    end
    y = sigma .* y

    # Return
    return pointer(y[1:nop])
end

"""
    ffGn_native(N)

Synthesize a sample of fractional Gaussian noise.

Adapted from the original MATLAB code avaialble at
https://www.urmc.rochester.edu/labs/carney/publications-code/auditory-models.aspx.

# Arguments
- `N::Int32`: length of the output sequence
- `tdres::Float64`: reciprocal of the sampling rate (1/Hz)
- `Hinput::Float64`: Hurst index. For 0 < H <= 1, we synthesize fractional Gaussian noise with Hurst index H. For 1 < H <= 2, we synthesize fractional Brownian motion with Hurst index H-1.
- `noiseType::Float64`: If noiseType == 0, we return zeros, else we return fractional Gaussian noise
- `mu::Float64`: Mean of the noise
- `sigma::Float64`: The variance of the noise
"""
function ffGn_native(
    N::Int64,
    tdres::Float64,
    Hinput::Float64,
    noiseType::Float64,
    mu::Float64,
    sigma::Float64;
    safety::Int64=4
)
    # Start by handling noiseType
    if noiseType == 0.0
        return zeros(N)
    end

    # If noiseType != 0, we're synthesizing fractional Gaussian noise
    # First, we downsample the number of points
    resamp = Int(ceil(1e-1 / tdres))
    nop = N
    N = Int(ceil(N / resamp) + 1)
    if N < 10
        N = 10
    end

    # Next, determine if fGn or fBm should be produced
    if Hinput < 1.0
        H = Hinput
        fBn = false
    else
        H = Hinput - 1
        fBn = true
    end

    # Calculate fGn
    if H == 0.5
        y = randn(N)
    else
        Nfft = Int(2 ^ ceil(log2(2*(N-1))))
        NfftHalf = Int(round(Nfft / 2))

        k = [0:(NfftHalf-1); NfftHalf:-1:1]
        Zmag = 0.5 * ( (k.+1) .^ (2*H) -2*k .^ (2*H) + abs.(k .- 1) .^ (2*H))

        Zmag = real.(fft(Zmag))
        Zmag = sqrt.(Zmag)

        Z = Zmag .* (randn(Nfft) + randn(Nfft) .* 1im)

        y = real.(ifft(Z)) * sqrt(Nfft)

        y = y[1:(N+safety)]
    end

    # Convert fGn to fBn if needed
    if fBn == 1.0
        y = cumsum(y)
    end

    # Resample to match AN model
    y = upsample(y, resamp)

    y = sigma .* y

    # Return
    return y[1:nop]
end

function ffGn_native(
    N::Int64,
    tdres::Float64,
    Hinput::Float64,
    noiseType::Float64,
    mu::Float64;
)
    if mu == 100.0
        ffGn_native(N, tdres, Hinput, noiseType, mu, 200.0)
    elseif mu == 5.0
        ffGn_native(N, tdres, Hinput, noiseType, mu, 30.0)
    elseif mu == 0.1
        ffGn_native(N, tdres, Hinput, noiseType, mu, 3.0)
    else
        error("The specified mu parameter is not recognized, must be 100, 5, or 0.1!")
    end
end


"""
    decimate(original_signal, k, resamp)

Downsamples a 1D signal of length k by a factor of 1/resamp.

This function is NOT intended to be called by a Julia user on a Julia array. Rather, this
function accepts a pointer to an array and returns a pointer to the output array. It is
intended to be passed as an argument to a function in C where it is called.
"""
function decimate(original_signal::Ptr{Cdouble}, k::Int32, resamp::Int32)
    temp_orig = unsafe_wrap(Array, original_signal, k)
    _resampled = resample(temp_orig, 1/resamp)
    return pointer(_resampled)
end

"""
    upsample(original_signal, resamp)

Upsamples a 1D signal by a factor of resamp.

# Warnings
- This function is (very marginally) affected by a known issue with DSP.resample. The function uses DSP.resample and pads a few zeros to compensate for the few samples gobbled by the upsampling filter (see issue at https://github.com/JuliaDSP/DSP.jl/issues/104).
"""
function upsample(original_signal::Vector{Float64}, resamp::Int64)
    # Directly upsample using DSP.resample
    upsampled_signal = DSP.resample(original_signal, resamp)
    upsampled_signal = [upsampled_signal; zeros(length(original_signal)*resamp - length(upsampled_signal))]
    # Return
    return upsampled_signal
end
