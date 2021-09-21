module AuditoryNerveFiber

# Handle imports
using DSP
using FFTW
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using libzbc2014_jll


"""
    random_numbers(length)

Generates random numbers like MATLAB's rand or Python's numpy.random.rand
"""
function random_numbers(n::Int32)
    return pointer(rand(Float64, n))
end


"""
    ffGn(N)

Return an array of zeros (as a standin for a real ffGn function).

# Warnings
- Note that this function currently just returns zeros, where it should return fractional Gaussian noise
"""
function ffGn(N::Int32)
    return zeros((N, ))
end


""" 
    ffGn(N)

Synthesize a sample of fractional Gaussian noise.

This is a direct translation of Python code written by Marek Rudnicki in the cochlea package
(https://github.com/mrkrd/cochlea).
"""
function ffGn(N::Int32, tdres::Float64, Hinput::Float64, noiseType::Float64, mu::Float64; safety::Int64=4)
    # Start by handling noiseType
    if noiseType == 0
        return pointer(zeros(N))
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
    upsample(original_signal, resamp)

Upsamples a 1D signal by a factor of resamp. 

# Warnings
- This function is (very marginally) affected by a known issue with DSP.resample. The function uses DSP.resample and pads a few zeros to compensate for the few samples gobbled by the upsampling filter (see issue at https://github.com/JuliaDSP/DSP.jl/issues/104).
"""
function upsample(original_signal::Array{Float64, 1}, resamp::Int64)
    # Directly upsample using DSP.resample
    upsampled_signal = DSP.resample(original_signal, resamp)
    upsampled_signal = [upsampled_signal; zeros(length(original_signal)*resamp - length(upsampled_signal))]
    # Return
    return upsampled_signal
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
    sim_ihc_zbc2014(input, cf; fs=10e4, cohc=1.0, cihc=1.0, species="cat")

Simulates inner hair cell potential for given acoustic input.

# Arguments
- `input::Array{Float64, 1}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the IHC in Hz
- `fs::Float64`: sampling rate in Hz
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::Int32`: species, either (1 = cat, 2 = humans with Shera tuning, 3 = humans with Glasberg tuning)

# Returns
- `output::Array{Float64, 1}`: inner hair cell potential output
"""
function sim_ihc_zbc2014(input::Array{Float64, 1}, cf::Float64; fs::Float64=10e4,
                         cohc::Float64=1.0, cihc::Float64=1.0, species::String="cat")
    # Map species string to species integer expected by IHCAN!
    species_flag = Dict([("cat", 1), ("human", 2), ("human_glasberg", 3)])[species]
    # Create empty array for output
    output = zeros((length(input), ))
    # Make call
    IHCAN!(input, cf, Int32(1), 1/fs, Int32(length(input)), cohc, cihc, Int32(species_flag), 
           output);
    # Return
    return output
end


"""
    sim_synapse_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, cohc=1.0, cihc=1.0)

Simulates synapse output for a given inner hair cell input

# Arguments
- `input::Array{Float64, 1}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `frac_noise::String`: controls whether we use true or approximate fractional Gaussian noise implementation, one of ("actual", "approximate")

# Returns
- `output::Array{Float64, 1}`: synapse output (unknown units?)
"""
function sim_synapse_zbc2014(input::Array{Float64, 1}, cf::Float64; fs::Float64=10e4,
                             fs_synapse::Float64=10e3, fiber_type::String="high", 
                             frac_noise::String="approximate", noise_type::String="ffGn")
    # Map fiber type string to float code expected by Synapse!
    spont = Dict([("low", 0.1), ("medium", 4.0), ("high", 100.0)])[fiber_type]
    # Map fractional noise implementation type to float code expected by Syanpse!
    implnt = Dict([("actual", 1.0), ("approximate", 0.0)])[frac_noise]
    # Map noise type to float code expected by Syanpse!
    noiseType = Dict([("ffGn", 1.0), ("Gaussian", 0.0)])[noise_type]
    # Create empty array for output
    output = zeros((length(input), ))
    # Make call
    Synapse!(input, 1.0/fs, cf, Int32(length(input)), Int32(1), spont, noiseType, implnt, 
             fs_synapse, output)
    # Return
    return output
end


"""
    sim_an_zbc2014(input, cf; fs=10e4, fiber_type="high", cohc=1.0, cihc=1.0)

Simulates auditory nerve output (spikes or firing rate) for a given inner hair cell input

# Arguments
- `input::Array{Float64, 1}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `frac_noise::String`: controls whether we use true or approximate fractional Gaussian noise implementation, one of ("actual", "approximate")

# Returns
- `meanrate::Array{Float64, 1}`: analytical estimate of instantaneous firing rate
- `varrrate::Array{Float64, 1}`: analytical estimate of instantaneous firing rate variance
- `psth::Array{Float64, 1}`: peri-stimulus time histogram (NOT IMPLEMENTED, SHOULD BE EMPTY)
"""
function sim_an_zbc2014(input::Array{Float64, 1}, cf::Float64; fs::Float64=10e4,
                        fiber_type::String="high", frac_noise::String="approximate",
                        noise_type::String="ffGn")

    # Map fiber type string to float code expected by Synapse!
    fibertype = Dict([("low", 1.0), ("medium", 2.0), ("high", 3.0)])[fiber_type]
    # Map fractional noise implementation type to float code expected by Syanpse!
    implnt = Dict([("actual", 1.0), ("approximate", 0.0)])[frac_noise]
    # Map noise type to float code expected by Syanpse!
    noiseType = Dict([("ffGn", 1.0), ("Gaussian", 0.0)])[noise_type]
    # Create empty array for output
    meanrate = zeros((length(input), ))
    varrate = zeros((length(input), ))
    psth = zeros((length(input), ))
    # Make call
    SingleAN!(input, cf, Int32(1), 1.0/fs, Int32(length(input)), fibertype, noiseType, implnt, 
              meanrate, varrate, psth)
    # Return
    return (meanrate, varrate, psth)
end



"""
    IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)

Direct binding to IHCAN C function in model_IHC.c

Passes arguments directly to IHCAN using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `px::Array{Float64, 1}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `nrep::Int32`: number of repetitions to simulate. Note that for the IHC simulation, one "true" simulation is conducted and then that simulation is copied and tiled (because there is no randomness in the IHC simulation) to simulate multiple times.
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `totalstim::Int32`: number of samples in simulation
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::Int32`: species, either (1 = cat, 2 = humans with Shera tuning, 3 = humans with Glasberg tuning)
- `ihcout::Array{Float64, 1}`: array of same size as `px`, used to store output from C
"""
function IHCAN!(px::Array{Float64, 1}, cf::Float64, nrep::Int32, tdres::Float64,
                totalstim::Int32, cohc::Float64, cihc::Float64, species::Int32,
                ihcout::Array{Float64, 1})
    ccall((:IHCAN, libzbc2014), Cvoid, (Ptr{Cdouble}, Cdouble, Cint, Cdouble, Cint,
                                    Cdouble, Cdouble, Cint, Ptr{Cdouble}),
          px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
end


"""
    Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, 
             synouttmp)

Direct binding to Synapse C function in model_Synapse.c

Passes arguments directly to Synapse using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `ihcout::Array{Float64, 1}`: output from IHC simulation (`IHCAN!`)
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `totalstim::Int32`: number of samples in simulation
- `nrep::Int32`: number of repetitions to simulate.
- `spont::Float64`: spontaneous rate, either (0.1 == low spont fiber, 4.0 == medium spont fiber, 100.0 == high spont fiber)
- `noiseType::Float64`: NOT CURRENTLY IMPLEMENTED
- `implnt::Float64`: whether or not to use exact implementation of fractional Gaussian noise, either (1.0 == use, 0.0 == approximate)
- `sampFreq::Float64`: sampling frequency of the power law stage in Hz. Simulations are decimated to sampFreq from 1/tdres before the power law stage and then upsampled back to the original sampling rate. The product of tdres and sampFreq, which indicates the amount to decimate by, must be an integer
- `synouttmp::Array{Float64, 1}`: array of same size as `ihcout`, used to store output from C
"""
function Synapse!(ihcout::Array{Float64, 1}, tdres::Float64, cf::Float64,
                  totalstim::Int32, nrep::Int32, spont::Float64,
                  noiseType::Float64, implnt::Float64, sampFreq::Float64,
                  synouttmp::Array{Float64, 1})
    ccall((:Synapse, libzbc2014), Cdouble, 
          (Ptr{Cdouble}, # ihcout
           Cdouble,      # tdres
           Cdouble,      # cf
           Cint,         # totalstim
           Cint,         # nrep
           Cdouble,      # spont
           Cdouble,      # noiseType
           Cdouble,      # implnt
           Cdouble,      # sampFreq
           Ptr{Cdouble}, # synouttmp
           Ptr{nothing}, # ffGn
           Ptr{nothing}),# decimate
          ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp, 
          @cfunction(ffGn, Ptr{Cdouble}, (Cint, Cdouble, Cdouble, Cdouble, Cdouble)), 
          @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)))
end


"""
    SingleAN!(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, 
              varrate, psth)

Direct binding to Synapse C function in model_Synapse.c

Passes arguments directly to Synapse using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `ihcout::Array{Float64, 1}`: output from IHC simulation (`IHCAN!`)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `nrep::Int32`: number of repetitions to simulate.
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `totalstim::Int32`: number of samples in simulation
- `fibertype::Float64`: fiber type, either (1.0 == low, 2.0 == med, 3.0 == high)
- `noiseType::Float64`: NOT CURRENTLY IMPLEMENTED
- `implnt::Float64`: whether or not to use exact implementation of fractional Gaussian noise, either (1.0 == use, 0.0 == approximate)
- `meanrate::Array{Float64, 1}`: array of same size as `ihcout`, used to store analytical firing rate output
- `varrate::Array{Float64, 1}`: array of same size as `ihcout`, used to store analytical firing rate variance output
- `psth::Array{Float64, 1}`: array of same size as `ihcout`, used to store empirical PSTH output
"""
function SingleAN!(ihcout::Array{Float64, 1}, cf::Float64, nrep::Int32,
                   tdres::Float64, totalstim::Int32, fibertype::Float64,
                   noiseType::Float64, implnt::Float64,
                   meanrate::Array{Float64, 1}, varrate::Array{Float64, 1},
                   psth::Array{Float64, 1})
    ccall((:SingleAN, libzbc2014), Cvoid, 
          (Ptr{Cdouble},  # ihcout
           Cdouble,       # cf
           Cint,          # nrep
           Cdouble,       # tdres
           Cint,          # totalstim
           Cdouble,       # fibertype
           Cdouble,       # noiseType
           Cdouble,       # implnt
           Ptr{Cdouble},  # meanrate
           Ptr{Cdouble},  # varrate
           Ptr{Cdouble},  # psth
           Ptr{nothing},  # ffGn
           Ptr{nothing},  # decimate
           Ptr{nothing}), # random_numbers
          ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt,
          meanrate, varrate, psth, 
          @cfunction(ffGn, Ptr{Cdouble}, (Cint, Cdouble, Cdouble, Cdouble, Cdouble)), 
          @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
          @cfunction(random_numbers, Ptr{Cdouble}, (Cint, )))
end

end # module

