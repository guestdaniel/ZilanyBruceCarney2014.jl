module AuditoryNerveFiber

using DSP
using FFTW
const libzbc2014 = "/home/oxenham0/guest121/thesis/AuditoryNerveFiber.jl/external/libzbc2014.so"

export sim_ihc_zbc2014, sim_synapse_zbc2014, sim_an_zbc2014, sim_anrate_zbc2014, sim_bm_zbc2014, sim_spikes_zbc2014


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
    ffGn(N)

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
function ffGn(
    N::Int32, 
    tdres::Float64, 
    Hinput::Float64, 
    noiseType::Float64, 
    mu::Float64; 
    safety::Int64=4
)
    # Start by handling noiseType
    if noiseType == 0.0
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


"""
    dispatch_vectorized_input(func)

Macro that defines a method for model functions dispatched over a vector of inputs 

# Arguments
- `func`: A function implementing an "input" the first position arg, and "cf" as the second, otherwise accepting kwargs
"""

macro dispatch_vectorized_input(func)
    :( $func(input::Vector{<:Vector{Float64}}, cf::Float64; kwargs...) = map(signal -> $func(signal, cf; kwargs...), input) )
end

"""
    dispatch_vectorized_cfs(func)

Macro that defines a method for model functions dispatched over a vector of cfs

# Arguments
- `func`: A function implementing an "input" the first position arg, and "cf" as the second, otherwise accepting kwargs
"""

macro dispatch_vectorized_cfs(func)
    :( $func(input::Vector{Float64}, cf::Vector{Float64}; kwargs...) = transpose(hcat(map(x -> $func(input, x; kwargs...), cf)...)) )
end

"""
    dispatch_vectorized_input_and_cfs(func)

Macro that defines a method for model functions dispatched over a vector of inputs and a vector of cfs

# Arguments
- `func`: A function implementing an "input" the first position arg, and "cf" as the second, otherwise accepting kwargs
"""

macro dispatch_vectorized_input_and_cfs(func)
    :( $func(input::Vector{<:Vector{Float64}}, cf::Vector{Float64}; kwargs...) = map(signal -> transpose(hcat(map(x -> $func(signal, x; kwargs...), cf)...)), input) )
end

"""
    dispatch_matrix_input(func)

Macro that defines a method for model functions dispatched over a matrix inputs and a vector of cfs

# Arguments
- `func`: A function implementing an "input" the first position arg, and "cf" as the second, otherwise accepting kwargs
"""

macro dispatch_matrix_input(func)
    :( 
    function $func(input::AbstractMatrix{Float64}, cf::Vector{Float64}; kwargs...) 
        n_cf = size(input)[1]
        output = Vector{Vector{Float64}}()
        for row in 1:size(input)[1]
            push!(output, $func(input[row, :], cf[row]; kwargs...))
        end
        return permutedims(hcat(output...))
    end
    )
end

"""
    dispatch_vector_of_matrix_input(func)

Macro that defines a method for model functions dispatched over a vector of matrix inputs and a vector of cfs

# Arguments
- `func`: A function implementing an "input" the first position arg, and "cf" as the second, otherwise accepting kwargs
"""

macro dispatch_vector_of_matrix_input(func)
    :( 
    function $func(input::Vector{<:AbstractMatrix{Float64}}, cf::Vector{Float64}; kwargs...) 
        map(neurogram -> $func(neurogram, cf; kwargs...), input)
    end
    )
end


"""
    sim_bm_zbc2014(input, cf; fs=10e4, cohc=1.0, cihc=1.0, species="human", n_rep=1)

Simulates basilar membrane vibration for given acoustic input.

# Arguments
- `input::Vector{Float64}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the IHC in Hz
- `fs::Float64`: sampling rate in Hz
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::String`: species, either ("cat" = cat, "human" = humans with Shera tuning, "human_glasberg" = humans with Glasberg tuning)

# Returns
- `output::Vector{Float64}`: basilar membrane output, same size as input
"""
function sim_bm_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=10e4,
    cohc::Float64=1.0, 
    cihc::Float64=1.0, 
    species::String="human",
)
    # Map species string to species integer expected by IHCAN!
    species_flag = Dict([("cat", 1), ("human", 2), ("human_glasberg", 3)])[species]
    # Create empty array for output
    output_ihc = zeros((length(input), ))
    output_bm = zeros((length(input), ))
    # Make call
    BM!(input, cf, Int32(1), 1/fs, Int32(length(input)), cohc, cihc, Int32(species_flag), output_ihc, output_bm)
    # Return
    return output_bm
end



"""
    sim_ihc_zbc2014(input, cf; fs=10e4, cohc=1.0, cihc=1.0, species="human", n_rep=1)

Simulates inner hair cell potential for given acoustic input.

# Arguments
- `input::Vector{Float64}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the IHC in Hz
- `fs::Float64`: sampling rate in Hz
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::String`: species, either ("cat" = cat, "human" = humans with Shera tuning, "human_glasberg" = humans with Glasberg tuning)
- `n_rep::Int64`: how many repetitions to perform. Because the inner hair cell stage has no randomness, n_rep > 1 simply concatenates n_rep copies of IHC simulation and returns 

# Returns
- `output::Vector{Float64}`: inner hair cell potential output, size of length(input)*n_rep
"""
function sim_ihc_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=10e4,
    cohc::Float64=1.0, 
    cihc::Float64=1.0, 
    species::String="human",
    n_rep::Int64=1
)
    # Map species string to species integer expected by IHCAN!
    species_flag = Dict([("cat", 1), ("human", 2), ("human_glasberg", 3)])[species]
    # Create empty array for output
    output = zeros((length(input)*n_rep, ))
    # Make call
    IHCAN!(input, cf, Int32(n_rep), 1/fs, Int32(length(input)), cohc, cihc, Int32(species_flag), output);
    # Return
    return output
end

@dispatch_vectorized_input(sim_ihc_zbc2014)
@dispatch_vectorized_cfs(sim_ihc_zbc2014)
@dispatch_vectorized_input_and_cfs(sim_ihc_zbc2014)
@dispatch_matrix_input(sim_ihc_zbc2014)
#@dispatch_vector_of_matrix_input(sim_ihc_zbc2014)


"""
    sim_synapse_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, fiber_type="high", power_law="approximate", fractional=false, n_rep=1)

Simulates synapse output for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, false)
- `n_rep::Int64`: number of repetititons to run. We assume that the input was also generated using n_rep=n_rep, hence we infer that the input acoustic waveform is of length length(input)/n_rep.

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is length(input)*n_rep
"""
function sim_synapse_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=10e4,
    fs_synapse::Float64=10e3, 
    fiber_type::String="high", 
    power_law::String="approximate", 
    fractional::Bool=false,
    n_rep::Int64=1
)
    # Map fiber type string to float code expected by Synapse!
    spont = Dict([("low", 0.1), ("medium", 4.0), ("high", 100.0)])[fiber_type]
    # Map power-law implementation type to float code expected by Syanpse!
    implnt = Dict([("actual", 1.0), ("approximate", 0.0)])[power_law]
    # Map fractional to float code expected by Syanpse!
    noiseType = Dict([(true, 1.0), (false, 0.0)])[fractional]
    # Create empty array for output
    output = zeros((length(input), ))
    # Make call
    Synapse!(input, 1.0/fs, cf, Int32(length(input)/n_rep), Int32(n_rep), spont, noiseType, implnt, fs_synapse, output)
    # Return
    return output
end

@dispatch_vectorized_input(sim_synapse_zbc2014)
@dispatch_vectorized_cfs(sim_synapse_zbc2014)
@dispatch_vectorized_input_and_cfs(sim_synapse_zbc2014)
@dispatch_matrix_input(sim_synapse_zbc2014)
#@dispatch_vector_of_matrix_input(sim_synapse_zbc2014)


"""
    sim_an_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

Simulates auditory nerve output (spikes and firing rate) for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, talse)
- `n_rep::Int64`: number of repetititons to run. We assume that the input was also generated using n_rep=n_rep, hence we infer that the input acoustic waveform is of length length(input)/n_rep.

# Returns
- `meanrate::Vector{Float64}`: analytical estimate of instantaneous firing rate, of size length(input)
- `varrrate::Vector{Float64}`: analytical estimate of instantaneous firing rate variance, of size length(input)
- `psth::Vector{Float64}`: peri-stimulus time histogram, of size length(input)

# Notes
- The behavior of `n_rep` in sim_an_zbc2014 and derived functions differs slightly from 
sim_synapse_zbc2014. sim_synapse_zbc2014 will return an array of length 
`length(input) * n_rep`, and will contain a single response to n_rep concatenated copies 
of the inner hair cell potential (assuming that sim_ihc_zbc2014 was also evaluated using
n_rep=n_rep). sim_an_zbc2014 will return an array of length `length(input)`. Underneath,
it generates the same synapse response as `sim_synapse_zbc2014`, but then averages over
the repetitions and averages the spike train into a peristimulus time histogram (PSTH).  
"""
function sim_an_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=10e4,
    fiber_type::String="high", 
    power_law::String="approximate",
    fractional::Bool=false, 
    n_rep::Int64=1
)
    # Calculate totalstim based on size of input
    totalstim = Int64(length(input)/n_rep)
    # Map fiber type string to float code expected by Synapse!
    fibertype = Dict([("low", 1.0), ("medium", 2.0), ("high", 3.0)])[fiber_type]
    # Map power-law implementation type to float code expected by Syanpse!
    implnt = Dict([("actual", 1.0), ("approximate", 0.0)])[power_law]
    # Map fractional to float code expected by Syanpse!
    noiseType = Dict([(true, 1.0), (false, 0.0)])[fractional]
    # Create empty array for output
    meanrate = zeros(totalstim)
    varrate = zeros(totalstim)
    psth = zeros(totalstim)
    # Make call
    SingleAN!(input, cf, Int32(n_rep), 1.0/fs, Int32(totalstim), fibertype, noiseType, implnt, 
              meanrate, varrate, psth)
    return (meanrate, varrate, psth)
end


"""
    sim_spikes_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

Simulates auditory nerve output (spikes only) for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, talse)
- `n_rep::Int64`: number of repetititons to run 

# Returns
- `psth::Vector{Float64}`: peri-stimulus time histogram, size of length(input) (see docs for sim_an_zbc2014 to understand why)
"""
function sim_spikes_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=10e4,
    fiber_type::String="high", 
    power_law::String="approximate",
    fractional::Bool=false, 
    n_rep::Int64=1
)
    sim_an_zbc2014(input, cf; fs=fs, fiber_type=fiber_type, power_law=power_law, fractional=fractional, n_rep=n_rep)[3]
end

@dispatch_vectorized_input(sim_spikes_zbc2014)
@dispatch_vectorized_cfs(sim_spikes_zbc2014)
@dispatch_vectorized_input_and_cfs(sim_spikes_zbc2014)
@dispatch_matrix_input(sim_spikes_zbc2014)
#@dispatch_vector_of_matrix_input(sim_spikes_zbc2014)


"""
    sim_anrate_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, power_law="approximate", fractional=false)

Simulates auditory nerve output (instantaneous firing rate) for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, talse)
- `n_rep::Int64`: number of repetititons to run 

# Returns
- `psth::Vector{Float64}`: analytical estimate of instantaneous firing rate, size of length(input) (see docs for sim_an_zbc2014 to understand why)
"""
function sim_anrate_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=10e4,
    fiber_type::String="high", 
    power_law::String="approximate",
    fractional::Bool=false, 
    n_rep::Int64=1
)
    sim_an_zbc2014(input, cf; fs=fs, fiber_type=fiber_type, power_law=power_law, fractional=fractional, n_rep=n_rep)[1]
end

@dispatch_vectorized_input(sim_anrate_zbc2014)
@dispatch_vectorized_cfs(sim_anrate_zbc2014)
@dispatch_vectorized_input_and_cfs(sim_anrate_zbc2014)
@dispatch_matrix_input(sim_anrate_zbc2014)
#@dispatch_vector_of_matrix_input(sim_anrate_zbc2014)


"""
    BM!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout, bmout)

Direct binding to BM C function in model_IHC.c

Passes arguments directly to BM using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `px::Vector{Float64}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `nrep::Int32`: number of repetitions to simulate. Note that for the IHC simulation, one "true" simulation is conducted and then that simulation is copied and tiled (because there is no randomness in the IHC simulation) to simulate multiple times.
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `totalstim::Int32`: number of samples in simulation
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::Int32`: species, either (1 = cat, 2 = humans with Shera tuning, 3 = humans with Glasberg tuning)
- `ihcout::Vector{Float64}`: array of same size as `px`, used to store IHC output from C
- `bmout::Vector{Float64}`: array of same size as `px`, used to store BM output from C
"""
function BM!(
    px::Vector{Float64}, 
    cf::Float64, 
    nrep::Int32, 
    tdres::Float64,
    totalstim::Int32, 
    cohc::Float64, 
    cihc::Float64, 
    species::Int32,
    ihcout::Vector{Float64},
    bmout::Vector{Float64},
)
    ccall(
            (:BM, libzbc2014),       # function call
            Cvoid,                   # return type
            (                        # arg types
                Ptr{Cdouble},        # px
                Cdouble,             # cf
                Cint,                # nrep
                Cdouble,             # tdres
                Cint,                # totalstim
                Cdouble,             # cohc
                Cdouble,             # cihc
                Cint,                # species
                Ptr{Cdouble},        # ihcout
                Ptr{Cdouble},        # bmout
            ),
            px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout, bmout  # pass arguments
        )
end



"""
    IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)

Direct binding to IHCAN C function in model_IHC.c

Passes arguments directly to IHCAN using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `px::Vector{Float64}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `nrep::Int32`: number of repetitions to simulate. Note that for the IHC simulation, one "true" simulation is conducted and then that simulation is copied and tiled (because there is no randomness in the IHC simulation) to simulate multiple times.
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `totalstim::Int32`: number of samples in simulation
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::Int32`: species, either (1 = cat, 2 = humans with Shera tuning, 3 = humans with Glasberg tuning)
- `ihcout::Vector{Float64}`: array of same size as `px`, used to store output from C
"""
function IHCAN!(
    px::Vector{Float64}, 
    cf::Float64, 
    nrep::Int32, 
    tdres::Float64,
    totalstim::Int32, 
    cohc::Float64, 
    cihc::Float64, 
    species::Int32,
    ihcout::Vector{Float64}
)
    ccall(
            (:IHCAN, libzbc2014),    # function call
            Cvoid,                   # return type
            (                        # arg types
                Ptr{Cdouble},        # px
                Cdouble,             # cf
                Cint,                # nrep
                Cdouble,             # tdres
                Cint,                # totalstim
                Cdouble,             # cohc
                Cdouble,             # cihc
                Cint,                # species
                Ptr{Cdouble}         # ihcout
            ),
            px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout  # pass arguments
        )
end


"""
    Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp)

Direct binding to Synapse C function in model_Synapse.c

Passes arguments directly to Synapse using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `ihcout::Vector{Float64}`: output from IHC simulation (`IHCAN!`)
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `totalstim::Int32`: number of samples in simulation
- `nrep::Int32`: number of repetitions to simulate.
- `spont::Float64`: spontaneous rate, either (0.1 == low spont fiber, 4.0 == medium spont fiber, 100.0 == high spont fiber)
- `noiseType::Float64`: whether we use ffGn or not (1.0 == ffGn, 0.0 == not)
- `implnt::Float64`: whether or not to use exact implementation of power-law adaptation, either (1.0 == actual, 0.0 == approximate)
- `sampFreq::Float64`: sampling frequency of the power law stage in Hz. Simulations are decimated to sampFreq from 1/tdres before the power law stage and then upsampled back to the original sampling rate. The product of tdres and sampFreq, which indicates the amount to decimate by, must be an integer
- `synouttmp::Vector{Float64}`: array of same size as `ihcout`, used to store output from C
"""
function Synapse!(
    ihcout::Vector{Float64}, 
    tdres::Float64, 
    cf::Float64,
    totalstim::Int32, 
    nrep::Int32, 
    spont::Float64,
    noiseType::Float64, 
    implnt::Float64, 
    sampFreq::Float64,
    synouttmp::Vector{Float64}
)
    ccall(
            (:Synapse, libzbc2014),  # function call
            Cdouble,                 # return type
            (
                Ptr{Cdouble}, # ihcout
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
                Ptr{nothing}  # decimate_fft
            ),
            ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp,   # input args
            @cfunction(ffGn, Ptr{Cdouble}, (Cint, Cdouble, Cdouble, Cdouble, Cdouble)),          # input cfunction
            @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint))                   # input cfunction
        )
end


"""
    SingleAN!(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, varrate, psth)

Direct binding to Synapse C function in model_Synapse.c

Passes arguments directly to Synapse using ccall. Arrays are converted to pointers, 
functions are converted to pointers, and all other types are converted directly to 
corresponding types in C. Note that while there are type checks enforced automatically by 
Julia, there are no sanity checks on any arguments.

# Arguments
- `ihcout::Vector{Float64}`: output from IHC simulation (`IHCAN!`)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `nrep::Int32`: number of repetitions to simulate.
- `tdres::Float64`: time-domain resolution (i.e., reciprocal of sampling rate)
- `totalstim::Int32`: number of samples in simulation
- `fibertype::Float64`: fiber type, either (1.0 == low, 2.0 == med, 3.0 == high)
- `noiseType::Float64`: whether we use ffGn or not (1.0 == ffGn, 0.0 == not)
- `implnt::Float64`: whether or not to use exact implementation of power-law adaptation, either (1.0 == actual, 0.0 == approximate)
- `meanrate::Vector{Float64}`: array of same size as `ihcout`, used to store analytical firing rate output
- `varrate::Vector{Float64}`: array of same size as `ihcout`, used to store analytical firing rate variance output
- `psth::Vector{Float64}`: array of same size as `ihcout`, used to store empirical PSTH output
"""
function SingleAN!(
    ihcout::Vector{Float64}, 
    cf::Float64, 
    nrep::Int32,
    tdres::Float64, 
    totalstim::Int32, 
    fibertype::Float64,
    noiseType::Float64, 
    implnt::Float64,
    meanrate::Vector{Float64}, 
    varrate::Vector{Float64},
    psth::Vector{Float64}
)
    ccall(
            (:SingleAN, libzbc2014),  # function call
            Cvoid,                    # return type
            (
                Ptr{Cdouble},  # ihcout
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
                Ptr{nothing}   # random_numbers
            ),
            ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, varrate, psth,  # input args
            @cfunction(ffGn, Ptr{Cdouble}, (Cint, Cdouble, Cdouble, Cdouble, Cdouble)),                 # input cfunction
            @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),                         # input cfunction
            @cfunction(random_numbers, Ptr{Cdouble}, (Cint, ))                                          # input cfunction
        )        
end

end # module

