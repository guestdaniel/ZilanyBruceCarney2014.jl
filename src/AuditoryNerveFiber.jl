module AuditoryNerveFiber

using DSP
using FFTW
#using libzbc2014_jll
#const libzbc2014 = "/home/daniel/AuditoryNerveFiber.jl/external/libzbc2014.so"
const libzbc2014 = joinpath(splitpath(Base.current_project())[1:(end-1)]..., "external", "libzbc2014.so")

export sim_ihc_zbc2014, sim_ihc_zbc2014!,
       sim_synapse_zbc2014, sim_synapse_zbc2014!,
       sim_anrate_zbc2014, sim_anrate_zbc2014!,
       sim_an_zbc2014, sim_spikes_zbc2014

include("utils.jl")

"""
    sim_ihc_zbc2014(input, cf; fs=100e3, cohc=1.0, cihc=1.0, species="human", n_rep=1)

Simulates inner hair cell potential for given acoustic input.

# Arguments
- `input::Vector{Float64}`: sound pressure waveform in pascals
- `cf::Float64`: characteristic frequency of the IHC in Hz
- `fs::Float64`: sampling rate in Hz
- `cohc::Float64`: outer hair cell survival (from 0 to 1)
- `cihc::Float64`: inner hair cell survival (from 0 to 1)
- `species::String`: species, either ("cat" = cat, "human" = humans with Shera tuning, "human_glasberg" = humans with Glasberg tuning)
- `n_rep::Int64`: how many repetitions to perform. Because the inner hair cell stage has no randomness, `n_rep` > 1 simply concatenates `n_rep` copies of IHC simulation and returns 

# Returns
- `output::Vector{Float64}`: inner hair cell potential output, size of `length(input)*n_rep`
"""
function sim_ihc_zbc2014(
    input::AbstractVector{Float64},
    cf::Float64; 
    fs::Float64=100e3,
    cohc::Float64=1.0, 
    cihc::Float64=1.0, 
    species::String="human",
    n_rep::Int64=1
)
    # Create empty array for output
    output = zeros(length(input)*n_rep)

    # Call in-place version function
    sim_ihc_zbc2014!(
        output,
        input,
        cf;
        fs=fs,
        cohc=cohc,
        cihc=cihc,
        species=species,
        n_rep=n_rep
    )

    # Return
    return output
end

function sim_ihc_zbc2014!(
    output::AbstractVector{Float64},
    input::AbstractVector{Float64},
    cf::Float64;
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    n_rep::Int64=1
)
    # Map species string to species integer expected by IHCAN!
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    # Make call to underlying C function
    IHCAN!(
        input,
        cf,
        Int32(n_rep),
        1/fs,
        Int32(length(input)),
        cohc,
        cihc,
        Int32(species_flag),
        output
    )
end

"""
    sim_synapse_zbc2014(input, cf; fs=100e3, fs_synapse=10e3, fiber_type="high", power_law="approximate", fractional=false, n_rep=1)

Simulates synapse output for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, false)
- `n_rep::Int64`: number of repetititons to run. We assume that the input was also generated using `n_rep=n_rep`, hence we infer that the input acoustic waveform is of length `length(input)/n_rep`.

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is `length(input)`
"""
function sim_synapse_zbc2014(
    input::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=100e3,
    fs_synapse::Float64=10e3, 
    fiber_type::String="high", 
    power_law::String="approximate", 
    fractional::Bool=false,
    n_rep::Int64=1
)
    # Create empty array for output
    output = zeros(length(input))

    # Make call to in-place function
    sim_synapse_zbc2014!(
        output,
        input,
        cf;
        fs=fs,
        fs_synapse=fs_synapse,
        fiber_type=fiber_type,
        power_law=power_law,
        fractional=fractional,
        n_rep=n_rep
    )

    # Return
    return output
end

function sim_synapse_zbc2014!(
    output::AbstractVector{Float64},
    input::AbstractVector{Float64},
    cf::Float64;
    fs::Float64=100e3,
    fs_synapse::Float64=10e3,
    fiber_type::String="high",
    power_law::String="approximate",
    fractional::Bool=false,
    n_rep::Int64=1
)
    # Map fiber type string to float code expected by Synapse!
    spont = Dict(
        "low" => 0.1,
        "medium" => 4.0,
        "high" => 100.0
    )[fiber_type]

    # Map power-law implementation type to float code expected by Syanpse!
    implnt = Dict(
        "actual" => 1.0,
        "approximate" => 0.0
    )[power_law]

    # Map fractional to float code expected by Syanpse!
    noiseType = Dict(
        true => 1.0,
        false => 0.0
    )[fractional]

    # Synthesize ffGn
    ffGn = ffGn_native(
        Int(ceil((length(input) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * fs_synapse)),
        1/fs_synapse,
        0.9,
        noiseType,
        spont,
    )

    # Make call to underlying C function
    Synapse!(
        input,
        ffGn,
        1.0/fs,
        cf,
        Int32(length(input)/n_rep),
        Int32(n_rep),
        spont,
        noiseType,
        implnt,
        fs_synapse,
        output
    )

    # Return
    return output
end


"""
    sim_anrate_zbc2014(input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false)

Simulates auditory nerve output (instantaneous firing rate) for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, false)
- `n_rep::Int64`: number of repetititons to run

# Returns
- `::Vector{Float64}`: analytical estimate of instantaneous firing rate, size of `length(input)` (see docs for `sim_an_zbc2014` to understand why)
"""
function sim_anrate_zbc2014(
    input,
    cf::Float64;
    fs::Float64=100e3,
    fiber_type::String="high",
    power_law::String="approximate",
    fractional::Bool=false,
    n_rep::Int64=1
)
    synout = sim_synapse_zbc2014(
        input,
        cf;
        fs=fs,
        fiber_type=fiber_type,
        power_law=power_law,
        fractional=fractional,
        n_rep=n_rep
    )

    return synout ./ (1.0 .+ 0.75e-3 .* synout)
end

function sim_anrate_zbc2014!(
    output::AbstractVector{Float64},
    input::AbstractVector{Float64},
    cf::Float64;
    fs::Float64=100e3,
    fiber_type::String="high",
    power_law::String="approximate",
    fractional::Bool=false,
    n_rep::Int64=1,
)
    sim_synapse_zbc2014!(
        output,
        input,
        cf;
        fs=fs,
        fiber_type=fiber_type,
        power_law=power_law,
        fractional=fractional,
        n_rep=n_rep
    )

    output .= output ./ (1.0 .+ 0.75e-3 .* output)
    return output
end

"""
    sim_an_zbc2014(input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

Simulates auditory nerve output (spikes and firing rate) for a given inner hair cell input

# Arguments
- `input::Vector{Float64}`: input hair cell potential (from sim_ihc_zbc2014)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, talse)
- `n_rep::Int64`: number of repetititons to run. We assume that the input was also generated using `n_rep=n_rep`, hence we infer that the input acoustic waveform is of length `length(input)/n_rep`.

# Returns
- `meanrate::Vector{Float64}`: analytical estimate of instantaneous firing rate, of size `length(input)/n_rep`
- `varrrate::Vector{Float64}`: analytical estimate of instantaneous firing rate variance, of size `length(input)/n_rep`
- `psth::Vector{Float64}`: peri-stimulus time histogram, of size `length(input)/n_rep`

# Notes
- The behavior of `n_rep` in `sim_an_zbc2014` and derived functions differs slightly from 
`sim_synapse_zbc2014`. `sim_synapse_zbc2014` will return an array of length 
`length(input)`, and will contain a single response to `n_rep` concatenated copies 
of the inner hair cell potential to a single copy of the underlying acoustic stimulus 
(assuming that `sim_ihc_zbc2014` was also evaluated using
`n_rep=n_rep`). `sim_an_zbc2014` will return an array of length `length(input)`. Underneath,
it generates the same synapse response as `sim_synapse_zbc2014`, but then averages over
the repetitions and averages the spike train into a peristimulus time histogram (PSTH).

# TODO
- Modify this function to match style of previous functions
"""
function sim_an_zbc2014(
    input::AbstractVector{Float64},
    cf::Float64; 
    fs::Float64=100e3,
    fiber_type::String="high", 
    power_law::String="approximate",
    fractional::Bool=false, 
    n_rep::Int64=1
)
    # Calculate totalstim based on size of input
    totalstim = Int64(length(input)/n_rep)
    # Create empty array for output
    meanrate = zeros(totalstim)
    varrate = zeros(totalstim)
    psth = zeros(totalstim)
    # Make call
    sim_an_zbc2014!(meanrate, varrate, psth, input, cf; fs=fs, fiber_type=fiber_type, power_law=power_law, fractional=fractional, n_rep=n_rep)
    return (meanrate, varrate, psth)
end

"""
    sim_an_zbc2014!(meanrrate, varrate, psth, input, cf; kwargs...)

Same as `sim_an_zbc2014` except modifies first three input args in-place
"""
function sim_an_zbc2014!(
    meanrate::AbstractVector{Float64},
    varrate::AbstractVector{Float64},
    psth::AbstractVector{Float64},
    input::AbstractVector{Float64},
    cf::Float64;
    fs::Float64=100e3,
    fiber_type::String="high",
    power_law::String="approximate",
    fractional::Bool=false,
    n_rep::Int64=1
)
    # Calculate totalstim based on size of input
    totalstim = Int64(length(input)/n_rep)
    # Map fiber type string to float code expected by Synapse!
    spont = Dict([("low", 0.1), ("medium", 4.0), ("high", 100.0)])[fiber_type]
    # Map fiber type string to float code expected by Synapse!
    fibertype = Dict([("low", 1.0), ("medium", 2.0), ("high", 3.0)])[fiber_type]
    # Map power-law implementation type to float code expected by Syanpse!
    implnt = Dict([("actual", 1.0), ("approximate", 0.0)])[power_law]
    # Map fractional to float code expected by Syanpse!
    noiseType = Dict([(true, 1.0), (false, 0.0)])[fractional]
    # Synthesize ffGn
    ffGn = ffGn_native(
        Int(ceil((length(input) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)),
        1/10e3,
        0.9,
        noiseType,
        spont,
    )
    # Make call
    SingleAN!(input, ffGn, cf, Int32(n_rep), 1.0/fs, Int32(totalstim), fibertype, noiseType, implnt,
              meanrate, varrate, psth)
    return (meanrate, varrate, psth)
end

sim_spikes_zbc2014(args...; kwargs...) = sim_an_zbc2014(args...; kwargs...)[3]

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
    ffGn::Vector{Float64},
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
                Ptr{Cdouble}, # ffGn
                Cdouble,      # tdres
                Cdouble,      # cf
                Cint,         # totalstim
                Cint,         # nrep
                Cdouble,      # spont
                Cdouble,      # noiseType
                Cdouble,      # implnt
                Cdouble,      # sampFreq
                Ptr{Cdouble}, # synouttmp
                Ptr{nothing}  # decimate
            ),
            ihcout, ffGn, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp,   # input args
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
    ffGn::Vector{Float64},
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
                Ptr{Cdouble},  # ffGn
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
                Ptr{nothing},  # decimate
                Ptr{nothing}   # random_numbers
            ),
            ihcout, ffGn, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, varrate, psth,  # input args
            @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),                         # input cfunction
            @cfunction(random_numbers, Ptr{Cdouble}, (Cint, ))                                          # input cfunction
        )        
end

end # module
