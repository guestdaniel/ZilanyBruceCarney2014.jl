module ZBC2014
include("../external/util.jl")

# Declare the location of the shared C library 
const libihc = "/home/daniel/ANF.jl/external/libihc.so"

# Declare some C functions
c_ffGn = @cfunction(ffGn, Vector{Cdouble}, (Cint, ))
c_decimate = @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint))

function translate_species_label(species::String)::Int32
    return Dict([("cat", 1), ("human", 2), ("human_glasberg", 3)])[species]
end


"""
    sim_ihc_zbc2014(input, cf::Float64; fs=10e4, cohc=1.0, cihc=1.0, species="cat")
"""
function sim_ihc_zbc2014(input::Array{Float64, 1}, cf::Float64; fs::Float64=10e4,
                         cohc::Float64=1.0, cihc::Float64=1.0, species::String="cat")
    # Map species string to species integer expected by IHCAN!
    species_flag = translate_species_label(species)
    # Create empty array for output
    output = zeros((length(input), ))
    # Make call
    IHCAN!(input, cf, Int32(1), 1/fs, Int32(length(input)), cohc, cihc, Int32(species_flag), output);
    # Return
    return output
end

"""
    sim_ihc_zbc2014(input, cf::Array{Float64, 1}; fs=10e4, cohc=1.0, cihc=1.0, species="cat")
"""
function sim_ihc_zbc2014(input::Array{Float64, 1}, cf::Array{Float64, 1}; fs::Float64=10e4,
                         cohc::Float64=1.0, cihc::Float64=1.0, species::String="cat")::Array{Float64, 2}
    # Map species string to species integer expected by IHCAN!
    species_flag = translate_species_label(species)
    # Create empty array for output
    _temp_output = zeros((length(input), ))  # pointer to send to C
    output = zeros((length(input), length(cf)))  # actual output
    # Make call
    for ii in 1:length(cf)
        IHCAN!(input, cf[ii], Int32(1), 1/fs, Int32(length(input)), cohc, cihc, Int32(species_flag),
               _temp_output);
        output[:, ii] = _temp_output
    end
    # Return
    return output
end


"""
    IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)

Direct binding to IHCAN C function in model_IHC.c

Passes arguments directly to IHCAN using ccall. Arrays are converted to pointers,
functions are converted to pointers, and all other types are converted directly
to corresponding types in C. Note that while there are type checks enforced
automatically by Julia, there are no sanity checks on any arguments.

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
    ccall((:IHCAN, libihc), Cvoid, (Ptr{Cdouble}, Cdouble, Cint, Cdouble, Cint,
                                    Cdouble, Cdouble, Cint, Ptr{Cdouble}),
          px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
end

"""
    Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp)

Direct binding to Synapse C function in model_Synapse.c

Passes arguments directly to Synapse using ccall. Arrays are converted to pointers,
functions are converted to pointers, and all other types are converted directly
to corresponding types in C. Note that while there are type checks enforced
automatically by Julia, there are no sanity checks on any arguments.

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
    ccall((:Synapse, libihc), Cdouble, (Ptr{Cdouble}, Cdouble, Cdouble, Cint,
                                        Cint, Cdouble, Cdouble, Cdouble, Cdouble,
                                        Ptr{Cdouble}, Ptr{nothing}, Ptr{nothing}),
          ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq,
          synouttmp, c_ffGn, c_decimate)
end

"""
    SingleAN!(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, varrate, psth) 

Direct binding to Synapse C function in model_Synapse.c

Passes arguments directly to Synapse using ccall. Arrays are converted to pointers,
functions are converted to pointers, and all other types are converted directly
to corresponding types in C. Note that while there are type checks enforced
automatically by Julia, there are no sanity checks on any arguments.

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
function SingleAN(ihcout::Array{Float64, 1}, cf::Float64, nrep::Int32,
                  tdres::Float64, totalstim::Int32, fibertype::Float64,
                  noiseType::Float64, implnt::Float64,
                  meanrate::Array{Float64, 1}, varrate::Array{Float64, 1},
                  psth::Array{Float64, 1})
    ccall((:SingleAN, libihc), Cvoid, (Ptr{Cdouble}, Cdouble, Cint, Cdouble,
                                       Cint, Cdouble, Cdouble, Cdouble,
                                       Ptr{Cdouble}, Ptr{Cdouble},
                                       Ptr{Cdouble}, Ptr{nothing}, Ptr{nothing}),
          ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt,
          meanrate, varrate, psth, c_ffGn, c_decimate)
end

end # module

