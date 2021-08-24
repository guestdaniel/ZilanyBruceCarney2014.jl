module Hearing
using Statistics

const libihc = "/home/daniel/Hearing/external/libihc.so"

"""
    amplify(signal, dB)

Amplifies a signal's power by a certain amount in dB (or attenuates the signal if dB is negative)
"""
function amplify(signal::Array{Float64, 1}, dB::Float64)::Array{Float64, 1}
    signal .* 10^(dB/20)
end


"""
    dbspl(signal)::Float64

Calculates the dB SPL value of a signal (assuming that the units of the signal are in Pascals)
"""
function dbspl(signal::Array{Float64, 1})::Float64
    20*log10(rms(signal)/20e-6)
end


"""
    pure_tone(freq, phase, dur, fs)

Synthesizes a pure tone with specified frequency, phase offset (in radians), duration, and sampling rate.
"""
function pure_tone(freq::Float64, phase::Float64, dur::Float64, fs::Int64)::Array{Float64, 1}
    sin.(2*pi*freq*LinRange(0, dur, Int64(dur*fs)) .+ phase)
end

function pure_tone(freq::Float64, phase::Irrational, dur::Float64, fs::Int64)::Array{Float64, 1}
    sin.(2*pi*freq*LinRange(0, dur, Int64(dur*fs)) .+ phase)
end


"""
    rms(signal)::Float64

Calculates the room-mean-square (RMS) of a signal
"""
function rms(signal::Array{Float64, 1})::Float64
    sqrt(mean(signal.^2))
end


"""
    scale_dbspl(signal, dB)

Adjusts a signals level to be a certain level in dB SPL
"""
function scale_dbspl(signal::Array{Float64, 1}, dB::Float64)::Array{Float64, 1}
    curr_dB = dbspl(signal)
    delta_dB = dB - curr_dB
    amplify(signal, delta_dB)
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
    Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp, ffGn, decimate)

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
- `ffGn`: function pointer to Julia function that synthesizes fractional Gaussian noise, prepared using `@cfunction` macro
- `decimate`: function pointer to Julia function that performs downsampling, prepared using `@cfunction` macro
"""
function Synapse!(ihcout::Array{Float64, 1}, tdres::Float64, cf::Float64,
                  totalstim::Int32, nrep::Int32, spont::Float64,
                  noiseType::Float64, implnt::Float64, sampFreq::Float64,
                  synouttmp::Array{Float64, 1}, ffGn, decimate)
    ccall((:Synapse, libihc), Cdouble, (Ptr{Cdouble}, Cdouble, Cdouble, Cint,
                                        Cint, Cdouble, Cdouble, Cdouble, Cdouble,
                                        Ptr{Cdouble}, Ptr{nothing}, Ptr{nothing}),
          ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq,
          synouttmp, ffGn, decimate)
end

"""
    SingleAN!(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, varrate, psth, ffGn, decimate)

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
- `ffGn`: function pointer to Julia function that synthesizes fractional Gaussian noise, prepared using `@cfunction` macro
- `decimate`: function pointer to Julia function that performs downsampling, prepared using `@cfunction` macro
"""
function SingleAN(ihcout::Array{Float64, 1}, cf::Float64, nrep::Int32,
                  tdres::Float64, totalstim::Int32, fibertype::Float64,
                  noiseType::Float64, implnt::Float64,
                  meanrate::Array{Float64, 1}, varrate::Array{Float64, 1},
                  psth::Array{Float64, 1}, ffGn, decimate)
    ccall((:SingleAN, libihc), Cvoid, (Ptr{Cdouble}, Cdouble, Cint, Cdouble,
                                       Cint, Cdouble, Cdouble, Cdouble,
                                       Ptr{Cdouble}, Ptr{Cdouble},
                                       Ptr{Cdouble}, Ptr{nothing}, Ptr{nothing}),
          ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt,
          meanrate, varrate, psth, ffGn, decimate)
end


end # module
