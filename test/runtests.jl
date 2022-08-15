"""
This testing suite tests tools and AN models implemented in ANF.jl. It makes
sure that models can be evaluated without error, makes sure arguments are
handled correctly and return correct output shapes, and it makes sure that
response properties are correct at a coarse level (exact response properties
are tested [visually] elsewhere).
"""

using Test
using AuditoryNerveFiber
using AuditorySignalUtils
using Statistics
using DSP

# Declare various constants that hold across all tests in this file
fs = 100e3
dur = 0.1
freq = 1000.0
pt = scale_dbspl(pure_tone(freq, 0.0, dur, fs), 50.0)
pt_up = scale_dbspl(pure_tone(freq, 0.0, dur, fs*5), 50.0)[1:50000]
tol = 1e-3  # absolute tolerance on comparisons of approximate equality

# Test ffGn
@testset "Fractional Gaussian noise" begin
    # Test that we can call the function
    @test begin
        sample = AuditoryNerveFiber.ffGn_native(10000, 1/fs, 0.75, 0.0, 1.0)
        isapprox(sample, zeros(10000))
    end
    # Test that the noiseType switch behaves as expect
    @test begin
        sample = AuditoryNerveFiber.ffGn_native(10000, 1/fs, 0.75, 1.0, 1.0)
        var(sample) > 0.0
    end
    # Test that raising the mean to the branch points in the code (0.5, 18.0
    # results in corresponding changes in output variance
    @test begin
        myfunc(mu) = std(AuditoryNerveFiber.ffGn_native(100000, 1/fs, 0.75, 1.0, mu))
        vars = map(myfunc, [0.4, 10.0, 20.0])
        vars[1] < vars[2] < vars[3]
    end
end

# Test upsampling function by upsampling pure tone and then verifying that points
# corresponding to original sample points match
@test begin
    pt_upsampled = AuditoryNerveFiber.upsample(pt, 5)
    maximum(abs.(pt_upsampled - pt_up)) < tol
end

# Start by testing the direct bindings and just make sure that they run!
@testset "C bindings: check callable" begin
  # First, we try testing the direct C binding to the IHC code (IHCAN!)
  @test begin
      px = pt
      cf = freq
      nrep = Int32(1)
      tdres = 1.0/fs
      totalstim = Int32(dur*fs)
      cohc = 1.0
      cihc = 1.0
      species = Int32(1)
      ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
      AuditoryNerveFiber.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      true
  end

  # Next, we try testing the direct C binding to the Synapse code (Synapse!)
  @test begin
      px = pt
      ffGn = zeros(Int(ceil((length(px) + 2 * floor(7500 / (freq / 1e3))) * 1/fs * 10e3)))
      tdres = 1.0/fs
      cf = freq
      totalstim = Int32(dur*fs)
      nrep = Int32(1)
      spont = 100.0
      noiseType = 1.0
      implnt = 0.0
      sampFreq = 50000.0
      cohc = 1.0
      cihc = 1.0
      species = Int32(1)
      ihcout = zeros(length(px))
      synouttmp = zeros(length(px))

      AuditoryNerveFiber.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      AuditoryNerveFiber.Synapse!(ihcout, ffGn, tdres, cf, totalstim, nrep, spont, noiseType, implnt,
                       sampFreq, synouttmp)
      true
  end

  # Finally, we try testing the direct C binding to the SingleAN!() function
  @test begin
      px = scale_dbspl(pure_tone(freq, 0.0, dur, fs), 50.0)
      tdres = 1.0/fs
      cf = freq
      totalstim = Int32(dur*fs)
      nrep = Int32(1)
      spont = 100.0
      fibertype = 3.0
      noiseType = 1.0
      implnt = 0.0
      sampFreq = 50000.0
      cohc = 1.0
      cihc = 1.0
      species = Int32(1)
      ihcout = zeros(length(px))
      synouttmp = zeros(length(px))
      ffGn = zeros(Int(ceil((length(px) + 2 * floor(7500 / (freq / 1e3))) * 1/fs * 10e3)))
      meanrate = zeros(length(px))
      varrate = zeros(length(px))
      psth = zeros(length(px))

      AuditoryNerveFiber.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      AuditoryNerveFiber.SingleAN!(ihcout, ffGn, cf, nrep, tdres, totalstim, fibertype, noiseType,
                       implnt, meanrate, varrate, psth)
      true
  end
end

# Next we test the "wrappers" that provide a nicer interface to the bindings, beginning with simple checks that they are callable
@testset "Wrappers: check callable" begin
    # Check if the wrapper can be evaluated
    @test begin
        sim_ihc_zbc2014(zeros(Int(dur*fs)), freq)
        true
    end
    @test begin
        sim_synapse_zbc2014(zeros(Int(dur*fs)), freq)
        true
    end
    @test begin
        sim_an_zbc2014(zeros(Int(dur*fs)), freq)
        true
    end
end

# Next we check that the outputs are sensible
@testset "Wrappers: outputs" begin
    # Check that calling with n_rep > 1 produces expected results (cloning and tiling of response along primary axis)
    # Note that we need to zero-pad the inputs, otherwise responses will trail over?
    @test begin
        output_single = sim_ihc_zbc2014([pt; zeros(20000)], freq)
        output_multiple = sim_ihc_zbc2014([pt; zeros(20000)], freq; n_rep=2)
        first_repeat = all(abs.(output_single[1:10000] - output_multiple[1:10000]) .< tol)
        second_repeat = all(abs.(output_single[1:10000] - output_multiple[(1:10000) .+ 30000]) .< tol)
    end
    # Check that inner hair cell response shows expected hallmarks (response to sinusoid at CF, partial rectification)
    @test begin
        output = sim_ihc_zbc2014(pt, freq)
        # Check that
        #1 at least one element is non-zero,
        #2 minima is smaller abs val than maxima (consistent with partial rectification in IHC)
        any(output .!= 0) && (abs(minimum(output)) < abs(maximum(output)))
    end
    # Check that each output of auditory nerve response is present and non-zero
    @test begin
        outs = sim_an_zbc2014(sim_ihc_zbc2014(pt, freq), freq)
        (sum(outs[1]) != 0) && (sum(outs[2]) != 0) && (sum(outs[3]) != 0)
    end
    # Check that auditory nerve response shows expected pattern of results with changing spontaneous rate
    @test begin
        # Check analytic firing rates
        (low, medium, high) = 
            map(fiber_type -> mean(sim_an_zbc2014(sim_ihc_zbc2014(pt, freq), freq; fiber_type=fiber_type)[1]), ["low", "medium", "high"])
        analytic = (low < medium < high)
        # Check spiking firing rate (note that this approach is not great because it's stochastic and could fail randomly)
        (low, medium, high) = 
            map(fiber_type -> mean(sim_an_zbc2014(sim_ihc_zbc2014(pt, freq), freq; fiber_type=fiber_type)[3]), ["low", "medium", "high"])
        spiking = (low < medium < high)
        analytic && spiking
    end
    # Check that auditory nerve shows reasonable rate-level function (e.g., range from 10 to 30 dB is constantly increasing)
    @test begin
        inputs = [scale_dbspl(pt, level) for level in [10.0, 15.0, 20.0, 25.0, 30.0]]
        ihcouts = [sim_ihc_zbc2014(input, freq) for input in inputs]
        outputs = [sim_an_zbc2014(ihcout, freq)[1] for ihcout in ihcouts]
        means = map(Statistics.mean, outputs)
        means[5] > means[4] > means[3] > means[2] > means[1]
    end
    # Check that IHC shows reasonable frequency tuning
    @test begin
        cfs = LogRange(200.0, 5000.0, 40)
        tuning_curve = map(cf -> mean(sim_ihc_zbc2014(pt, cf)), cfs)
        abs(cfs[findmax(tuning_curve)[2]] - freq) < 500
    end
    # Check that synapse shows reasonable frequency tuning
    @test begin
        cfs = LogRange(200.0, 5000.0, 40)
        tuning_curve = map(cf -> mean(sim_synapse_zbc2014(sim_ihc_zbc2014(pt, cf), cf)), cfs)
        abs(cfs[findmax(tuning_curve)[2]] - freq) < 500
    end
    # Check that auditory nerve shows reasonable frequency tuning
    @test begin
        cfs = LogRange(200.0, 5000.0, 40)
        tuning_curve = map(cf -> mean(sim_an_zbc2014(sim_ihc_zbc2014(pt, cf), cf)[1]), cfs)
        abs(cfs[findmax(tuning_curve)[2]] - freq) < 500
    end
end

# Next we test that n_rep is handled appropriately
@testset "Wrappers: check n_rep handling" begin
    # Verify that sim_ihc_zbc2014 handles varying n_rep and returns correct lengths
    @test begin
        results = map(1:5:100) do n_rep
            return length(sim_ihc_zbc2014(pt, 1000.0; n_rep=n_rep)) == (length(pt) * n_rep)
        end
        all(results)
    end

    # Verify that sim_synapse_zbc2014 handles varying n_rep and returns correct lengths
    @test begin
        results = map(1:5:100) do n_rep
            return length(sim_synapse_zbc2014(sim_ihc_zbc2014(pt, 1000.0; n_rep=n_rep), 1000.0; n_rep=n_rep)) == (length(pt) * n_rep)
        end
        all(results)
    end

    # Verify that sim_an_zbc2014 handles varying n_rep and returns correct lengths
    @test begin
        results = map(1:5:100) do n_rep
            return length(sim_an_zbc2014(sim_ihc_zbc2014(pt, 1000.0; n_rep=n_rep), 1000.0; n_rep=n_rep)[1]) == (length(pt))
        end
        all(results)
    end

    # Verify that sim_anrate_zbc2014 handles varying n_rep and returns correct lengths
    @test begin
        results = map(1:5:100) do n_rep
            return length(sim_anrate_zbc2014(sim_ihc_zbc2014(pt, 1000.0; n_rep=n_rep), 1000.0; n_rep=n_rep)) == (length(pt) * n_rep)
        end
        all(results)
    end

    # Verify that sim_spikes_zbc2014 handles varying n_rep and returns correct lengths
    @test begin
        results = map(1:5:100) do n_rep
            return length(sim_spikes_zbc2014(sim_ihc_zbc2014(pt, 1000.0; n_rep=n_rep), 1000.0; n_rep=n_rep)) == (length(pt))
        end
        all(results)
    end
end

# Next, we will check that fractional Gaussian noise produces randomness in synaptic outputs
@testset "Behavior: fractional Gaussian noise" begin
    # Verify that if we request fGn that responses are random
	  @test begin
	      x1 = sim_anrate_zbc2014(zeros(10000), 1000.0; fractional=true)
	      x2 = sim_anrate_zbc2014(zeros(10000), 1000.0; fractional=true)
        sum(abs.(x1 .- x2)) > tol
    end
end
