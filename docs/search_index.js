var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The best way to learn to use the package is to build from examples.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"CurrentModule = AuditoryNerveFiber","category":"page"},{"location":"examples/#Pure-tone-response-from-one-ANF","page":"Examples","title":"Pure-tone response from one ANF","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Simulating responses from a single auditory-nerve fiber is easy. The first element returned from sim_an_zbc2014 is the analytic firing rate approximation, so we can pass a pure-tone stimulus through sim_ihc_zbc2014 to sim_an_zbc2014 and then extract the first element to get the instantaneous firing rate of the auditory nerve responding to the pure tone.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AuditoryNerveFiber\nconst ANF = AuditoryNerveFiber\nusing AuditorySignalUtils\nconst ASU = AuditorySignalUtils\nusing Plots\n\n# Define variables\nfreq = 1000.0   # frequency and CF, Hz\nphase = 0.0     # starting phase, rads\ndur = 0.2       # duration, seconds\nfs = 10e4       # sampling rate, Hz\nlevel = 50.0    # level, dB SPL\n\n# Synthesize a pure tone\nx = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);\n# Simulate response \ny = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(x, freq), freq)[1];\nplot(y[1:2000], ylabel=\"Firing rate (sp/s)\", xlabel=\"Samples\")","category":"page"},{"location":"examples/#Iso-level-tuning-curve","page":"Examples","title":"Iso-level tuning curve","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Simulating an iso-level tuning curve requires a bit more work. Here, we define simple functions that synthesize a pure tone and simulate a response to that pure tone. Then, we use map to simulate an average response at several pure-tone frequencies. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AuditoryNerveFiber\nconst ANF = AuditoryNerveFiber\nusing AuditorySignalUtils\nconst ASU = AuditorySignalUtils\nusing Plots\nusing Statistics\n\n# Define variables\ncf = 1000.0     # CF, Hz\nphase = 0.0     # starting phase, rads\ndur = 0.2       # duration, seconds\nfs = 10e4       # sampling rate, Hz\nlevel = 50.0    # level, dB SPL\n\n# Define a function to synthesize a pure tone\npure_tone(freq) = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);\n# Define a function to simulate a single response\nanf_response(x) = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(x, cf), cf)[1];\n# Generate a log-spaced frequency axis\nfreqs = ASU.LogRange(200.0, 20000.0, 50)\n# Synthesize tone and simulate response at each freq\nresults = map(freq -> mean(anf_response(pure_tone(freq))), freqs)\nplot(freqs, results, ylabel=\"Firing rate (sp/s)\", xlabel=\"Frequency (Hz)\", xscale=:log)","category":"page"},{"location":"examples/#Generating-spike-trains","page":"Examples","title":"Generating spike trains","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"sim_an_zbc2014 has multiple outputs.  The first is an analytic approximation of the underlying instantaneous firing rate. The second is an analytic approximation of the variance of the underlying instantaneous variance. The third is a spike train.  Here, we can see how to extract and analyze each of these outputs.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AuditoryNerveFiber\nconst ANF = AuditoryNerveFiber\nusing AuditorySignalUtils\nconst ASU = AuditorySignalUtils\nusing Plots\n\n# Define variables\nfreq = 1000.0   # frequency and CF, Hz\nphase = 0.0     # starting phase, rads\ndur = 0.2       # duration, seconds\nfs = 10e4       # sampling rate, Hz\nlevel = 50.0    # level, dB SPL\n\n# Synthesize a pure tone\nx = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);\n# Simulate response \n(mean, var, spikes) = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(x, freq), freq);\n\n# Plot\nl = @layout [a ; b]\np1 = plot(mean, ylabel=\"Firing rate (sp/s)\", xlabel=\"Samples\")\np2 = plot(spikes, ylabel=\"Firing rate (sp/s)\", xlabel=\"Samples\")\nplot(p1, p2, layout=l)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The top row shows the analytic firing rate. The bottom row shows a single example spike train. ","category":"page"},{"location":"examples/#Extending-functions-with-multiple-dispatch","page":"Examples","title":"Extending functions with multiple dispatch","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The available interface may not always satisfy your needs. For example, you may need to simulate responses at several CFs but the default method for sim_ihc_zbc2014 is only defined for a scalar CF parameter. With Julia's multiple dispatch system, you can readily define your own extensions to the methods provided by AuditoryNerveFiber.jl. Rather than define methods for common use cases in this package, we defer to the user or packages that extend this package to define methods that suit their needs (i.e., this package is merely a thin implementation of underlying models and not a modeling toolbox). ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AuditoryNerveFiber\nconst ANF = AuditoryNerveFiber\nusing AuditorySignalUtils\nconst ASU = AuditorySignalUtils\nusing Plots\n\n# Define extension to vector-valued CF parameters\nfunction ANF.sim_ihc_zbc2014(input::Array{Float64, 1}, cf::Array{Float64, 1})\n    map(_cf -> ANF.sim_ihc_zbc2014(input, _cf), cf)\nend\n\n# Define variables\nfreq = 1000.0   # freq, Hz\ncfs = [500.0, 1000.0, 1500.0]  # CFs, Hz\nphase = 0.0     # starting phase, rads\ndur = 0.2       # duration, seconds\nfs = 10e4       # sampling rate, Hz\nlevel = 50.0    # level, dB SPL\n\n# Define a function to synthesize a pure tone\npure_tone = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);\n# Simulate IHC response at several CFs\nresults = ANF.sim_ihc_zbc2014(pure_tone, cfs; species=\"human\")\nplot([result[1:1000] for result in results], layout=3, labels=\"CF = \" .* string.(hcat(cfs...)))","category":"page"},{"location":"examples/#Plotting-neurograms","page":"Examples","title":"Plotting neurograms","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Neurograms are likewise easy to generate, so long as we organize the simulations correctly. Here, we extend sim_an_zbc2014 in a similar way as above.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AuditoryNerveFiber\nconst ANF = AuditoryNerveFiber\nusing AuditorySignalUtils\nconst ASU = AuditorySignalUtils\nusing Plots\n\n# Define extension to vector-valued CF parameters\nfunction ANF.sim_an_zbc2014(input::Array{Float64, 1}, cf::Array{Float64, 1})\n    map(_cf -> ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(input, _cf), _cf)[1], cf)\nend\n\n# Define variables\nfreq = 1000.0   # freq, Hz\ncfs = collect(ASU.LogRange(200.0, 20000.0, 100))  # CFs, Hz\nphase = 0.0     # starting phase, rads\ndur = 0.2       # duration, seconds\nfs = 10e4       # sampling rate, Hz\nlevel = 50.0    # level, dB SPL\n\n# Define a function to synthesize a pure tone\npure_tone = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);\n# Simulate IHC response at several CFs\nresults = ANF.sim_an_zbc2014(pure_tone, cfs)\n\n# Plot\nheatmap(transpose(hcat(results...))[:, 1:3000], xlabel=\"Samples\", ylabel=\"CF (#)\")","category":"page"},{"location":"examples/#Simulating-hearing-loss","page":"Examples","title":"Simulating hearing loss","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can extend our neurogram simulation to simulate what happens in the case of broad loss of outer hair cells by passing a new value to the cohc parameter in sim_ihc_zbc2014.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AuditoryNerveFiber\nconst ANF = AuditoryNerveFiber\nusing AuditorySignalUtils\nconst ASU = AuditorySignalUtils\nusing Plots\n\n# Define extension to vector-valued CF parameters\nfunction ANF.sim_an_zbc2014(input::Array{Float64, 1}, cf::Array{Float64, 1}; cohc=1.0)\n    map(_cf -> ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(input, _cf; cohc=cohc), _cf)[1], cf)\nend\n\n# Define variables\nfreq = 1000.0   # freq, Hz\ncfs = collect(ASU.LogRange(200.0, 20000.0, 100))  # CFs, Hz\nphase = 0.0     # starting phase, rads\ndur = 0.2       # duration, seconds\nfs = 10e4       # sampling rate, Hz\nlevel = 50.0    # level, dB SPL\n\n# Define a function to synthesize a pure tone\npure_tone = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);\n# Simulate IHC response at several CFs\nresults = ANF.sim_an_zbc2014(pure_tone, cfs; cohc=0.25)\n\n# Plot\nheatmap(transpose(hcat(results...))[:, 1:3000], xlabel=\"Samples\", ylabel=\"CF (#)\")","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Here, you can see the documentation for each type and method provided by AuditoryNerveFiber.jl. Higher-level functions designed for end users are documented first, while lower-level functions designed mostly for internal are documented later.","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"CurrentModule = AuditoryNerveFiber","category":"page"},{"location":"reference/#Wrappers","page":"Reference","title":"Wrappers","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Wrappers are functions that \"wrap\" around lower-level bindings and/or implementations and provide more convenience than the bindings themselves. For example, when calling wrappers a user should never have to worry about handling pointers, pre-allocating arrays, etc. These are the functions that most users will want to use in their own code. Presently, the wrappers provide access to the Zilany, Bruce, and Carney (2014) auditory-nerve model and all share the same basic signature:","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"func(input, cf; kwargs)","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"where the input is a 1D array, cf is a characteristic frequency in Hz, and kwargs are various other parameters that have default values (i.e., you can set them if you want to, but if you don't a default value will be used). Below, you can see each of the wrappers documented below, but you can also head over to Examples to see examples of how these functions can be used.","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"sim_ihc_zbc2014\nsim_synapse_zbc2014\nsim_an_zbc2014","category":"page"},{"location":"reference/#AuditoryNerveFiber.sim_ihc_zbc2014","page":"Reference","title":"AuditoryNerveFiber.sim_ihc_zbc2014","text":"sim_ihc_zbc2014(input, cf; fs=10e4, cohc=1.0, cihc=1.0, species=\"cat\")\n\nSimulates inner hair cell potential for given acoustic input.\n\nArguments\n\ninput::Array{Float64, 1}: sound pressure waveform in pascals\ncf::Float64: characteristic frequency of the IHC in Hz\nfs::Float64: sampling rate in Hz\ncohc::Float64: outer hair cell survival (from 0 to 1)\ncihc::Float64: inner hair cell survival (from 0 to 1)\nspecies::String: species, either \"cat\", \"human\" (Shera tuning), or \"human_glasberg\" (Glasberg tuning)\n\nReturns\n\noutput::Array{Float64, 1}: inner hair cell potential output\n\n\n\n\n\n","category":"function"},{"location":"reference/#AuditoryNerveFiber.sim_synapse_zbc2014","page":"Reference","title":"AuditoryNerveFiber.sim_synapse_zbc2014","text":"sim_synapse_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, fiber_type=\"high\", frac_noise=\"approximate\", noise_type=\"ffGn\", n_rep=1)\n\nSimulates synapse output for a given inner hair cell input\n\nArguments\n\ninput::Array{Float64, 1}: input hair cell potential (from simihczbc2014)\ncf::Float64: characteristic frequency of the fiber in Hz\nfs::Float64: sampling rate of the input in Hz\nfs_synapse::Float64: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.\nfiber_type::String: fiber type, one of (\"low\", \"medium\", \"high\") spontaneous rate\nfrac_noise::String: controls whether we use true or approximate fractional Gaussian noise implementation, one of (\"actual\", \"approximate\")\nnoise_type::String: whether we use ffGn or simply Gaussian noise, one of (\"ffGn\", \"Gaussian\")\nn_rep::Int64: number of repetititons to run (note that this does not appear to work correctly for the time being)\n\nReturns\n\noutput::Array{Float64, 1}: synapse output (unknown units?)\n\n\n\n\n\n","category":"function"},{"location":"reference/#AuditoryNerveFiber.sim_an_zbc2014","page":"Reference","title":"AuditoryNerveFiber.sim_an_zbc2014","text":"sim_an_zbc2014(input, cf; fs=10e4, fs_synapse=10e3, frac_noise=\"approximate\", noise_type=\"ffGn\", n_rep=1)\n\nSimulates auditory nerve output (spikes or firing rate) for a given inner hair cell input\n\nArguments\n\ninput::Array{Float64, 1}: input hair cell potential (from sim_ihc_zbc2014)\ncf::Float64: characteristic frequency of the fiber in Hz\nfs::Float64: sampling rate of the input in Hz\nfs_synapse::Float64: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.\nfiber_type::String: fiber type, one of (\"low\", \"medium\", \"high\") spontaneous rate\nfrac_noise::String: controls whether we use true or approximate fractional Gaussian noise implementation, one of (\"actual\", \"approximate\")\nnoise_type::String: whether we use ffGn or simply Gaussian noise, one of (\"ffGn\", \"Gaussian\")\nn_rep::Int64: number of repetititons to run (note that this does not appear to work correctly for the time being)\n\nReturns\n\nmeanrate::Array{Float64, 1}: analytical estimate of instantaneous firing rate\nvarrrate::Array{Float64, 1}: analytical estimate of instantaneous firing rate variance\npsth::Array{Float64, 1}: peri-stimulus time histogram (with bin width = 1/fs)\n\n\n\n\n\n","category":"function"},{"location":"reference/#Bindings","page":"Reference","title":"Bindings","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Bindings are functions that directly interface with external implementations of models in other languages.  All bindings here work generally the same way:","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Inputs are passed directly to external functions with minimal or no input checking or preprocessing\nNo defaults are provided\nNames and types of variables match those of the original-language source code","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Generally, it is expected that end-users will never need to call these functions.  Instead, they should use a wrapper or another higher-level function to access the same functionality in a a more \"user-friendly\" way.","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"IHCAN!\nSynapse!\nSingleAN!","category":"page"},{"location":"reference/#AuditoryNerveFiber.IHCAN!","page":"Reference","title":"AuditoryNerveFiber.IHCAN!","text":"IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)\n\nDirect binding to IHCAN C function in model_IHC.c\n\nPasses arguments directly to IHCAN using ccall. Arrays are converted to pointers,  functions are converted to pointers, and all other types are converted directly to  corresponding types in C. Note that while there are type checks enforced automatically by  Julia, there are no sanity checks on any arguments.\n\nArguments\n\npx::Array{Float64, 1}: sound pressure waveform in pascals\ncf::Float64: characteristic frequency of the fiber in Hz\nnrep::Int32: number of repetitions to simulate. Note that for the IHC simulation, one \"true\" simulation is conducted and then that simulation is copied and tiled (because there is no randomness in the IHC simulation) to simulate multiple times.\ntdres::Float64: time-domain resolution (i.e., reciprocal of sampling rate)\ntotalstim::Int32: number of samples in simulation\ncohc::Float64: outer hair cell survival (from 0 to 1)\ncihc::Float64: inner hair cell survival (from 0 to 1)\nspecies::Int32: species, either (1 = cat, 2 = humans with Shera tuning, 3 = humans with Glasberg tuning)\nihcout::Array{Float64, 1}: array of same size as px, used to store output from C\n\n\n\n\n\n","category":"function"},{"location":"reference/#AuditoryNerveFiber.Synapse!","page":"Reference","title":"AuditoryNerveFiber.Synapse!","text":"Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, \n         synouttmp)\n\nDirect binding to Synapse C function in model_Synapse.c\n\nPasses arguments directly to Synapse using ccall. Arrays are converted to pointers,  functions are converted to pointers, and all other types are converted directly to  corresponding types in C. Note that while there are type checks enforced automatically by  Julia, there are no sanity checks on any arguments.\n\nArguments\n\nihcout::Array{Float64, 1}: output from IHC simulation (IHCAN!)\ntdres::Float64: time-domain resolution (i.e., reciprocal of sampling rate)\ncf::Float64: characteristic frequency of the fiber in Hz\ntotalstim::Int32: number of samples in simulation\nnrep::Int32: number of repetitions to simulate.\nspont::Float64: spontaneous rate, either (0.1 == low spont fiber, 4.0 == medium spont fiber, 100.0 == high spont fiber)\nnoiseType::Float64: whether we use ffGn or Gaussian noise (1.0 == ffGn, 0.0 == Gaussian)\nimplnt::Float64: whether or not to use exact implementation of fractional Gaussian noise, either (1.0 == actual, 0.0 == approximate)\nsampFreq::Float64: sampling frequency of the power law stage in Hz. Simulations are decimated to sampFreq from 1/tdres before the power law stage and then upsampled back to the original sampling rate. The product of tdres and sampFreq, which indicates the amount to decimate by, must be an integer\nsynouttmp::Array{Float64, 1}: array of same size as ihcout, used to store output from C\n\n\n\n\n\n","category":"function"},{"location":"reference/#AuditoryNerveFiber.SingleAN!","page":"Reference","title":"AuditoryNerveFiber.SingleAN!","text":"SingleAN!(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, \n          varrate, psth)\n\nDirect binding to Synapse C function in model_Synapse.c\n\nPasses arguments directly to Synapse using ccall. Arrays are converted to pointers,  functions are converted to pointers, and all other types are converted directly to  corresponding types in C. Note that while there are type checks enforced automatically by  Julia, there are no sanity checks on any arguments.\n\nArguments\n\nihcout::Array{Float64, 1}: output from IHC simulation (IHCAN!)\ncf::Float64: characteristic frequency of the fiber in Hz\nnrep::Int32: number of repetitions to simulate.\ntdres::Float64: time-domain resolution (i.e., reciprocal of sampling rate)\ntotalstim::Int32: number of samples in simulation\nfibertype::Float64: fiber type, either (1.0 == low, 2.0 == med, 3.0 == high)\nnoiseType::Float64: whether we use ffGn or Gaussian noise (1.0 == ffGn, 0.0 == Gaussian)\nimplnt::Float64: whether or not to use exact implementation of fractional Gaussian noise, either (1.0 == use, 0.0 == approximate)\nmeanrate::Array{Float64, 1}: array of same size as ihcout, used to store analytical firing rate output\nvarrate::Array{Float64, 1}: array of same size as ihcout, used to store analytical firing rate variance output\npsth::Array{Float64, 1}: array of same size as ihcout, used to store empirical PSTH output\n\n\n\n\n\n","category":"function"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"AuditoryNerveFiber.jl is a Julia package that provides access to implementations of auditory-nerve models.  This is the documentation for AuditoryNerveFiber.jl, and here you can learn about how to use the package and access information about details of the underlying implementations.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = AuditoryNerveFiber","category":"page"}]
}
