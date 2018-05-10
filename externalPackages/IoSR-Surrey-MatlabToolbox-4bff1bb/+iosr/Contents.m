% +IOSR
% 
%   Contents file for +IOSR and its subfolders.
%   
%   +IOSR
%   iosr.install                         - Set search paths, and download and install dependencies
%   
%   +IOSR/+ACOUSTICS
%   iosr.acoustics.irStats               - Calculate RT, DRR, Cte, and EDT for impulse response file
%   iosr.acoustics.rtEst                 - Estimate reverberation time based on room size and absorption
%   
%   +IOSR/+AUDITORY
%   iosr.auditory.azimuth2itd            - Convert azimuth in degrees to ITD
%   iosr.auditory.binSearch              - Conduct a binary search
%   iosr.auditory.calcIld                - Calculate normalised interaural level difference
%   iosr.auditory.chXcorr                - Calculate cross-correlograms with a wide range of options
%   iosr.auditory.chXcorr2               - Calculate cross-correlograms with a range of options
%   chXcorr2_c.c
%   chXcorr_c.c
%   iosr.auditory.createWindow           - Create a Hann or exp. window with specified onsets/offsets
%   iosr.auditory.dupWeight              - Calculate duplex weighting coefficients for ITD and ILD
%   iosr.auditory.erbRate2hz             - Convert ERB rate to Hz
%   iosr.auditory.freqMulti              - Calculate frequency coefficient for ITD-azimuth warping
%   iosr.auditory.gammatoneFast          - Produce an array of responses from gammatone filters via FFT
%   iosr.auditory.hz2erbRate             - Convert Hz to ERB rate
%   iosr.auditory.iso226                 - ISO 226:2003 Normal equal-loudness-level contours
%   iosr.auditory.itd2azimuth            - Convert ITD to azimuth
%   iosr.auditory.lindemannInh           - Signal pre-processing for Lindemann's cross-correlation
%   iosr.auditory.loudWeight             - Calculate loudness weighting coefficients
%   iosr.auditory.makeErbCFs             - Make a series of center frequencies equally spaced in ERB-rate
%   iosr.auditory.meddisHairCell         - Calculate Ray Meddis' hair cell model for a number of channels
%   iosr.auditory.perceptualCentroid     - Perceptual spectral centroid
%   iosr.auditory.perceptualCentroid2    - Alternative perceptual spectral centroid
%   iosr.auditory.xcorrLindemann         - Cross-correlation based on Lindemann's precedence model
%   xcorrLindemann_c.c
%   
%   +IOSR/+BSS
%   iosr.bss.applyIdealMasks             - Calculate and apply ideal masks via STFT
%   iosr.bss.applyMask                   - Apply a time-frequency mask to an STFT
%   iosr.bss.calcImr                     - Calculates the Ideal Mask Ratio (IMR)
%   iosr.bss.calcSnr                     - Calculate the separation SNR
%   iosr.bss.cfs2fcs                     - Calculate gammatone crossover frequencies
%   iosr.bss.example                     - Determine STFT parameters
%   iosr.bss.generateMixtures            - Generate arrays of mixtures from targets and interferers
%   iosr.bss.getFullMask                 - Convert frame rate mask to a sample-by-sample mask
%   iosr.bss.idealMasks                  - Calculate ideal time-frequency masks from STFTs
%   iosr.bss.mixture                     - Class of sound source separation mixture
%   iosr.bss.resynthesise                - Resynthesise a target from a time-frequency mask
%   iosr.bss.source                      - Class of sound source separation source
%   
%   +IOSR/+DSP
%   iosr.dsp.audio                       - Abstract superclass providing audio-related properties and methods
%   iosr.dsp.autocorr                    - Perform autocorrelation via FFT
%   iosr.dsp.convFft                     - Convolve two vectors using FFT multiplication
%   iosr.dsp.istft                       - Calculate the Inverse Short-Time Fourier Transform
%   iosr.dsp.lapwin                      - Laplace window
%   iosr.dsp.localpeaks                  - Find local peaks and troughs in a vector
%   iosr.dsp.ltas                        - Calculate the long-term average spectrum of a signal
%   iosr.dsp.matchEQ                     - Match the LTAS of a signal to an arbitrary spectral magnitude
%   iosr.dsp.rcoswin                     - Raised cosine window
%   iosr.dsp.rms                         - Calculate the rms of a vector or matrix
%   iosr.dsp.sincFilter                  - Apply a near-ideal low-pass or band-pass brickwall filter
%   iosr.dsp.smoothSpectrum              - Apply 1/N-octave smoothing to a frequency spectrum
%   iosr.dsp.stft                        - Calculate the short-time Fourier transform of a signal
%   iosr.dsp.vsmooth                     - Smooth a vector using mathematical functions
%   
%   +IOSR/+FIGURES
%   iosr.figures.chMap                   - Create a monochrome-compatible colour map
%   iosr.figures.cmrMap                  - Create a monochrome-compatible colour map
%   iosr.figures.multiwaveplot           - Stacked line plots from a matrix or vectors
%   iosr.figures.subfigrid               - Create axis positions for subfigures
%   
%   +IOSR/+GENERAL
%   iosr.general.cell2csv                - Output a cell array to a CSV file
%   iosr.general.checkMexCompiled        - Check if mex file is compiled for system
%   iosr.general.getContents             - Get the contents of a specified directory
%   iosr.general.updateContents          - Create a Contents.m file including subdirectories
%   iosr.general.urn                     - Generate random number sequence without duplicates
%   
%   +IOSR/+STATISTICS
%   iosr.statistics.boxPlot              - Draw a box plot
%   iosr.statistics.functionalBoxPlot    - Draw a functional boxplot
%   iosr.statistics.functionalPlot       - Abstract superclass for functional plots
%   iosr.statistics.functionalSpreadPlot - Draw a functional plot showing data spread
%   iosr.statistics.getRmse              - Calculate the root-mean-square error between input data
%   iosr.statistics.laprnd               - Pseudorandom numbers drawn from the Laplace distribution
%   iosr.statistics.qqPlot               - Quantile-quantile plot with patch option
%   iosr.statistics.quantile             - Quantiles of a sample via various methods
%   iosr.statistics.statsPlot            - An abstract superclass for classes that plot statistics
%   iosr.statistics.tab2box              - Prepare tabular data for boxPlot function
%   iosr.statistics.trirnd               - Pseudorandom numbers drawn from the triangular distribution
%   
%   +IOSR/+SVN
%   iosr.svn.buildSvnProfile             - Read data from files tagged with SVN keywords
%   iosr.svn.headRev                     - Retrieve the head revision for specified files
%   iosr.svn.readSvnKeyword              - Read data from a file tagged with an SVN keyword
%    
%   This file was generated by updateContents.m on 31 May 2017 at 17:06:32.
