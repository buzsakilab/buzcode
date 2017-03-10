% Specialized analyses for FMAToolbox.
%
% Firing map and curves, etc. (place cells, grid cells, head direction cells).
%
%   FiringCurve            - Compute spatial firing curve (e.g. for a head direction cell).
%   FiringMap              - Compute spatial firing map (e.g. for a place cell).
%   PhaseCurve             - Compute phase curve.
%   PhaseMap               - Compute phase map.
%   PhaseMango             - Compute phase as a function of spike rate and acceleration.
%   DefineZone             - Define a restricted zone to which analyses can subsequently be circumscribed.
%   IsInZone               - Test when the animal is in a given zone.
%   FieldShift             - Estimate firing field shifts between two conditions.
%   NormalizeFields        - Normalize one or more firing fields in space and rate.
%   ReconstructPosition    - Bayesian reconstruction of positions from spike trains.
%   TestRemapping          - Test if firing fields remap (or shift) between two conditions.
%   TestSkewness           - Test if firing field skewness changes between two conditions.
%   CompareDistributions   - Compare two N-dimensional distributions (e.g. prospective place fields).
%
% General purpose mapping.
%
%   Map                    - Map Z on (X,Y) where X, Y and Z are time-varying variables (samples).
%   MapStats               - Compute statistics for a map Z = f(X,Y), or a curve Z = f(X).
%
% Evoked potentials, peri-event time histograms, rasters.
%
%   Sync                   - Make sample timestamps relative to synchronizing events.
%   SyncHist               - Compute a histogram on event-synchronized samples (e.g. a PSTH).
%   SyncMap                - Create a map from successive event-synchronized data.
%   PETHTransition         - Find a transition point in a peri-event time histogram.
%
% Local field potentials.
%
%   FilterLFP              - Filter the local field potentials, e.g. in the theta band.
%   Phase                  - Compute instantaneous phase in LFP.
%   MTCoherence            - Compute LFP coherence by multi-taper estimation.
%   MTCoherogram           - Compute LFP coherogram by multi-taper estimation.
%   MTSpectrum             - Compute LFP spectrum by multi-taper estimation.
%   MTSpectrogram          - Compute LFP spectrogram by multi-taper estimation.
%   SpectrogramBands       - Determine running power in physiological bands.
%   CoherenceBands         - Determine running coherence in physiological bands.
%   CSD                    - Compute current source density.
%   FindRipples            - Find hippocampal ripples (100~200Hz oscillations).
%   RippleStats            - Compute descriptive stats for ripples (100~200Hz oscillations).
%   TuneArtefactTimes      - Fine-tune artefact (e.g. stimulation) times.
%   RemoveArtefacts        - Remove artefacts from local field potentials.
%   FieldPSP               - Measure field EPSPs (waveforms, slope, amplitude) across time.
%
% Spike trains.
%
%   Frequency              - Compute instantaneous frequency for a spike train.
%   CCG                    - Compute multiple cross- and auto-correlograms.
%   CCGParameters          - Reformat time series for CCG computation.
%   FitCCG                 - Fit dampened sinewave to the cross-correlogram of a pair of theta-modulated cells.
%   CountSpikesPerCycle    - Count number of spikes per LFP cycle.
%   CV                     - Compute coefficient of variation for a point process.
%   IsolationDistance      - Determine the isolation quality for one or more clusters.
%   MTPointSpectrogram     - Compute point process spectrogram by multi-taper estimation.
%   MTPointSpectrum        - Compute point process spectrum by multi-taper estimation.
%   PhasePrecession        - Compute spike phase precession.
%   SelectSpikes           - Discriminate bursts vs single spikes.
%   ShortTimeCCG           - Time-varying auto/cross-correlograms of point processes.
%   ThresholdSpikes        - Post-hoc threshold correction for spike detection.
%
% Behaviour
%
%   AngularVelocity        - Compute instantaneous angular velocity.
%   LinearVelocity         - Compute instantaneous linear velocity.
%   QuietPeriods           - Find extended periods of immobility.
%   BrainStates            - Determine brain state using LFP, EMG and quiescence.
%   Distance               - Compute instantaneous distance to a reference point.
%   RadialMaze             - Compute simple statistics for the 8-arm radial maze task.
%   RadialMazeTurns        - Compute distribution of successive turn angles in a radial maze.
%
% Data Overview
%
%   SurveyFiringMaps       - Compute and plot firing maps for all subsessions.
%   SurveyPhasePrecession  - Compute and plot phase precession plots for all subsessions.
%   SurveyTuningCurves     - Compute and plot tuning curves for all subsessions.
%   SurveyCCG              - Compute and plot cross-correlograms (CCGs) for all subsessions.
%
%
