% Specialized analyses for FMAToolbox.
%
% Firing map and curves, etc. (place cells, grid cells, head direction cells).
%
%   FiringCurve            - Compute spatial firing curve (e.g. for a head direction cell).
%   FiringMap              - Compute spatial firing map (e.g. for a place cell).
%   PhaseCurve             - Compute phase curve.
%   PhaseMap               - Compute phase map.
%   DefineZone             - Define a restricted zone to which analyses can subsequently be circumscribed.
%   IsInZone               - Test when the animal is in a given zone.
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
%   PhaseDistribution      - Compute phase distribution.
%   MTSpectrogram          - Compute LFP spectrogram by multi-taper estimation (chronux version).
%   MTSpectrum             - Compute LFP spectrum by multi-taper estimation.
%   SpectrogramBands       - Determine running power in physiological bands.
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
%   CCG                    - Compute multiple cross- and auto-correlograms
%   CountSpikesPerCycle    - Count number of spikes per LFP cycle.
%   CV                     - Compute coefficient of variation for a point process.
%   IsolationDistance      - Determine the isolation quality for one or more clusters.
%   MTPointSpectrum        - Compute point process spectrum by multi-taper estimation.
%   PhasePrecession        - Compute spike phase precession.
%   SelectSpikes           - Discriminate bursts vs single spikes.
%   ShortTimeCCG           - Time-varying auto/cross-correlograms of point processes.
%
% Behaviour
%
%   AngularVelocity        - Compute instantaneous angular velocity.
%   LinearVelocity         - Compute instantaneous linear velocity.
%   QuietPeriods           - Find extended periods of immobility.
%   BrainStates            - Determine brain state using LFP, EMG and quiescence.
%   RadialMaze             - Compute simple statistics for the 8-arm radial maze task.
%   RadialMazeTurns        - Compute distribution of successive turn angles in a radial maze.
%