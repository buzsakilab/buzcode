% Input/output functions (including batch handling) for FMAToolbox.
%
% NOTE
%
%   Although individual data files can be handled by these 'low-level' functions,
%   it is generally easier to use <a href="matlab:help SetCurrentSession">SetCurrentSession</a> and related <a href="matlab:help Data">Get...</a> functions
%   instead.
%
%
% Creating new event files
%
%   NewEvents            - Create events structure.
%   SaveEvents           - Save events to evt file.
%   SaveRippleEvents     - Save ripple events to evt file.
%
% Low-level binary input/output
%
%   LoadBinary           - Load data from a binary file.
%   SaveBinary           - Save data to binary file.
%   ResampleBinary       - Resample binary file (e.g. dat->lfp).
%   ChangeBinaryGain     - Change gain in a multiplexed binary file.
%
% Low-level data input/output
%
%   LoadEvents           - Load events from an evt file.
%   LoadParameters       - Load parameters from an xml file.
%   LoadPositions        - Load position data from a pos file.
%   LoadSpikeTimes       - Load spike times from a res file.
%   LoadSpikeFeatures    - Load spike features from file.
%   LoadSpikeWaveforms   - Load spike waveforms from a spk file.
%