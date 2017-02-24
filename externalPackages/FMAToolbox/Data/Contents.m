% Session data handling functions for FMAToolbox.
%
% FMAToolbox provides a number of functions to handle data files stored
% in the open and hardware-independent <a href="http://neuroscope.sourceforge.net/UserManual/data-files.html">formats</a> used by <a href="http://klusters.sourceforge.net">Klusters</a>, <a href="http://neuroscope.sourceforge.net">NeuroScope</a>
% and <a href="http://ndmanager.sourceforge.net">NDManager</a>. To read data stored in other formats, one can either use
% format converter tools (such as <a href="http://ndmanager.sourceforge.net">NDManager</a> plugins), use dedicated Matlab
% toolboxes (provided by hardware vendors), or develop custom Matlab code.
%
% Data files for a typical recording session include several spike timestamp
% (.res.n) and cluster (.clu.n) files (one per tetrode), a position file (.pos),
% one or more event files (.evt), and a local field potential file (.lfp).
% Additional information can be provided by spike waveform files (.spk.n) and
% wideband data files (.dat).
%
% Although any single file can be read individually using <a href="matlab:help IO">low-level</a> functions,
% the easiest way to handle all data for a given recording session is to call
% <a href="matlab:help SetCurrentSession">SetCurrentSession</a> (which automatically reads the session parameters and loads
% all relevant data) and then use the Get... functions listed below.
%
% NOTE
%
%   In order to avoid repeatedly reading data from disk, spikes, positions and
%   events are loaded once and stored in the global structure DATA. However,
%   LFP data is usually too large to keep in memory and must be read from disk
%   on-demand.
%
% Loading data
%
%   SetCurrentSession        - Load all data for a given recording session.
%   GetCurrentSession        - Get information about current session.
%
%
% Accessing data
%
%   GetEvents                - Get events.
%   GetLFP                   - Get local field potentials.
%   GetPositions             - Get position samples.
%   GetSpikes                - Get spike timestamps.
%   GetSpikeTimes            - Get spike timestamps.
%   GetSpikeFeatures         - Get spike features.
%   GetSpikeWaveforms        - Get spike waveforms.
%   GetUnits                 - Get list of units.
%   GetWidebandData          - Get wideband data.
%
% Data overview
%
%   SurveyFiringMaps         - Plot firing maps for all conditions in the same figure.
%   SurveyPhasePrecession    - Plot phase precession plots for all conditions.
%
% Default values
%
%   CustomDefaults           - User-defined custom default values for function properties.
%   GetCustomDefaults        - Get custom default value for a given function property.
%   Settings                 - Global default values.
%