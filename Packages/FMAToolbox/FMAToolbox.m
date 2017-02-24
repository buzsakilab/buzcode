% Freely Moving Animal (FMA) Toolbox
%
%  INTRODUCTION
%
%    The purpose of the Freely Moving Animal (FMA) Toolbox is to help analyze
%    electrophysiological and behavioral data recorded from freely moving animals.
%
%    These data typically include broadband brain signals (e.g. local field potentials),
%    spike data (timestamps as well as full waveforms), moment-to-moment position of
%    the animal, and behavioral events (stimulus presentation, photodetector crossing,
%    brain stimulation, etc.)
%
%    <a href="http://fmatoolbox.sourceforge.net/">FMAToolbox</a> is part of a larger data analysis framework, including <a href="http://klusters.sourceforge.net/">Klusters</a> (a powerful
%    and easy-to-use cluster cutting application), <a href="http://neuroscope.sourceforge.net/">NeuroScope</a> (an advanced viewer for
%    electrophysiological and behavioral data), and <a href="http://ndmanager.sourceforge.net/">NDManager</a> (a simple application to
%    manage recording parameters and data preprocessing).
%
%  INSTALLATION
%
%    The easiest way to install FMAToolbox is using the <a href="matlab:run pathtool">pathtool</a>. Click 'Add with subfolders...'
%    and select the 'FMAToolbox' directory.
%
%    For efficiency reasons, the toolbox includes a few functions in C/C++. These can be
%    easily compiled for your platform. In the command window, type:
%
%      >> compilefma;
%
%    or simply click <a href="matlab:run compilefma">here</a>.
%
%    Some functions in the toolbox depend on other Matlab toolboxes:
%
%      * <a href="http://www.artefact.tk/software/matlab/xml/">xmltree</a> to read XML parameter files
%      * <a href="http://www.chronux.org">chronux</a> to perform time-frequency analyses
%      * <a href="http://sourceforge.net/projects/mym/">mYm</a> to provide database functionality
%      * <a href="http://sites.google.com/site/oliverwoodford/software/export_fig">export_fig</a> (optional) to export figures as PNG, EPS, PDF, etc.
%
%  LICENSE
%
%    This toolbox was developed by <a href="mailto:michael.zugaro@college-de-france.fr?subject=FMAToolbox" title="FMAToolbox">Michaël Zugaro</a> at the <a href="http://www.lppa.college-de-france.fr/EN/equipes/people/Zugaro/index.htm">LPPA</a> (CRNS-Collège de France, Paris,
%    France). It is free software distributed under the <a href="http://www.gnu.org/licenses/gpl.html">General Public License (GPL)</a>.
%
%  CONTENTS
%
%    The toolbox is organized in distinct functional categories:
%
%      Data       - Experimental data handling
%      General    - General-purpose data processing and statistical tests.
%      Analyses   - Specialized analyses (PSTH, place maps, tuning curves...)
%      Plot       - Enhanced and specialized plotting functions.
%      Database   - Database access and storage
%      IO         - Input/output helper functions, including batch handling
%
%  OVERVIEW
%
%    Functions provided in the <a href="matlab:help Data">Data</a> category can be used to read experimental data files.
%    The file format is the same as in other programs such as <a href="http://klusters.sourceforge.net">Klusters</a>, <a href="http://neuroscope.sourceforge.net">NeuroScope</a> and
%    <a href="http://ndmanager.sourceforge.net">NDManager</a> (see e.g. <a href="matlab:help http://neuroscope.sourceforge.net/UserManual/data-files.html">http://neuroscope.sourceforge.net/UserManual/data-files.html</a>).
%    However, the rest of the toolbox is independent from these file formats, so you
%    are free to use any format you like (of course, you would then have to add your own
%    code to read your files).
%
%    Upon loading, the data are stored in simple matrices where the first column contains
%    timestamps and the remaining columns contain the corresponding values. For instance,
%    moment-to-moment positions of the animal are stored as
%
%              0        100.625         74.375
%         0.0256         99.375             75
%         0.0512         99.375             75
%         0.0768        100.625             75
%         0.1024        100.625             75
%            ...            ...            ...
%
%    where the first column is time (in seconds), the second is the X coordinate of the
%    animal, and the third is its Y coordinate. Positions stored in such a matrix are
%    referred to as position 'samples'.
%
%    Many functions in FMAToolbox provide basic processing of 'samples', such as time
%    selection, interpolation or smoothing. These are grouped in the <a href="matlab:help General">General</a> category,
%    together with statistical measures and tests, interval handling, binning, etc.
%
%    More specialized analysis functions are grouped in the <a href="matlab:help Analyses">Analyses</a> category, including
%    peri-event time histograms (PETH), place maps, tuning curves, spectrograms, current
%    source density, ripple detection, etc.
%
%    Results can be plotted using enhanced or specialized functions provided in the <a href="matlab:help Plot">Plot</a>
%    category, and stored in a local or network database using functions in the <a href="matlab:help Database">Database</a>
%    category.
%
%    Last, the <a href="matlab:help IO">IO</a> category contains low-level disk access functions (which should normally
%    not be used directly), as well as batch-processing handling functions.
%
%  GETTING STARTED
%
%    A typical data analysis session would start by loading the experimental data:
%
%      >> SetCurrentSession
%
%    A dialog pops up, where you can navigate your disk to choose the recording session to
%    analyze. Files are then loaded and kept in memory. You can now access the data very
%    easily. For instance, if you need the position samples, do:
%
%      >> p = GetPositions;
%
%    To plot the firing map of the place cell corresponding to cluster 5 on tetrode 3:
%
%      >> s = GetSpikes([3 5]);
%      >> map = FiringMap(p,s);
%      >> figure;
%      >> PlotColorMap(map.rate,map.time); % time is used to dimm the colors
%
%    To plot theta phase precession information (phase as a function of position, phase
%    distribution as a function of position, and phase as a function of firing rate) for
%    the same neuron:
%
%      >> lfp = GetLFP(20); % theta is measured on channel 20
%      >> theta = FilterLFP(lfp,'passband','theta');
%      >> phi = Phase(theta);
%      >> precession = PhasePrecession(p,s,phi);
%      >> figure;
%      >> PlotPhasePrecession(precession);
%
%    Now let us plot a spike raster synchronized on stimulation events:
%
%      >> stims = GetEvents('stimulation');
%      >> [sync,i] = Sync(s,stims);
%      >> figure;
%      >> PlotSync(sync,i);
%
%    and the corresponding PSTH and time-varying spike distribution:
%
%      >> [h,t] = SyncHist(sync,i);
%      >> figure;
%      >> bar(t,h); % PSTH
%      >> [h,t,n] = SyncHist(sync,i,'mode','dist');
%      >> figure;
%      >> PlotColorMap(h,1,'x',t,'y',n); % distribution
%
%    For more information, read the help topics listed above (CONTENTS).
%
%  NOTE
%
%    You may wish to add the following line to your <a href="matlab:edit startup.m">startup</a> file:
%
%      Browse('on');
%
%    This will automatically activate <a href="matlab:help Browse">browsing features</a> in every figure.

help FMAToolbox;