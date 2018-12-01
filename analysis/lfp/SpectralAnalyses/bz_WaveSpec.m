function [wavespec] = bz_WaveSpec(lfp,varargin)
%[wavespec] = bz_WaveSpec(lfp) calculates the 
%wavelet transform of a signal with nfreqs frequencies in the range frange 
%[fmin fmax]. Spacing between frequencies can be 'lin' or 'log'.
%Time-frequency resolution is defined by ncyc, the number of cycles in each
%wavelet. Uses Morlet (Gabor) wavelet.
%
%
%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       'frange'    [low frequency, high frequency]     (default: [1 128])
%       'nfreqs'    number of frequencies               (default: 100
%       'roundfreqs' round freqs to unique integer vals (default: false
%                       *Note this may decrease number of freqs
%       'nfreqs'    number of frequencies               (default: 100
%       'ncyc'      number of cycles in your wavelet    (default: 5)
%       'fvector'   predefined vector of frequencies 
%       'space'     'log' or 'lin'  spacing of f's      (default: 'log')
%       'samplingRate' (only if input is not a buzcode structure)
%       'intervals'  ADD THIS - ability to spec intervals
%       'showprogress' true/false (default:false)
%       'saveMat '   put the basePath to save an LFP file
%       'MatNameExtraText'   text(X) to add to name as in: 'basename.wavespec(text).lfp.mat'
%    =========================================================================
%
%OUTPUT
%   wavespec            buzcode-style structure
%       .data           [t x nfreqs] your spectrogram
%       .timestamps     [t x 1] timestamps
%       .freqs          frequencies of each column
%       .samplingRate   (Hz)
%       .channels       Channel indices of channels filtered... taken from lfp input
%       .filterparms    a structure that holds the parameters used for
%                       filtering, for future reference
%

%TO DO:
%   -Matricise For loop?
%   -Don't need to FFT every freq...
%   -update t output for if LFP is a cell array
%
%
%Dependencies
%   WaveFilt
%   MorletWavelet
%   FConv
%
%
%Last Updated: 10/9/15
%DLevenstein
%Modified by Antonio FR, 7/18/18

%% Parse the inputs

%Parameters
parms = inputParser;
addParameter(parms,'frange',[1 128],@isnumeric);
addParameter(parms,'nfreqs',100,@isnumeric);
addParameter(parms,'ncyc',5,@isnumeric);
addParameter(parms,'space','log');
addParameter(parms,'samplingRate',[]);
addParameter(parms,'showprogress',false,@islogical);
addParameter(parms,'roundfreqs',false,@islogical);
addParameter(parms,'saveMatPath',[]);
addParameter(parms,'MatNameExtraText',[]);
addParameter(parms,'fvector',[]);

parse(parms,varargin{:})
frange = parms.Results.frange;
nfreqs = parms.Results.nfreqs;
ncyc = parms.Results.ncyc;
space = parms.Results.space;
samplingRate = parms.Results.samplingRate;
showprogress = parms.Results.showprogress;
roundfreqs = parms.Results.roundfreqs;
saveMatPath = parms.Results.saveMatPath;
MatNameExtraText = parms.Results.MatNameExtraText;
fvector = parms.Results.fvector;


%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif isempty(lfp)
    wavespec = lfp;
    return
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

si = 1./samplingRate;


%%
if ~isa(data,'single') || ~isa(data,'double')
    data = single(data);
end

%Frequencies
if ~isempty(fvector)
    freqs = fvector;
    nfreqs = length(fvector);
else
    fmin = frange(1);
    fmax = frange(2);
    if strcmp(space,'log')
        assert(fmin~=0,'Log-spaced frequencies cannot have min of 0')
        freqs = logspace(log10(fmin),log10(fmax),nfreqs);
    elseif strcmp(space,'lin')
        freqs = linspace(fmin,fmax,nfreqs);
    else
        display('Frequency spacing must be "lin" or "log".')
    end    
end

if roundfreqs
    freqs = unique(round(freqs));
end
nfreqs = size(freqs,2);

%Filter with wavelets
spec = [];
for cidx = 1:size(data,2)
    tspec = zeros(length(timestamps),nfreqs);
    for f_i = 1:nfreqs
        if showprogress
            if mod(f_i,10) == 1
                display(['freq ',num2str(f_i),' of ',num2str(nfreqs)]);
            end  
        end
        wavelet = MorletWavelet(freqs(f_i),ncyc,si);
        tspec(:,f_i) = FConv(wavelet',data(:,cidx));
    end

    if exist('celllengths','var')
        tspec = mat2cell(tspec,nfreqs,celllengths);
    end
    
    spec = cat(3,spec,tspec);
    clear tspec
end

%% Output in buzcode format
wavespec.data = spec;
wavespec.timestamps = timestamps;
wavespec.freqs = freqs;
wavespec.nfreqs = nfreqs;
wavespec.samplingRate = samplingRate;
if isstruct(lfp) && isfield(lfp,'channels')
    wavespec.channels = lfp.channels;
end
wavespec.filterparms.ncyc = ncyc;
wavespec.filterparms.nfreqs = nfreqs;
wavespec.filterparms.frange = frange;
wavespec.filterparms.space = space;

if saveMatPath
    baseName = bz_BasenameFromBasepath(saveMatPath);
    if ~isempty(MatNameExtraText)
        lfpfilename = fullfile(saveMatPath,[baseName,'.wavespec' MatNameExtraText '.lfp.mat']);    
    else
        lfpfilename = fullfile(saveMatPath,[baseName,'.wavespec.lfp.mat']);
    end
    
    s = whos('wavespec');
    if s.bytes>=1073741824%if greater than 2GB
        save(lfpfilename,'wavespec','-v7.3')
    else
        disp('wavespec variable greater than 2GB, saving as v7.3 .mat file')
        save(lfpfilename,'wavespec')
    end
end

end

