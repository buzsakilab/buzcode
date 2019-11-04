function [ comod ] = bz_Comodulogram(lfp,specparms,figparms)
%[ comod ] = bz_Comodulogram(lfp,specparms,figparms) calculates the
%comodulogram (i.e. power-power correlation) for an lfp file.
%
%
%INPUT
%   LFP         LFP structure in buzcode format:
%               a structure with fields lfp.data, lfp.timestamps, lfp.refCh
%               (optional). If only one lfp channel is provided, that is unsed
%               for corr computation. If two, pairwise correlation is
%               performed. If one refCh and n data channels are provided
%               the correlation of the refCh against all data channels is
%               performed. 

%   specparms   structure of parameters for the spectrogram
%       .frange     [min max] frequency
%       .nfreqs     number of frequencies 
%       .space      spacing between freqs, 'lin' or 'log'
%       .fvector    predefined vector of frequencies 
%       .specnorm   normalization for spectral power,
%                   options: 'mean','logmean','log' (default: log)
%       .numvarbins number of bins for your external variable
%       .varnorm    normalization for the variable,
%                   options: 'percentile', 'none' (default: 'none')
%       .type       'wavelet' or 'FFT'      (default: 'wavelet')
%      -if type: 'wavelet'-
%       .ncyc       number of wavelet cycles (recommended: ~5)
%      -if type: 'FFT'-
%       .winsize (s)
%       .overlap
%   figparms    (optional) parameters for the figure
%       .plotname   nametag for the saved figure
%       .figfolder  folder to save the figure in
%       .baseName   baseName of the recording for the figure
%
% DLevenstein 2017
% Modified by Antonio FR, 7/18/18

%TO DO 

%% Parse the inputs

%Parameters
parms = inputParser;
addParameter(parms,'frange',[1 128],@isnumeric);
addParameter(parms,'nfreqs',100,@isnumeric);
addParameter(parms,'ncyc',5,@isnumeric);
addParameter(parms,'space','log');
addParameter(parms,'samplingRate',[]);
addParameter(parms,'showprogress',false,@islogical);
addParameter(parms,'saveMat',false);
addParameter(parms,'fvector',[]);
addParameter(parms,'specnorm','log');
addParameter(parms,'type','wavelet');
addParameter(parms,'winsize',1,@isnumeric);
addParameter(parms,'overlap',0.5,@isnumeric);


parse(parms,specparms)
specparms.frange = parms.Results.frange;
specparms.nfreqs = parms.Results.nfreqs;
specparms.ncyc = parms.Results.ncyc;
specparms.space = parms.Results.space;
specparms.samplingRate = parms.Results.samplingRate;
specparms.showprogress = parms.Results.showprogress;
specparms.saveMat = parms.Results.saveMat;
specparms.fvector = parms.Results.fvector;
specparms.specnorm = parms.Results.specnorm;
specparms.type = parms.Results.type;
specparms.winsize = parms.Results.winsize;
specparms.overlap = parms.Results.overlap;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
    if isfield(lfp,'refCh')
       refCh = lfp.refCh;
    else
        refCh = [];
    end
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


%% Calculate the spectrogram - FFT or WVLT

switch specparms.type
    
    case 'wavelet'
        %Calcualte the Wavelet Transform
        if size(lfp.data,2) == 1
            [wavespec] = bz_WaveSpec(single(data),...
                'frange',specparms.frange,'nfreqs',specparms.nfreqs,'ncyc',specparms.ncyc,...
                'samplingRate',lfp.samplingRate,'space',specparms.space,'fvector',specparms.fvector);
            spec = wavespec.data';
        elseif size(lfp.data,2) > 1 
            for i = 1:size(lfp.data,2)
                [wavespec] = bz_WaveSpec(single(data(:,i)),...
                    'frange',specparms.frange,'nfreqs',specparms.nfreqs,'ncyc',specparms.ncyc,...
                    'samplingRate',lfp.samplingRate,'space',specparms.space,'fvector',specparms.fvector);
                spec{i} = wavespec.data';                
            end
        end
        if ~isempty(refCh)
            [wavespec] = bz_WaveSpec(single(refCh),...
                'frange',specparms.frange,'nfreqs',specparms.nfreqs,'ncyc',specparms.ncyc,...
                'samplingRate',lfp.samplingRate,'space',specparms.space,'fvector',specparms.fvector);
            specRef = wavespec.data';
        else
            specRef = [];
        end
            spectimestamps = timestamps; %Wavelet timestamp are same as LFP        
            comod.freqs = wavespec.freqs;
            
	case 'FFT'
        %Calculate the frequences to use
        if ~isempty(specparms.fvector)
            comod.freqs = fvector;
        else
            switch specparms.space
                case 'log'
                    comod.freqs = logspace(log10(specparms.frange(1)),...
                        log10(specparms.frange(2)),specparms.nfreqs);
                case 'lin'
                    comod.freqs = linspace(specparms.frange(1),...
                        specparms.frange(2),specparms.nfreqs);  
            end
        end
        %Calculate the FFT spectrogram parameters - covert from s to sf
        winsize = specparms.winsize*samplingRate;
        noverlap = specparms.noverlap*samplingRate;
        %Calculate the FFT spectrogram
        if size(lfp.data,2) == 1
            [spec,~,spectimestamps] = spectrogram(single(data),...
                winsize,noverlap,comod.freqs,samplingRate);
        elseif size(lfp.data,2) == 2
            for i = 1:2
                [specT,~,spectimestamps] = spectrogram(single(data(:,i)),...
                    winsize,noverlap,comod.freqs,samplingRate);
                spec{i} = specT;
            end
        end
        if ~isempty(refCh)
            [specRef,~,spectimestamps] = spectrogram(single(refCh),...
                winsize,noverlap,comod.freqs,samplingRate);            
        else
            specRef = [];
        end
            spectimestamps = spectimestamps'+timestamps(1); %account for any time offset
end


%% Calculate the power-power correlations
if isempty(refCh)
    if size(lfp.data,2) == 1
       spec = log10(abs(spec)); %Log-transform power
       comod.corrs = corr(spec','type','spearman');
    elseif size(lfp.data,2) == 2
       for i = 1:2
           spec{i} = log10(abs(spec{i})); %Log-transform power
       end
       comod.corrs = corr(spec{1}',spec{2}','type','spearman');   
    end
elseif exist('refCh')
    specRef = log10(abs(specRef)); %Log-transform power
    if size(lfp.data,2) == 1
       spec = log10(abs(spec)); %Log-transform power
       comod.corrs = corr(specRef',spec','type','spearman');
    elseif size(lfp.data,2) > 1
       for i = 1:size(lfp.data,2)
           spec{i} = log10(abs(spec{i})); %Log-transform power
           comod.corrs{i} = corr(specRef',spec{i}','type','spearman');     
       end
    end    
end


%% Plot
% needs fix for multiple channels

if exist('figparms','var') && isempty(refCh)  %This whole figure thing can be better.
corrcolor= [makeColorMap([1 1 1],[0 0 0.8],[0 0 0]);...
    makeColorMap([0 0 0],[0.8 0 0],[1 1 1])];
figure
colormap(corrcolor)
    imagesc(log2(comod.freqs),log2(comod.freqs),comod.corrs)
    colorbar
    ColorbarWithAxis([-0.4 0.4],'Power-Power Correlation (rho)')
    LogScale('xy',2)
    xlabel('f (Hz)');ylabel('f (Hz)')
    
NiceSave(['Comodulogram',figparms.plotname],figparms.figfolder,figparms.baseName)

end

