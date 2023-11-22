
function [ica] = bz_RunIca(varargin)
% [ica] = bz_RunIca(varargin)
%
% Perform Independent Component Analysis (ICA) decomposition of input data 
% using the logistic infomax ICA algorithm of Bell & Sejnowski (1995) with 
% the natural gradient feature of Amari, Cichocki & Yang. 
%
% INPUTS
% <optional>
% basepath      Default pwd 
% lfp           a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels.
%               If not provided, runs bz_getLFP('all') on basepath 
% passband      Prefiltering passband interval, default [30 300]
% saveMat       Save results, default true.
% force         Force analysis (disable loading option if already computed, 
%                   default false)
%
% TO DO: INCLUDE IMPORTANTS IMPUTS TO runica.mat as additional arguments!
% 
% OUTPUT
% ica           a buzcode structure with the following fields:
% .data         independent components organized by explanatory variance.
% .timestamps
% .sphere       data sphering matrix (chans,chans)
% .weights      ICA weight matrix (comps,chans)
% .meanvar      Explained variance.
%
% Manu Valero 2020
% This functions is a buzcode wrapper for the function 'runica.mat' from Scott
% Makeig (CNL/The Salk Institute, La Jolla, 1996-)
% Copyright (C) 2004-2011 by 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'passband',[30 300],@isnumeric)
addParameter(p,'nICs',8,@isnumeric)

parse(p,varargin{:});
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
force = p.Results.force;
passband = p.Results.passband;
nICs = p.Results.nICs;

% Deal with inputs
prevBasepath = pwd;
cd(basepath);

targetFile = dir('*.ica.channelInfo.mat');
if ~isempty(targetFile) && ~force
    disp('ICA already computed! Loading file.');
    load(staFile.name);
    return
end

if isempty(lfp)
    try lfp = bz_GetLFP('all');
    catch
        error('LFP not found!');
    end
end
% Filtering (Schomburg et al, 2014)
disp('Filtering...');
lfpFilt = bz_Filter(lfp,'passband',passband,'order',4,'filter','butter');


% Run runica
[weights,sphere,meanvar,bias,signs,lrates,data,y] = runica(double(lfpFilt.data'),'pca',nICs);

ica.data = data';
ica.timestamps = lfp.timestamps;
ica.sphere = sphere;
ica.weights = weights;
ica.meanvar = meanvar;
ica.samplingRate = lfp.samplingRate;
ica.channels = lfp.channels;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.ica.channelInfo.mat'],'ica','-v7.3');
end

cd(prevBasepath);
end