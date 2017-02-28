function [returnVar,msg] = RemoveDCfromDat_AllDat(varargin)

% USAGE:
%     RemoveDCfromDat_AllDat(fbasename)
% 
% removes DC offset in all channels from all dat files starting with
% 'fbasename'. Assumes that there is at least one xml file associated to
% one the dat file (must be created before with Neuroscope and/or NDmanger)
%     
% Adrien Peyrache, 2012

fbasename = '';

if ~isempty(varargin)
    fbasename = varargin{1};
end
d = dir([fbasename '*.dat']);
ii = 1;
fname = d(ii).name;
fname = fname(1:end-4);
while ~exist([fname '.xml'],'file')
    if ii <= length(d);
        fname = d(ii).name;
        fname = fname(1:end-4);
        ii=ii+1;
    else
        break
    end
end
if ~exist([fname '.xml'],'file')
    [fname,PathName] = uigetfile('*.xml','Find .xml containing number of channels in these .dats');
end

try 
    syst = LoadXml(fname);
catch
    error('Error - no valid .xml found')
end

for ii=1:length(d)
    [returnVar,msg] = RemoveDCfromDat(d(ii).name,syst.nChannels);
    if ~returnVar
        warning(msg)
    else
%         printf('...done\n');
    end
end