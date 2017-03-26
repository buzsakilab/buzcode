function MakeXMLFromProbeMaps(basepath,varargin)
%MakeXMLFromProbeMaps - Generate a .xml file to accompany .dat files for a
%recording in the neuroscope/klusters/ndmanager system.  Uses a library of
%probe map layouts.
%
%  USAGE
%
%    MakeXMLFromProbeMaps(basepath,ProbeFileName1,ProbeFileName2...)
%       Example:
%    MakeXMLFromProbeMaps(cd,'NRX_Buzsaki64_8X8','NRX_Buzsaki64_6X10');
%
%    Writes a standardized .xml file based on a user-selection of probe
%    maps and in a sequence specified by the user (ie 64site probe first
%    then 32site probe second).  Probe maps can be found at:
%    /buzcode/tree/master/generalComputation/geometries
%
%  INPUT
%
%    basepath       path to directory to which to write output xml file and
%                   where to potentially find .rhd file 
%    ProbeFileNames names of .xlsx files specifying probe geometries, must 
%                   be on the path, ie from
%                   buzcode/GeneralComputation/Geometries).
%                   Must have the following format features:
%                       - A column with row 1 having text "BY VERTICAL POSITION/SHANK (IE FOR DISPLAY)".
%                       In this row must be cells specifying SHANK 1, 
%                       SHANK 2, etc.  
%                       - Two rows right of that a column with channel
%                       numbers listed from superficial to deep in the
%                       shank group specified two to the left.
%
%  OUTPUT
%
%    (.xml file written to disk at basepath)
%
%  SEE
%
%    See also 
% Copyright (C) 2017 Brendon Watson
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Defaults
numchans = 1;
samprate = 20000;
bitspersample = 16;
voltageRange = 20;
amplification = 1000;
lfpsamprate = 1250;
pointsperwaveform = 32;
peakpointinwaveform = 16;
numfeatures = 4;

%% Define text components to assemble later
chunk1 = {'<?xml version=''1.0''?>';...
'<parameters version="1.0" creator="neuroscope-2.0.0">';...
' <acquisitionSystem>';...
['  <nBits>' num2str(bitspersample) '</nBits>']};

channelcountlinestart = '  <nChannels>';
channelcountlineend = '</nChannels>';

chunk2 = {['  <samplingRate>' num2str(samprate) '</samplingRate>'];...
['  <voltageRange>' num2str(voltageRange) '</voltageRange>'];...
['  <amplification>' num2str(amplification) '</amplification>'];...
'  <offset>0</offset>';...
' </acquisitionSystem>';...
' <fieldPotentials>';...
['  <lfpSamplingRate>' num2str(lfpsamprate) '</lfpSamplingRate>'];...
' </fieldPotentials>';...
' <files>';...
'  <file>';...
'   <extension>lfp</extension>';...
['   <samplingRate>' num2str(lfpsamprate) '</samplingRate>'];...
'  </file>';...
'  <file>';...
'   <extension>whl</extension>';...
'   <samplingRate>39.0625</samplingRate>';...
'  </file>';...
' </files>';...
' <anatomicalDescription>';...
'  <channelGroups>'};

anatomygroupstart = '   <group>';%repeats w every new anatomical group
anatomychannelnumberline_start = ['    <channel skip="0">'];%for each channel in an anatomical group - first part of entry
anatomychannelnumberline_end = ['</channel>'];%for each channel in an anatomical group - last part of entry
anatomygroupend = '   </group>';%comes at end of each anatomical group

chunk3 = {' </channelGroups>';...
  '</anatomicalDescription>';...
 '<spikeDetection>';...
  ' <channelGroups>'};%comes after anatomical groups and before spike groups

spikegroupstart = {'  <group>';...
        '   <channels>'};%repeats w every new spike group
spikechannelnumberline_start = ['    <channel>'];%for each channel in a spike group - first part of entry
spikechannelnumberline_end = ['</channel>'];%for each channel in a spike group - last part of entry
spikegroupend = {'   </channels>';...
   ['    <nSamples>' num2str(pointsperwaveform) '</nSamples>'];...
   ['    <peakSampleIndex>' num2str(peakpointinwaveform) '</peakSampleIndex>'];...
   ['    <nFeatures>' num2str(numfeatures) '</nFeatures>'];...
    '  </group>'};%comes at end of each spike group

chunk4 = {' </channelGroups>';...
 '</spikeDetection>';...
 '<neuroscope version="2.0.0">';...
  '<miscellaneous>';...
   '<screenGain>0.2</screenGain>';...
   '<traceBackgroundImage></traceBackgroundImage>';...
  '</miscellaneous>';...
  '<video>';...
   '<rotate>0</rotate>';...
   '<flip>0</flip>';...
   '<videoImage></videoImage>';...
   '<positionsBackground>0</positionsBackground>';...
  '</video>';...
  '<spikes>';...
  '</spikes>';...
  '<channels>'};

channelcolorstart = ' <channelColors>';...
channelcolorlinestart = '  <channel>';
channelcolorlineend = '</channel>';
channelcolorend = {'  <color>#0080ff</color>';...
    '  <anatomyColor>#0080ff</anatomyColor>';...
    '  <spikeColor>#0080ff</spikeColor>';...
   ' </channelColors>'};

channeloffsetstart = ' <channelOffset>';
channeloffsetlinestart = '  <channel>';
channeloffsetlineend = '</channel>';
channeloffsetend = {'  <defaultOffset>0</defaultOffset>';...
   ' </channelOffset>'};

chunk5 = {   '</channels>';...
 '</neuroscope>';...
'</parameters>'};

lineend = '\n';

% % basic scaffolding
% thischan = 0;
% s = chunk1;
%     s = cat(1,s,anatomygroupstart);
%         s = cat(1,s,[anatomychannelnumberline_start, num2str(thischan) anatomychannelnumberline_end]);
%     s = cat(1,s,anatomygroupend);
% s = cat(1,s,chunk2);
%     s = cat(1,s,spikegroupstart);
%         s = cat(1,s,[spikechannelnumberline_start, num2str(thischan) spikechannelnumberline_end]);
%     s = cat(1,s,spikegroupend);
% s = cat(1,s, chunk3);

%% Gather probe maps
if ~isempty(varargin)
    probemaplist = varargin;
else %make gui
    % mapfolder = fileparts(which 'NRX_Buzsaki64_8X8.xlsx');
    probemaplist = [];
end

grouplist_all = [];
groupchans_all = [];
channelcountoffset = 0;
for pmidx = 1:size(probemaplist,2)
    tpf = probemaplist{pmidx};
    if ~strcmp(probemaplist{pmidx}(end-4:end),'.xlsx')
        tpf = strcat(tpf,'.xlsx');
    end
    tpp = which(tpf);
    [tpnum,tptxt,tpraw] = xlsread(tpp);
    
    %find groups
    groupcolumn = strmatch('BY VERTI',tptxt(1,:));
    groupdenoterows = strmatch('SHANK ',tptxt(:,groupcolumn));
    groupperchannel = [];
    for ridx = 1:length(groupdenoterows);
       groupperchannel = cat(1,groupperchannel,str2num(tptxt{groupdenoterows(ridx),groupcolumn}(7:end))); 
    end
    grouplist_byprobe{pmidx} = unique(groupperchannel);
    grouplist_all = cat(1,grouplist_all,unique(groupperchannel));
    
    %find channels
    channelcolumn = groupcolumn+2;%may definitely to change this
%     groupcolumn = strmatch('Neuroscope channel',tptxt(1,:));
    tc = tpraw(groupdenoterows,channelcolumn);
    for cidx = 1:length(tc)
        channelnums(cidx,1) = tc{cidx};
    end
    
    %for each group, find the channels in it, save in sequence
    for gidx = 1:length(grouplist_byprobe{pmidx})
       tgidx = groupperchannel==grouplist_byprobe{pmidx}(gidx);%find rows whith this group/shank denotation 
       groupchans_byprobe{pmidx,gidx} = channelnums(tgidx)+channelcountoffset;
       groupchans_all = cat(1,groupchans_all,channelnums(tgidx))+channelcountoffset; 
    end
    numchansthisprobe = length(channelnums);
    channelcountoffset = numchansthisprobe;
    numchans = length(groupchans_all);
end

if isempty(probemaplist)
    groupchans_byprobe{1,1} = [0:numchans-1];
end

%% Make basic text 
s = chunk1;

s = cat(1,s,[channelcountlinestart, num2str(length(groupchans_all)) channelcountlineend]);

s = cat(1,s,chunk2);

%add channel count here

for pidx = 1:size(groupchans_byprobe,1)%for each probe
    for gidx = 1:size(groupchans_byprobe,2)%for each spike group
        s = cat(1,s,anatomygroupstart);
        tchanlist = groupchans_byprobe{pidx,gidx};
        for chidx = 1:length(tchanlist)
            thischan = tchanlist(chidx);
            s = cat(1,s,[anatomychannelnumberline_start, num2str(thischan) anatomychannelnumberline_end]);
        end
        s = cat(1,s,anatomygroupend);
    end
end

s = cat(1,s,chunk3);

for pidx = 1:size(groupchans_byprobe,1)%for each probe
    for gidx = 1:size(groupchans_byprobe,2)%for each spike group
        s = cat(1,s,spikegroupstart);
        tchanlist = groupchans_byprobe{pidx,gidx};
        for chidx = 1:length(tchanlist)
            thischan = tchanlist(chidx);
            s = cat(1,s,[spikechannelnumberline_start, num2str(thischan) spikechannelnumberline_end]);
        end
        s = cat(1,s,spikegroupend);
    end
end

s = cat(1,s, chunk4);

for pidx = 1:size(groupchans_byprobe,1)%for each probe
    for gidx = 1:size(groupchans_byprobe,2)%for each spike group
        tchanlist = groupchans_byprobe{pidx,gidx};
        for chidx = 1:length(tchanlist)
            thischan = tchanlist(chidx);
            s = cat(1,s,channelcolorstart);
            s = cat(1,s,[channelcolorlinestart, num2str(thischan) channelcolorlineend]);
            s = cat(1,s,channelcolorend);
            s = cat(1,s,channeloffsetstart);
            s = cat(1,s,[channeloffsetlinestart, num2str(thischan) channeloffsetlineend]);
            s = cat(1,s,channeloffsetend);
        end
    end
end

s = cat(1,s, chunk5);


%% Output
charcelltotext(s,'name.xml')


function charcelltotext(charcell,filename)
%based on matlab help.  Writes each row of the character cell (charcell) to a line of
%text in the filename specified by "filename".  Char should be a cell array 
%with format of a 1 column with many rows, each row with a single string of
%text.

[nrows,ncols]= size(charcell);

fid = fopen(filename, 'w');

for row=1:nrows
    fprintf(fid, '%s \n', charcell{row,:});
end

fclose(fid);
