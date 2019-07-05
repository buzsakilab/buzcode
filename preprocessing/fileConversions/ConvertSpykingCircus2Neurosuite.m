function ConvertSpykingCircus2Neurosuite(basepath,basename)
% Converts hdf5 files from SpykingCircus into klusters-compatible
% fet,res,clu,spk files.  Works on all shanks of a recording.  Assumes a
% 16bit .dat and an .xml file is present in "basepath" (home folder) and 
% that they are named basename.dat and basename.xml.  Also assumes that
% subdirectories in that basepath are made for each shank with names
% specified by numbers (ie 1,2,3,4..8).  In each shank folder should be
% .kwik and .kwx files made by klusta with names as follows:
% basename_sh[shankumber].kwik/kwx.  This 
% 
% Inputs:
% basepath - directory path to the main recording folder with .dat and .xml
% as well as shank folders made by makeProbeMapKlusta2.m (default is
% current directory matlab is pointed to)
% basename - shared file name of .dat and .xml (default is last part of
% current directory path, ie most immediate folder name)
%
% Brendon Watson 2016

if ~exist('basepath','var');
    [~,basename] = fileparts(cd);
    basepath = cd;
end

datpath = fullfile(basepath,[basename '.dat']);
xmlpath = fullfile(basepath,[basename '.xml']);
spydir = fullfile(basepath,basename);

% move some files to take care of issues with numeric names being rejected
% by h5 functions
d = dir(fullfile(spydir,'*.result.hdf5'));
oldresultspath = fullfile(spydir,d(end).name);
newresultspath = fullfile(spydir,'result.result.hdf5');
eval(['!cp ' oldresultspath ' ' newresultspath])

d = dir(fullfile(spydir,'*.templates.hdf5'));
oldtemplatespath = fullfile(spydir,d(end).name);
newtemplatespath = fullfile(spydir,'templates.templates.hdf5');
eval(['!cp ' oldtemplatespath ' ' newtemplatespath])

%% Get basic channel info
par = LoadPar(xmlpath);
% par = LoadXml(xmlpath);
totalch = par.nChannels;
% sbefore = 17;%samples before/after for spike extraction
% safter = 16;%... could read from SpkGroups in xml - was originally 16
if isfield(par.SpkGrps,'nSamples')
    if ~isempty(par.SpkGrps(1).nSamples);
        if isfield(par.SpkGrps,'PeakSample')
            if ~isempty(par.SpkGrps(1).PeakSample);
                sbefore = par.SpkGrps(1).PeakSample;
                safter = par.SpkGrps(1).nSamples - par.SpkGrps(1).PeakSample;
            end
        end
    end
end

% sbefore = 16;
% safter = 16;



grouplookup = zeros(totalch,1);
% for a= 1:par.nElecGps
for a = 1:length(par.SpkGrps) % fixed
    grouplookup(par.SpkGrps(a).Channels+1) = a;
end
allgroups = unique(grouplookup);

%Grp 0 contain discared channels
allgroups(allgroups==0) = [];

%% get number of clusters
t = h5info(newresultspath);
for a = 1:length(t.Groups);
    n{a} = t.Groups(a).Name;
end
n = strmatch('/spiketimes',n);
numclus = size(t.Groups(n).Datasets,1);

%% get all spike times
clutimes = [];
clugrps = [];
for a = 1:numclus
    eval(['clutimecell{' num2str(a) '}  = double(h5read(newresultspath,''/spiketimes/temp_' num2str(a-1) '''));'])
    clutimes = cat(1,clutimes,clutimecell{a});
    clugrps = cat(1,clugrps,a*ones(size(clutimecell{a})));
end
[spktimes,sidx] = sort(clutimes);
clu = clugrps(sidx);

%% get templates
templates_size = double(h5read(newtemplatespath, '/temp_shape'));
N_e = templates_size(2);
N_t = templates_size(1);
temp_x = double(h5read(newtemplatespath, '/temp_x') + 1);
temp_y = double(h5read(newtemplatespath, '/temp_y') + 1);
temp_z = double(h5read(newtemplatespath, '/temp_data'));
tmps = sparse(temp_x, temp_y, temp_z, templates_size(1)*templates_size(2), templates_size(3));
templates_size = [templates_size(1) templates_size(2) templates_size(3)/2];
for a = 1:numclus
    templates(:,:,a) = full(reshape(tmps(:, a), templates_size(2), templates_size(1)));
end

%% assign templates to shanks
m = max(abs(templates),[],1);%find the most deviated value of each waveform on each channel
[~,m] = max(m,[],2);%find which channel has most deviated value for each templnate
m = squeeze(m);%squeeze to 1d vector

templateshankassignments = grouplookup(m);%for the list of maximal channels, which group is each in 

%% write shank-wise information
for groupidx = 1:length(allgroups)
    tgroup          = allgroups(groupidx);%shank number
    ttemplateidxs   = find(templateshankassignments==tgroup);%which templates/clusters are in that shank

    tidx            = ismember(clu,ttemplateidxs);%find spikes indices in this shank
    tclu            = clu(tidx);%extract template/cluster assignments of spikes on this shank
    tspktimes       = spktimes(tidx);
    
    %% 
    gidx            = find(grouplookup == tgroup);%find all channels in this group
    channellist     = [];
    for ch = 1:length(par.SpkGrps)
        if ismember(gidx(1),par.SpkGrps(ch).Channels+1)
            channellist = par.SpkGrps(ch).Channels+1;
            break
        end
    end
    if isempty(channellist)
        disp(['Cannot find spkgroup for group ' num2str(groupidx) ])
        continue
    end
    
%     %% spike extraction from dat - old version
%     if groupidx == 1
%         dat = memmapfile(datpath,'Format','int16');
%     end
%     tsampsperwave   = (sbefore+safter);
%     ngroupchans     = length(channellist);
%     valsperwave     = tsampsperwave * ngroupchans;
%     wvforms_all     = zeros(length(tspktimes)*tsampsperwave*ngroupchans,1,'int16');
% %     wvranges        = zeros(length(tspktimes),ngroupchans);
% %     wvpowers        = zeros(1,length(tspktimes));
%     
%     for j=1:length(tspktimes)
%         try
%             w = dat.data((double(tspktimes(j))-sbefore).*totalch+1:(double(tspktimes(j))+safter).*totalch);
%             wvforms=reshape(w,totalch,[]);
%             %select needed channels
%             wvforms = wvforms(channellist,:);
%     %         % detrend
%     %         wvforms = floor(detrend(double(wvforms)));
%             % median subtract
%             wvforms = wvforms - repmat(median(wvforms, 2),1,sbefore+safter);
%             wvforms = wvforms(:);
%             
%         catch
%             disp(['Error extracting spike at sample ' int2str(double(tspktimes(j))) '. Saving as zeros']);
%             disp(['Time range of that spike was: ' num2str(double(tspktimes(j))-sbefore) ' to ' num2str(double(tspktimes(j))+safter) ' samples'])
%             wvforms = zeros(valsperwave,1);
%         end
% 
%         %some processing for fet file
% %         wvaswv = reshape(wvforms,tsampsperwave,ngroupchans);
% %         wvranges(j,:) = range(wvaswv);
% %         wvpowers(j) = sum(sum(wvaswv.^2));
% 
%         lastpoint = tsampsperwave*ngroupchans*(j-1);
%         wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms';
%     %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
%         if rem(j,100000) == 0
%             disp([num2str(j) ' out of ' num2str(length(tspktimes)) ' done'])
%         end
%     end

        %% spike extraction from dat - new version
    dat=memmapfile(datpath,'Format','int16');
    tsampsperwave = (sbefore+safter);
    ngroupchans = length(channellist);
    valsperwave = tsampsperwave * ngroupchans;
    wvforms_all=zeros(length(tspktimes)*tsampsperwave*ngroupchans,1,'int16');
    wvranges = zeros(length(tspktimes),ngroupchans);
    wvpowers = zeros(1,length(tspktimes));
    for j=1:length(tspktimes)
        try
            w = dat.data((double(tspktimes(j))-sbefore).*totalch+1:(double(tspktimes(j))+safter).*totalch);
            wvforms=reshape(w,totalch,[]);
            %select needed channels
            wvforms = wvforms(channellist,:);
    %         % detrend
    %         wvforms = floor(detrend(double(wvforms)));
            % median subtract
            wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
            wvforms = wvforms(:);
        catch
            disp(['Error extracting spike at sample ' int2str(double(tspktimes(j))) '. Saving as zeros']);
            disp(['Time range of that spike was: ' num2str(double(tspktimes(j))-sbefore) ' to ' num2str(double(tspktimes(j))+safter) ' samples'])
            wvforms = zeros(valsperwave,1);
        end

        %some processing for fet file
        wvaswv = reshape(wvforms,tsampsperwave,ngroupchans);
        wvranges(j,:) = range(wvaswv);
        wvpowers(j) = sum(sum(wvaswv.^2));

        lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
    %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,50000) == 0
            disp([num2str(j) ' out of ' num2str(length(tspktimes)) ' done'])
        end
    end
    clear dat
    wvranges = wvranges';
    
    
    
    
    %% not doing fets for now
    
    %% writing to clu, res, fet, spk
    cluname = fullfile(basepath,[basename '.clu.' num2str(tgroup)]);
    resname = fullfile(basepath,[basename '.res.' num2str(tgroup)]);
%     fetname = fullfile(basepath,[basename '.fet.' num2str(tgroup)]);
    spkname = fullfile(basepath,[basename '.spk.' num2str(tgroup)]);
  %fet
%     SaveFetIn(fetname,fets);

    %clu
    tclu = cat(1,length(unique(tclu)),double(tclu));
    fid=fopen(cluname,'w'); 
    fprintf(fid,'%.0f\n',tclu);
    fclose(fid);
    clear fid

    %res
    fid=fopen(resname,'w'); 
    fprintf(fid,'%.0f\n',tspktimes);
    fclose(fid);
    clear fid

    %spk
    fid=fopen(spkname,'w'); 
    fwrite(fid,wvforms_all,'int16');
    fclose(fid);
    clear fid 

    disp(['Shank ' num2str(tgroup) ' done'])
    %end
    %end    
    
    
end
clear dat

%% Fets made the traditional way here
MakeClassicFet(basename,basepath)

%% Save original clus away, in case
mkdir(fullfile(basepath,'OriginalClus'))
copyfile(fullfile(basepath,[basename,'.clu.*']),fullfile(basepath,'OriginalClus'))


function Par = LoadParIn(FileName, varargin)
% LoadPar(FileName)
% loads the specified par file and returns a structure with these elements:
%
% .FileName      -> name of file loaded from
% .nChannels     -> number of total channels
% .nBits         -> number of bits of the file
% .SampleTime    -> time, in microseconds, of 1 sample (ie 1e6 / sample rate)
% .HiPassFreq    -> High pass filter frequency
% .nElecGps      -> number of electrodes (i.e. electrode groups)
% .ElecGp        -> a cell array giving the channels in the electrodes
%                    e.g. if .ElectrodeGroup{3} = [2 3 4 5], electrode 3
%                    is a tetrode for channels 2 3 4 and 5. 
% channel numbers here are from 0. be carefull.
[SpecInfo] = DefaultArgsIn(varargin,{1});

if ~isempty(strfind(FileName,'.par'))
    FileBase = FileName(1:strfind(FileName,'.par')-1);
elseif ~isempty(strfind(FileName,'.xml'))
    FileBase = FileName(1:strfind(FileName,'.xml')-1);
else 
    FileBase = FileName;
end


if exist([FileBase '.xml'],'file') %& ~isempty(strfind(FileName,'.xml'))
    Par = LoadXmlIn(FileBase);
elseif exist([FileBase '.par'],'file')


    % open file

    fp = fopen([FileBase '.par'], 'r');
    Par.FileName = FileBase;

    % read in nChannels and nBits
    Line = fgets(fp);
    A = sscanf(Line, '%d %d');
    Par.nChannels = A(1);
    Par.nBits = A(2);

    % read in SampleTime and HiPassFreq
    Line = fgets(fp);
    A = sscanf(Line, '%d %f', 2);
    Par.SampleTime = A(1);
    Par.HiPassFreq = A(2);

    % read in nElectrodes
    Line = fgets(fp);
    if Line==-1
        fclose(fp);
        return;
    end
    A = sscanf(Line, '%d', 1);
    Par.nElecGps = A(1);

    % read in ElectrodeGroup
    for i=1:Par.nElecGps
        Line = fgets(fp);
        A = sscanf(Line, '%d');
        Par.ElecGp{i} = A(2:end);
    end;
    fclose(fp);
else
    error('Par or Xml file do not exist!');
end

if SpecInfo
    if FileExists([FileBase '.eeg.par'])
        if ~isfield(Par,'nElecGps')
            ParTmp = LoadPar([FileBase '.par']);
        else
            ParTmp = Par;
        end

        EegPar=LoadEegPar(FileBase);
        for el=1:ParTmp.nElecGps
            for eegel=1:EegPar.nElec
                if ~isempty(intersect(ParTmp.ElecGp{el},EegPar.ElecChannels{eegel}))
                    Par.ElecLoc{el} = EegPar.ElecLoc{eegel};
                end
            end
        end
   
    end
end

function varargout = DefaultArgsIn(Args, DefArgs)
% auxillary function to replace argument check in the beginning and def. args assigment
% sets the absent or empty values of the Args (cell array, usually varargin)
% to their default values from the cell array DefArgs. 
% Output should contain the actuall names of arguments that you use in the function

% e.g. : in function MyFunction(somearguments , varargin)
% calling [SampleRate, BinSize] = DefaultArgs(varargin, {20000, 20});
% will assign the defualt values to SampleRate and BinSize arguments if they
% are empty or absent in the varargin cell list 
% (not passed to a function or passed empty)
if isempty(Args)
    Args ={[]};
end

% if iscell(Args) & isstr(Args{1}) & length(Args)==1
%     Args = Args{1};
% end
    
if ~iscell(DefArgs)
    DefArgs = {DefArgs};
end
nDefArgs = length(DefArgs);
nInArgs = length(Args);
%out = cell(nDefArgs,1);
if (nargout~=nDefArgs)
    error('number of defaults is different from assigned');
    %keyboard
end
for i=1:nDefArgs
    
    if (i>nInArgs | isempty(Args{i}))
        varargout(i) = {DefArgs{i}};
    else 
        varargout(i) = {Args{i}};
    end
end

%function [xml, rxml] = LoadXml_SleepScore(FileBase)
%loads the xml file using xmltools_ss (have to have it in the path)
% rxml returns it's original layout - very messy structure but contains all
% the xml file contents.
% xml - is the ouput structure which is backwards compatible to LoadPar
% output, so you can use it instead ..also loads some usefull stuff -
% Anatomoical groups with Skips , Spike electrode groups
% more can be added later (e.g. parameters of the process scripts)
% this script is written for xml version 1.1 .. older version doesn't work.
% additions are welcome
% 
% By default, fbasename is the name of the directory where the data are, so
% the exact xml file path should not be put in argument
% By calling
% xml = LoadXml(xmlfile,'raw')
% the xml file specified in argument will be loaded
% (added by A Peyrache)

function [xml, rxml] = LoadXmlIn(fbasename,varargin)
xml = struct;

%if isempty(varargin)
%    [fbasename mergedir rootdir] = extractfbasename(fbasename);
%elseif strcmpi(varargin,'raw')
%    fprintf('Loading directly xml file')
%else
%    error('Unrecognized option')
%end
%if xml was in the filebase by chance
xmli = strfind(fbasename,'.xml');
if ~isempty(xmli)
    fbasename = fbasename(1:xmli-1);
end
rxml = xmltoolsIn([fbasename '.xml']);

rxml = rxml.child(2);

% from this level all children are the different parameters fields

xml.FileName = fbasename;

for i=1:length(rxml.child)

    switch rxml.child(i).tag
        
        case 'generalInfo'
            xml.Date = rxml.child(i).child(1).value; % date of xml file creation?

        case 'acquisitionSystem'
            xml.nBits = str2num(rxml.child(i).child(1).value); % number of bits of the file
            xml.nChannels = str2num(rxml.child(i).child(2).value);
            xml.SampleRate = str2num(rxml.child(i).child(3).value);
            xml.SampleTime = 1e6/xml.SampleRate; %to make backwards compatible
            xml.VoltageRange = str2num(rxml.child(i).child(4).value);
            xml.Amplification = str2num(rxml.child(i).child(5).value);
            xml.Offset =  str2num(rxml.child(i).child(6).value);
            
        case 'fieldPotentials'
            xml.lfpSampleRate = str2num(rxml.child(i).child.value);
            
        case 'anatomicalDescription'
            tmp = rxml.child(i).child.child;
            for grpI =1:length(tmp)
                for chI=1:length(tmp(grpI).child)
                    xml.AnatGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(chI).value);
                    xml.AnatGrps(grpI).Skip(chI) = str2num(tmp(grpI).child(chI).attribs.value);
                end
            end
            
        case 'spikeDetection'
            if ~isempty(rxml.child(i).child)
                tmp =rxml.child(i).child.child;
                for grpI =1:length(tmp)
                    for chI=1:length(tmp(grpI).child(1).child)
                        xml.SpkGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(1).child(chI).value);
                    end
                    if length(tmp(grpI).child)>1
                        xml.SpkGrps(grpI).nSamples = str2num(tmp(grpI).child(2).value);
                        xml.SpkGrps(grpI).PeakSample = str2num(tmp(grpI).child(3).value);
                        xml.SpkGrps(grpI).nFeatures = str2num(tmp(grpI).child(4).value);
                    end
                    %backwards compatibility
                    xml.nElecGps = length(tmp);
                    xml.ElecGp{grpI} = xml.SpkGrps(grpI).Channels;
                end
            else
                xml.nElecGps = 0;
            end


        case 'programs'
            tmp = rxml.child(i).child;
            for i=1:length(tmp)
                if strcmp(tmp(i).child(1).value,'process_mhipass')
                    for j=1:length(tmp(i).child(2).child )
                        if strcmp(tmp(i).child(2).child(j).child(1).value,'frequency')
                            xml.HiPassFreq = str2num(tmp(i).child(2).child(j).child(2).value);
                            break
                        end
                    end
                end
            end
    end
end


% general recursive parsing will have to wait.
function z = xmltoolsIn( arg, out_file, varargin)
% xmltools_ss - tools for managing xml data sets
%      - if arg is a string : arg is an XML file to convert into MATLAB struct
%      - if arg is a variable : it is a MATLAB struct to write into XML, to stdout if out_file is not given
% use :
%  z = xmltools_ss('filename.xml'); read an xml file and store it into z
%  xmltools_ss(z,'filename.xml'); write z into the file
%  xmltools_ss(z,'get','tag-name'); returns only subset of z child which name is tag-name
%  xmltools_ss(z,'get-attrib', 'attrib-name')
%  xmltools_ss(z,'get','tag-name', 'attribs'|'value'|'child');
%
% project 'File parsing'
% title    'XML parsing'
% author  'Charles-Albert Lehalle'
% mailto  'charles.lehalle@miriadtech.com'
% version '2.5'
% date    'mar2003--sept2003'

version = '2.5';

%%** XML TOOLS FOR MATLAB
% This is an OCAMAWEB (http://ocamaweb.sourceforge.net) generated documentation.\\
% This function manage the exchange of data between XML files and a MATLAB structure.
% The MATLAB structure is this : a {\bf node} is a struct with fields {\it tag}, {\it value}, {\it attribs}, {\it child}.
% Where {\it tag} is the tag name, {\it value} its contents, {\it attribs} an array of structure with fields {\it name} and {\it value},
% and {\it child} an array of such nodes.\\
% All those fields are always present and can be empty except for the first node (root) that has only children.\\
%

global Verbose;
Verbose=0;
%%* READ AN XML FILE
if isstr(arg)
  
  %< Rï¿½cupï¿½ration du fichier dans un string
  fid = fopen(arg, 'r');
  F = fread(fid);
  s = char(F');
  fclose(fid);
  %>
  
  %< Parsing
  z = parse_xml(s);
  %>
  return
end

if ~isstr(arg)

  %<* SELECT A SUBSET OF z CHILD
  if length(varargin) >= 1 %& ~FileExists_ss(arg)
    
    cmode = upper( out_file);
    z     = arg;
    
    switch upper( cmode)
      
     case 'GET'
      %< Get the subset
      % warning: I will have to change the value of next at some places
      next = 'child';
    
      if ~isfield(z, next)
        error('xmltools_ss:GET', 'For child selection, structured first argument is needed');
      end
      tag_name = (varargin{1});
      
      z = get_childs(z, next, tag_name);
      %>
      
      %< get values
      % xmltools_ss( z, 'get', 'tag-name', 'attribs'|'value'|'child')
      if length(varargin) == 2
        switch upper( varargin{2})
         case 'VALUE'
          z = [z.value];
         case 'ATTRIBS'
          z = [z.attribs];
         case 'CHILD'
          z = [z.child];
         otherwise
          error('xmltools_ss:GET', 'get mode <%s> is not defined, use one of <attribs>|<value>|<child>', upper(varargin{2}));
        end
      end
      %>
      
     case 'GET-ATTRIB'
      %< get attrib
      % xmltools_ss(z, 'get-attrib', 'attrib-name')
      s = z.attribs;
      for i=1:length(s)
        if strcmp((s(i).name), ( varargin{1}) )
          z = s(i).value;
          return
        end
      end
      error('xmltools_ss:GET-ATTRIB', 'no attribute found'); %, varargin{1});
      %>
      
    end
           
    return
  end
  %>*
  
  %<* WRITE AN XML STRUCTURE
  
  %< Selection de la cible
  if nargin < 2
    fid = 1;
  else
    fid = fopen(out_file, 'w');
  end
  %>
  
  %< Ecriture proprement dite
  write_xml(fid, arg, -1);
  %>
  
  %< Fermeture
  if nargin > 1
    fclose(fid);
  end
  %>
  
  %>*
end

%%** Fonctions internes

%<* parser un string xml
function [z, str] = parse_xml( str, current_tag, current_value, attribs, idx)
global Verbose;
next = 'child';

if nargin < 2
  current_tag   = '';
  current_value = '';
  attribs       = '';
  idx           = 0;
end
z = [];

eot = 0;

while ~eot & ~isempty(udeblank(deblank(str)))
  
  f_end = strfind(str, '</');
  f_beg = strfind(str, '<');
  
  %< Si je n'ai plus de tag dans mon document
  if isempty(f_end) & isempty(f_beg)
    
    if ~strcmp(lower(current_tag), '?xml') & ~isempty(current_tag)
      error('xmltools_ss:parse_xml', 'malformed xml string (current [%s])', current_tag);
    else
      fprintf('end parsing at level %d\n',idx);
      eot = 1;
      return
    end
  end
  %>
  
  if isempty(f_end)
    f_end = length(str)
  else
    f_end = f_end(1);
  end
  if isempty(f_beg)
    f_beg = length(str)
  else
    f_beg = f_beg(1);
  end
  
  if f_end <= f_beg
    %< je rencontre une fermeture
    new_tag = str((f_end+2):end);
    str_t   = str(1:f_end-1);
    f_end = strfind(new_tag,'>');
    if isempty(f_end)
      error('xmltools_ss:parse_xml', 'malformed xml string : never ending tag [%s] encountered', current_tag);
    end
    f_end = f_end(1);
    str     = new_tag(f_end+1:end); % reste
    new_tag = new_tag(1:f_end-1);
    if ~strcmp((new_tag), (current_tag))
      error('xmltools_ss:parse_xml', 'malformed xml string : [%s] not properly closed (closing [%s] encountered)', current_tag, new_tag);
    end
    if Verbose
       fprintf('%sclose [%s]\n', repmat(' ', 2*(idx-1),1), current_tag);
    end
    z.tag     = (current_tag);
    z.attribs = parse_attribs(attribs);
    z.value   = udeblank(deblank(sprintf('%s %s',current_value, str_t)));
    eot       = 1;
    %>
  else
    %< je rencontre une ouverture
    % je vais appeler le mï¿½me code sur ce qu'il y a aprï¿½s moi
    current_value = sprintf('%s %s', current_value, str(1:f_beg-1));
    new_tag   = str(f_beg+1:end);
    f_end = strfind(new_tag,'>');
    if isempty(f_end)
      error('xmltools_ss:parse_xml', 'malformed xml string : never ending tag encountered');
    end
    f_end   = f_end(1);
    str_t   = new_tag(f_end+1:end);
    new_tag = new_tag(1:f_end-1);
    if (new_tag(end) == '/')|(new_tag(end) == '?')
      %< Self closing tag
      % Je met (temporairement!) eot ï¿½ 1, cela me permet de passer quelques lignes
      % de code tranquilement
      eot = 1;
      %>
    end
    %< Attributs
    f_beg   = strfind(new_tag, ' ');
    if isempty(f_beg)
      new_attribs = '';
      if eot
	new_tag = new_tag(1:end-1);
      end
    else
      new_attribs = new_tag(f_beg+1:end);
      if eot
	new_attribs = new_attribs(1:end-1);
      end
      new_tag     = new_tag(1:f_beg-1);
    end
    %>
    if Verbose
        fprintf('%sopen  [%s]\n', repmat(' ', 2*idx,1), new_tag);
    end
    if eot
      %< If self-colsing tag
      if Verbose
          fprintf('%sclose [%s]\n', repmat(' ', 2*idx,1), new_tag);
      end
      new_attribs = parse_attribs( new_attribs);
      if isfield(z, next)
        nxt = getfield(z, next);
        nxt(end+1) = struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []);
        z   = setfield(z, next, nxt);
	%z.(next)(end+1) = struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []);
      else
        z = setfield(z, next, struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []) );
	%z.(next) = struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []);
      end
      str = str_t;
      eot = 0;
      %>
    else
      %< Appel du mï¿½me code sur la suite

      % et stockage du resultat dans mes children.
      % Le code met aussi ï¿½ jour le string courant |str|,
      % il en enlï¿½ve la partie correspondant au string que je viens de trouver.
      [t,str] = parse_xml(str_t, new_tag, '', new_attribs, 1+idx);
      if isfield(t, next)
	nx = getfield( t, next);
        %nx = t.(next);
      else
	nx = [];
      end
      if isfield(z, next)
        nxt = getfield(z, next);
        nxt(end+1) = struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx);
        z   = setfield(z, next, nxt);
	%z.(next)(end+1) = struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx);
      else
	z = setfield(z, next, struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx));
	%z.(next) = struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx);
      end

      %>
    end
  end
  %>
end
%>

%< Parse attribs
function z =  parse_attribs( a)
if isempty(a)
  z = struct( 'name', '', 'value', '');
  return
end
b = tokens(a, ' ');
j = 1;
   
for i=1:length(b)
  if ~isempty(b{i})
    t = tokens(b{i}, '=');
    if length(t)==2
      u = t{2};
      if u(1)=='"'
	u = u(2:end);
      end
      if u(end)=='"'
	u = u(1:end-1);
      end
      z(j) = struct( 'name', (t{1}), 'value', u);
    else
        %here is the possible problem: if there is space between name and
        %'=' and value this will be counted 3 times as different names and
        %empty value
      z(j) = struct( 'name', (a), 'value', '');
    end
    j = j +1;
  end
end
%>


%<* Ecriture d'une structure xml
function z = write_xml(fid, xml_struct, idx)

next = 'child';

if nargin < 3
  idx = 0;
end

margin = repmat(' ',idx,1);

closed_tag = 1;
%< Ouverture du tag
if isfield(xml_struct, 'tag')
  closed_tag = 0;
  fprintf(fid, '%s<%s', margin, xml_struct.tag);
  %< Ecriture des attributs
  if ~isfield(xml_struct, 'attribs')
    error('xmltools_ss:write_xml', 'malformed MATLAB xml structure : tag without attribs');
  end
  for i=1:length(xml_struct.attribs)
    if ~isempty(xml_struct.attribs(i).name)
        if xml_struct.tag(1) == '?'
            fprintf(fid, ' %s=%s', xml_struct.attribs(i).name, xml_struct.attribs(i).value);
        else
          fprintf(fid, ' %s="%s"', xml_struct.attribs(i).name, xml_struct.attribs(i).value);
        end
    end
  end
  %>
  
  %< Gestion des Auto closed tags
  % Si le tag n'est pas auto fermï¿½, alors |closed_tag| est ï¿½ zï¿½ro
  if ~isfield(xml_struct, next)
    error('xmltools_ss:write_xml', 'malformed MATLAB xml structure : tag without %s', next);
  end
  if ~isfield(xml_struct, 'value')
    error('xmltools_ss:write_xml', 'malformed MATLAB xml structure : tag without value');
  end
  if xml_struct.tag(1) == '?'
    fprintf(fid, '?>\n');
    closed_tag = 1;
  elseif isempty(getfield(xml_struct, next)) & isempty(xml_struct.value)
  %elseif isempty(xml_struct.(next)) & isempty(xml_struct.value)
    fprintf(fid, '/>\n');
    closed_tag = 1;
  elseif ~isempty(xml_struct.value)
      fprintf(fid, '>');
  else
      fprintf(fid, '>\n');
  end
  %>
end
%>

%< Ecriture de la value
if isfield(xml_struct, 'value')
  if ~isempty(xml_struct.value)
    fprintf(fid, '%s',  xml_struct.value);
  end
end
%>

%< Ecriture des enfants
if ~isfield(xml_struct, next)
  error('xmltools_ss:write_xml', 'malformed MATLAB xml structure : tag without %s', next);
end
those_children = getfield(xml_struct, next);
%those_children = xml_struct.(next);

for i=1:length(those_children)
  write_xml(fid, those_children(i), idx+1);
end
%>

%< Fermeture du tag
if ~closed_tag
    if ~isempty(xml_struct.value)
        fprintf(fid, '</%s>\n',  xml_struct.tag);
    else
      fprintf(fid, '%s</%s>\n',  margin,xml_struct.tag);
    end
end
%>
%>*


%<* get childs with a specific tag name
function z = get_childs(z, next, tag_name);
u = getfield(z, next);
zo = [];
for i=1:length(u)
  v = u(i);
  if strcmp((v.tag), (tag_name))
    if isempty(zo)
      zo.anext= v;
    else
      zo.anext(end+1) = v;
    end
  end
end
if ~isstruct( zo)
  if isfield(z, 'tag')
    tn = z.tag;
  else
    tn = 'root?';
  end
  error('xmltools_ss:GET-TEG', 'problem in finding tag <%s> under one <%s>', tag_name, tn);
end
z = [ zo.anext ];
%>*

%< udeblank
function s = udeblank(str)
s = deblank(str(end:-1:1));
s = s(end:-1:1);
if length(s)==0
  s = '';
end
%>

%< emptystruct
function z = emptystruct(next)
z = struct( 'tag', [], 'value', [], 'attribs', [], next, []);
%>

%< Tokens
function l = tokens(str,del)
l={} ; 
% Boucle sur les tokens.
del = sprintf(del) ;

space_pos = strfind(str,del);
eq_pos = strfind(str,'=');
index2squash = [];
for k=1:length(eq_pos)
    if sum(ismember(eq_pos(k)+[-1 1], space_pos))==2
        index2squash = [index2squash eq_pos(k)+[-1 1]];
    end
end
str(index2squash)='';

while ~isempty(str)
  [tok,str] = strtok(str,del) ;
  l{end+1} = tok ;
end

%fixing possible spaces btw '=' and name/value
% if length(l)==3 & strcmp(l{2},'=')
%     keyboard
% end

%>
