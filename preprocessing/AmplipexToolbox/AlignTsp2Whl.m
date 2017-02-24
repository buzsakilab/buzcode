function AlignTsp2Whl(fbasename,colorvec,varargin)

% USAGE:
%     AlignTsp2Whl(fbasename,colorvec)
% 
% convert .tsp file to .whl file. Requires the .meta file
% 
% INPUT:
% 'fbasename': the file base names ('fbasename.tsp', etc.)
% 'colorvec': a 2 value vector [colFront colRear] which defines which
% color from [R G B] is the front and rear LEDs. Default: [1 3]
% optionnal:
%    AlignTsp2Whl(fbasename,colorvec,frommpg)
%    if frommpg is set to 1 (default 0), the tracking is done from the video
% 
% Adrien Peyrache and Antal Berényi, 2012
% This particular version is corrected by Adrien Peyrache, 2012.
% It also includes 2 new things: NaN instead of -1 for interpolation and
% video tracking for bad tsp files.
% Oct 2013: Brendon Watson made somewhat robust to instances of recording system
% crashes with blank .meta files

frommpg = 0;
if ~isempty(varargin)
    frommpg = varargin{1};
end

if strcmp('.tsp',fbasename(end-3:end));
    tspfile = fbasename;
    fbasename = fbasename(1:end-4);
else
   tspfile = [fbasename '.tsp'];
end

% get timestamps and positions from tsp file
if exist(tspfile,'file')
    tspdata=load(tspfile);
else
    error('Couldn''t do anything without the tsp file')
end

if size(tspdata,2)==1 || frommpg==1
    suffix = [];
    if exist([tspfile '.old'])
        suffix = 1;
        while exist([tspfile '.old. ' num2str(suffix)])
            suffix = suffix+1;
        end
        suffix = ['.' num2str(suffix)];
    end
    eval(['!cp ' tspfile ' ' tspfile '.old' suffix]);
    RecreateTspFile(fbasename);
    tspdata=load(tspfile);
end

%get length of dat file
%infoFile = dir([fbasename '.dat']);
%filelength=infoFile.bytes/(512*2);
%clear infoFile
%get start and end timestamp of dat file
EndTimestamp = [];%set these up so can detect if problem w10389984000ith getting them from .meta
StartTimestamp = [];
ChanNum = [];
DatSize = [];

fid=fopen([fbasename '.meta']);
tline= fgetl(fid);
while ischar(tline)
    try
    if strcmp(tline(1:20),'TimeStamp of the end')
        tline=tline(59:end);
        EndTimestamp=sscanf(tline,'%d',1);
    end
    end
    try
    if strcmp(tline(1:22),'TimeStamp of the start')
        tline=tline(61:end);
        StartTimestamp=sscanf(tline,'%d',1);
    end
    end
    try
    if strcmp(tline(1:9),'Number of')
        tline=tline(31:end);
        ChanNum=sscanf(tline,'%d',1);
    end
    end
    try
    if strcmp(tline(1:9),'File size')
        tline=tline(21:end);
        DatSize=sscanf(tline,'%lu',1);
    end
    end
    
    tline= fgetl(fid);
end
fclose(fid);

%making up for problems reading parameters (basically if .meta isn't good,
%but will address one at a time
if isempty(EndTimestamp)
    EndTimestamp = tspdata(end,1);
end

if isempty(StartTimestamp)
    StartTimestamp = tspdata(1,1);
end

if isempty(ChanNum);
    fid=fopen([fbasename '.ini']);
    tline= fgetl(fid);
    while ischar(tline)
        try
        if strcmp(tline(1:14),'ChannelsToSave')
            tline=tline(16:end);
            ChanNum = sum(str2num(tline));
        end
        end
        tline= fgetl(fid);
    end
    fclose(fid);
end

if isempty(DatSize)
    DatSize = dir([fbasename,'.dat']);
    DatSize = DatSize.bytes;
else
    DatSize = double(DatSize);
end



%Calculate file length from dat file sample rat
DatLength=DatSize/(ChanNum*2*20); %Dat file size in ms

% TspLength=EndTimestamp-StartTimestamp;  %Unused > problematic if crashes
% yielding blank .meta files

colorvec = sort([2*colorvec-1 2*colorvec])+1;
%remove lines from tspdata which has the same timestamp as the previous
%does
repeatingts=find(tspdata(1:end-1,1)==tspdata(2:end,1))+1;
for r=size(repeatingts):-1:1
    tspdata(repeatingts(r),:)=[];
end

%Putting all bad tracked points at NaN
tspdata(tspdata==-1)=NaN;

%interpolate to 1 kHz - computer clock
t=tspdata(1,1):tspdata(end,1);
interpTsp1kHz = zeros(length(t),7);
interpTsp1kHz(:,1)=t;
warning off %Doesn't like NaN but does well
interpTsp1kHz(:,2:end)=interp1(tspdata(:,1),tspdata(:,2:7),interpTsp1kHz(:,1));
warning on

for i=1:3%length(colorvec)
    minusonesegments(:,1)=tspdata(isnan(tspdata((2:end-1),i*2)),1);
    minusonesegments(:,2)=tspdata(find(isnan(tspdata((2:end-1),i*2)))+2,1);
    for j=1:size(minusonesegments,1)
        %Replace bad values by NaN so that interpolation ignores those
        %points.
        interpTsp1kHz(minusonesegments(j,1)-interpTsp1kHz(1,1)+1:minusonesegments(j,2)-interpTsp1kHz(1,1),i*2:i*2+1)=NaN;
    end
    clear minusonesegments
end

%align the beginning
if (StartTimestamp<tspdata(1,1)) 
    i=1:tspdata(1,1)-StartTimestamp;
    tspnew(i,1)=StartTimestamp+i-1;
    tspnew(i,2:7)=NaN;
    tspnew=[tspnew; interpTsp1kHz];
else
    tspnew=interpTsp1kHz(StartTimestamp-tspdata(1,1)+1:end,:);
end
clear i

% Tony's version
if (EndTimestamp<tspnew(end,1)) 
    tspnew=tspnew(1:find(tspnew(:,1)==EndTimestamp),:);
else
    t = tspnew(end,1)+1:EndTimestamp;
    %tspnew=[tspnew; [t' repmat(NaN,[DatLength-size(tspnew,1) length(colorvec)])]];
    tspnew=[tspnew; [t' NaN(size(t,2),6)]];
end

clear t interpTsp1kHz

t = tspnew(1,1):size(tspnew,1)/DatLength:tspnew(end,1);
interpTsp = zeros(length(t),7);
interpTsp(:,1)=tspnew(1,1):size(tspnew,1)/DatLength:tspnew(end,1);
warning off
interpTsp(interpTsp==-1)=NaN;
interpTsp(:,2:end)=interp1(tspnew(:,1),tspnew(:,2:end),interpTsp(:,1));
warning on
interpTsp(isnan(interpTsp))=-1;

if 0 %Do we want a whl1k file? Personnally, no...
    fid=fopen([file '.whl1k'],'w');
    for i=1:size(interpTsp,1)
        fprintf (fid,'%i %i\n',interpTsp(i,2:3));
    end
    disp(['Samples in .dat file per channel: ' int2str(DatLength)]);
    disp(['Lines in .whl1k file:' int2str(size(interpTsp,1))]);
    fclose(fid);
end
%resample to 39.0625 Hz
whldata=interpTsp(floor(1:25.6:size(interpTsp,1)),2:end);

if 0
    figure(1),clf
    hold on
    if    sum(colorvec == 1)
        scatter(whldata(:,1),whldata(:,2),'r.')
    end
    if    sum(colorvec == 2)
        scatter(whldata(:,3),whldata(:,4),'b.')
    end
    if    sum(colorvec == 3)
        scatter(whldata(:,5),whldata(:,6),'g.')
    end
    %pause
    %keyboard
end

%dlmwrite([fbasename '.whlall'],whldata,'delimiter','\t');
dlmwrite([fbasename '.whl'],whldata(:,colorvec-1),'delimiter','\t');

end
