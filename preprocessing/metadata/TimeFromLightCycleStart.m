function [SecondsAfterLightCycleStart_PerFile,SecondsAfterLightCycleStart] =  TimeFromLightCycleStart(basepath,lightsonhours,lightsonminutes,lightsonseconds)
% Zeitgeber times of recording files... based on file names and/or metadata
% files.  
% Inputs: 
% - basepath - per buzcode standard, where all files for a given session
% are kept
% - lightsonhours - hour of the day, in clock time when lights turn on
% for your animal, to determine zeitgeber time.  For example this is 6 if 
% light on at 6:34:15am.
% - lightsonminutes - hour of the day, in clock time when lights turn on
% for your animal, to determine zeitgeber time.  For example this is 34 if 
% light on at 6:34:15am.
% - lightsonseconds - hour of the day, in clock time when lights turn on
% for your animal, to determine zeitgeber time.  For example this is 15 if 
% light on at 6:34:15am.
% 
% Brendon Watson, 2016

%% Constants and settings

%% Input and directory handling 
if ~exist('basepath','var')
 basepath = cd;   
end
basename = bz_BasenameFromBasepath(basepath);

disp(basepath)
disp(basename)

savepath = fullfile(basepath,[basename '_SecondsFromLightsOn.mat']);

% lightson = 07:00:00;
lightsonhours_default = 6;
lightsonminutes_default = 0;
lightsonseconds_default = 0;
if ~exist('lightsonhours','var')
 lightsonhours = lightsonhours_default;   
end
if ~exist('lightsonminutes','var')
 lightsonminutes = lightsonminutes_default;   
end
if ~exist('lightsonseconds','var')
 lightsonseconds = lightsonseconds_default;   
end


%% First find out what system was used to record this... based on files present
% check if amplipex (.metas present) or intan (rhd.infos present)... or neither
recsys = '';

%Amplipex case
dmeta = dir(fullfile(basepath,'*.meta'));
if ~isempty(dmeta)
    recsys = 'Amplipex';
end

%Intan case
drhd = dir(basepath);
for didx = 1:length(drhd)
   if drhd(didx).isdir
       t = dir(fullfile(basepath,drhd(didx).name,'info.rhd'));
       if ~isempty(t)%only use folders with info.rhd inside
            recsys = 'Intan';
            break
       end
   end
end

%% Based on recording system, set up some basics for the function
switch recsys
    case 'Amplipex'%if from an amplipex
        bunderscore = strfind(basename,'_');
        basedate = basename(bunderscore+1:end);
        bdateday = str2num(basedate(3:4));
        bdateyear = 2000+str2num(basedate(5:6));
        bdatemonth = str2num(basedate(1:2));
    case 'Intan' %if not amplipex, look for intan
        bunderscore = strfind(basename,'_');
        basedate = basename(bunderscore(end)+1:end);
        bdateday = str2num(basedate(5:6));
        bdateyear = 2000+str2num(basedate(1:2));
        bdatemonth = str2num(basedate(3:4));
end

lightsondatenum = datenum(bdateyear,bdatemonth,bdateday);

switch recsys
    case 'Amplipex'%if from an amplipex
       d = dir(fullfile(basepath,[basename(1:end-3) '*' '.meta']));%if I record over new year's eve I'll have to handle it :)
       fname1 = fullfile(basepath,[basename '-01.meta']);%explicitly look for -01
       seconds = getmetafilestarttime(fname1,lightsondatenum,lightsonhours,lightsonminutes,lightsonseconds);
       SecondsAfterLightCycleStart = seconds;
       SecondsAfterLightCycleStart_PerFile = nan(1,length(d));
       for a = 1:length(d)
           fname = fullfile(basepath,d(a).name);  
           seconds = getmetafilestarttime(fname,lightsondatenum,lightsonhours,lightsonminutes,lightsonseconds);
           SecondsAfterLightCycleStart_PerFile(a) = seconds;
       end
    case 'Intan' %if not amplipex, look for intan
        d = dir(fullfile(basepath,[basename(1:end-3) '*']));%if I record over new year's eve I'll have to handle it :)
        for a = length(d):-1:1
            if ~d(a).isdir
                d(a) = [];
            end
        end
        SecondsAfterLightCycleStart_PerFile = nan(1,length(d));
        for a = 1:length(d)
            seconds = getintanfilestarttime(d(a).name,lightsondatenum,lightsonhours,lightsonminutes,lightsonseconds);
            if a == 1 
               SecondsAfterLightCycleStart = seconds;
            end
            SecondsAfterLightCycleStart_PerFile(a) = seconds;
        end
end
    

save(savepath,'SecondsAfterLightCycleStart','SecondsAfterLightCycleStart_PerFile')
% SecondsAfterLightCycleStart_PerFile

function seconds = getmetafilestarttime(fname,lightsondatenum,lightsonhours,lightsonminutes,lightsonseconds)

try
    ttime = ReadMetaAspects(fname,'starttime');
    h = str2num(ttime(1:2));
    m = str2num(ttime(4:5));
    s = str2num(ttime(7:8));
    dashes = strfind(fname,'-');

    tdate = ReadMetaAspects(fname,'startdate');
    spaces = strfind(tdate,' ');
    datemonth = tdate(1:spaces(1)-1);
    dateday = tdate(spaces(1)+1:spaces(2)-1);
    if length(dateday) == 1;
        dateday = ['0' dateday];
    end
    dateyear = tdate(spaces(2)+1:end);
    dn = datenum([dateday '-' datemonth '-' dateyear],'dd-mmm-yyyy');

    hoursecs = (h-lightsonhours)*3600;
    minsecs = (m-lightsonminutes)*60;
    secs = s-lightsonseconds;
    daysecs = 24*3600*(dn-lightsondatenum);

    seconds = daysecs+hoursecs+minsecs+secs;
catch 
    disp(['reading time/date from ' fname ' did not work'])
    seconds = NaN;
end



function seconds = getintanfilestarttime(fname,lightsondatenum,lightsonhours,lightsonminutes,lightsonseconds)

underscores = strfind(fname,'_');
ttime = fname(underscores(end)+1:end);
h = str2num(ttime(1:2));
m = str2num(ttime(3:4));
s = str2num(ttime(5:6));

tdate = fname(underscores(end-1)+1:underscores(end)-1);
dateday = tdate(5:6);
dateyear = num2str(2000+str2num(tdate(1:2)));
datemonth = tdate(3:4);
dn = datenum([dateday '-' datemonth '-' dateyear],'dd-mmm-yyyy');

hoursecs = (h-lightsonhours)*3600;
minsecs = (m-lightsonminutes)*60;
secs = s-lightsonseconds;
daysecs = 24*3600*(dn-lightsondatenum);

seconds = daysecs+hoursecs+minsecs+secs;