function Par = LoadParXML(FileName)
% LoadPar(FileName)
% loads the specified par file and returns a structure with these elements:
%
% .FileName      -> name of file loaded from
% .nChannels     -> number of total channels
% .nBits         -> number of bits of the file
% .SampleTime    -> time, in microseconds, of 1 sample (ie 1e6 / sample rate)
% .HiPassFreq    -> High pass filter frequency
% .nElecGps   -> number of electrodes (i.e. electrode groups)
% .ElecGp -> a cell array giving the channels in the electrodes
%                    e.g. if .ElectrodeGroup{3} = [2 3 4 5], electrode 3
%                    is a tetrode for channels 2 3 4 and 5.
%
% can use xpathReader:
%cmdstr =['xpathReader ' FileBase '.xml "//anatomicalDescription/channelGroups/group/channel[@skip=0]" > nsch '];
%system(cmdstr);
% or the xml toolbox

% open file
fp = fopen(FileName, 'r');
Par.FileName = FileName;

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
fnbase = FileName(1:strfind(FileName,'.par')-1);

if FileExists([fnbase '.eeg1.par'])
   
    EegPar=LoadEegPar(fnbase);
    for el=1:Par.nElecGps
        for eegel=1:EegPar.nElec
            if ~isempty(intersect(Par.ElecGp{el},EegPar.ElecChannels{eegel}))
                   Par.ElecLoc{el} = EegPar.ElecLoc{eegel};
            end
        end
	end
elseif Par.SampleTime ==40 % Peter's files
    for el=1:Par.nElecGps
        Par.ElecLoc{el} = 'c';
    end
end
