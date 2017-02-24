function Par = LoadPar(FileName, varargin)
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
[SpecInfo] = DefaultArgs(varargin,{1});

if ~isempty(strfind(FileName,'.par'))
    FileBase = FileName(1:strfind(FileName,'.par')-1);
elseif ~isempty(strfind(FileName,'.xml'))
    FileBase = FileName(1:strfind(FileName,'.xml')-1);
else 
    FileBase = FileName;
end


if FileExists([FileBase '.xml']) %& ~isempty(strfind(FileName,'.xml'))
    Par = LoadXml(FileBase);
elseif FileExists([FileBase '.par'])


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