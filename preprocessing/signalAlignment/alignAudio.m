function [] = alignAudio()

audioFolder = dir('*audio*');
if isempty(audioFolder)
    error('no audio folder found');
elseif length(audioFolder)>1
    error('too many folders?');
end
cd(audioFolder(1).name)
audioFiles = dir('*WAV');

audio = [];
for i=1:length(audioFiles)
    [wav fs{i}]=audioread(audioFiles(i).name,'native');
    audio = [audio;wav];
end
cd ..

if exist('analoginin.dat')
    analogin = uint16(bz_LoadBinary('analoginin.dat','nChannels',2,'channels',1));
end
analogin(analogin<mean(analogin)+std(single(analogin)))=mean(analogin);

%% take the first and last 5 minutes of the recording

analoginFS = 20000;

% audio envelope filtered at 5k





end