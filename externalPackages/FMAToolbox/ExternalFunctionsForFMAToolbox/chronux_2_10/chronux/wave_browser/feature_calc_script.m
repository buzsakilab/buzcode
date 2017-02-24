% feature_calc_script

tmp=dir;
nfiles=length(tmp);
names=cell(nfiles,1); for j=1:nfiles, names(j)=cellstr(tmp(j).name); end; 
wavnames=regexp(names,'.*wav$','match');
ind=~cellfun('isempty',wavnames);
wavnames=wavnames(ind);
segnames=regexp(names,'.*wav.seg.txt$','match');
ind=~cellfun('isempty',segnames);
segnames=segnames(ind);

nwaves=length(wavnames);

movingwin=[0.01 0.005];
params.tapers=[2 3]; params.pad=1; params.Fs=44100; params.fpass=[5000 20000];
for i=1:nwaves, 
    data=wavread(wavnames{i}{1});
    [Feat S t f]= acoustic_features_MB(data,movingwin,params);
    sname=[wavnames{i}{1} '.sp'];
    save(char(sname),'S','t','f','Feat','-mat');
end
