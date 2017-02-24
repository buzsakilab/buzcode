% segment_calc_script

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
F1=[]; F2=[]; F3=[]; F4=[]; F0=[]; C1=[]; C2=[]; 
for i=1:nwaves, 
    i
%    data=wavread(wavnames{i}{1});
%    [Feat S t f]= acoustic_features_MB(data,movingwin,params);
    sname=[wavnames{i}{1} '.sp'];
    load(char(sname),'-mat');
    C=fft(log(S),[],2); cep=real(C); C1=[C1 cep(:,2)']; C2=[C2 cep(:,3)'];  
    x0=log(mean(S,2))/log(10)*10; x0=x0-max(x0); F0=[F0 x0'];
    x=log(Feat(:,1)); w=Feat(:,2)-x; x=x/log(10)*10; x=x-max(x);  
    F1=[F1 x'];
    F2=[F2 w'];
    F3=[F3 Feat(:,3)'];
    tmp=load(char(segnames{i}(1))); sz=size(tmp);
    seg=zeros(length(t),1); 
    if sz(1)>0, 
        for j=1:sz(1), ind=find((t>tmp(j,1)) & (t<tmp(j,2))); seg(ind)=1; end
    end; 
    F4=[F4 seg'];
end
tb=find(F4<1); ts=find(F4>0);  
F5=zeros(length(F0),1);
F5(find(F0>-20))=1.1;
F6=medfilt1(F5,5);