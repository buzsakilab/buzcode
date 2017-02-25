% function testscript(pname,direction,movingwin,segave,params,fscorr)
%
% This script runs a sequence of analysis steps using the test
% data contained in data. The data consists of a single tetrode
% recording from macaque area LIP during a memory saccade experiment
% a la Pesaran et al (2002). The data is already separated into
% spikes and LFPs. LFPs are contained in variable dlfp, spikes from two
% neurons are in in a struct array dsp, and event information is in
% the following set of variables:
%
% trialtimes - start times of trials
% fixon - fixation light comes on
% fixacq - fixation acquired
% targon - target light on
% targoff - target light off
% fixoff - fixation off
% saccade - saccade
%
% Note that spikes and event times are in seconds and the sampling
% frequency for the LFP in this experiment was 1kHz.
%
% the script takes the following input argument - 
% pname - path name on your computer where the data file LIPdata is stored.
% direction - target direction to be analysed (0-7)
%
% The remaining parameters control various computations and are discussed
% in chronux.m - type Help chronux.m for more information.
% 
% if nargin < 4;
%     error('Need 6 input parameters - see help');
% end;
% if nargin < 5 | isempty(params);
%    [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
% end;
% if nargin < 6 | isempty(fscorr);
%     fscorr=1;
% end;
pname='data';
params.Fs=1000; % sampling frequency
params.fpass=[10 100]; % band of frequencies to be kept
params. tapers=[3 5]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[2 0.05];
params.trialave=1;
movingwin=[0.5 0.05];
segave=1;
direction=5;

wintrig=[5*movingwin(1) 5*movingwin(1)];
winseg=2*movingwin(1);
%
% Load data
%
eval(['load LIPdata.mat']);
% %
% % Create rearranged data blocks for further analysis: We are going to
% % extract segments of data centered on the target off times from the first channel of LFP data and the from one of the two spike trains
% % 
% %
% indx1=find(targets==5);indx2=find(targets==1); % trials to preferred and antipreferred direction
% E1=targon(indx1); E2=targon(indx2); % target on times trials to preferred adn anti-preferred directions
% dlfp1=createdatamatc(dlfp(:,1),E1,Fs,wintrig);dlfp2=createdatamatc(dlfp(:,1),E2,Fs,wintrig); % extract event triggered segments of the first LFP channel
% dsp1=createdatamatpt(dsp(1),E1,wintrig); dsp2=createdatamatpt(dsp(1),E2,wintrig); % the same for one of the spike trains

% compute spectrum of the first few seconds of LFP channels 1-2
NT=round(params.Fs*10*movingwin(1));
data=dlfp(1:NT,:); data1=data(:,1:2);
[S,f,Serr]=mtspectrumc(data1,params);
figure;
plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:))); xlabel('Frequency Hz'); ylabel('Spectrum');
%%% pause

% compute derivative of the spectrum for the same data
phi=[0 pi/2];
[dS,f]=mtdspectrumc(data1,phi,params);
figure;
plot(f,dS(1,:),f,dS(2,:)); xlabel('frequency Hz'); ylabel('Derivatives'); legend('Time','Frequency');
%%% pause

% compute coherency between  channels 1-2 and  3-4
data1=data(:,1:2);data2=data(:,3:4);
[C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(data1,data2,params); 
figure;plot(f,C,f,Cerr(1,:),f,Cerr(2,:));xlabel('frequency'); ylabel('Coherency');
%%% pause

% coherency matrix of data1
[C,phi,S12,f,confC,phierr,Cerr]=cohmatrixc(data1,params);

% compute spectrogram for 1-2
[S,t,f,Serr]=mtspecgramc(data1,movingwin,params);
figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar
%%% pause

% compute time-frequency derivative of the spectrogram for 1-2
phi=[0 pi/2];
[dS,t,f]=mtdspecgramc(data1,movingwin,phi,params);
% pause
% figure;subplot(211);imagesc(t,f,squeeze(dS(1,:,:))'); axis xy; colorbar;
% subplot(212);imagesc(t,f,squeeze(dS(2,:,:))'); axis xy; colorbar;
% %%% pause

% compute coherogram between 1-2 and 3-4
NT=round(movingwin(1)*Fs);
[C,phi,S12,S1,S2,t,f,confC,phierr,Cerr]=cohgramc(data1,data2,movingwin,params);
figure;imagesc(t,f,C'); axis xy; colorbar;
%%% pause

% compute segmented spectrum of 1 
NT=10*round(winseg*Fs);
data1=dlfp(1:NT,1);
[S,f,varS,C,Serr]=mtspectrumsegc(data1,winseg,params,segave);
figure; subplot(211);plot(f,10*log(S));
imagesc(f,f,C); axis xy;colorbar;
%%% pause

% compute segmented coherency between 1 and 2
NT=10*round(winseg*Fs);
data1=dlfp(1:NT,1); data2=dlfp(1:NT,2);
[C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencysegc(data1,data2,winseg,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1)); subplot(313);plot(f,10*log10(S2))
%%% pause

% compute spectrum of channel 1 triggered to events E
E1=targon(find(targets==direction)); 
data1=dlfp(:,1); 
[S,f,Serr]=mtspectrumtrigc(data1,E1,wintrig,params);
figure;plot(f,10*log10(S)'); axis xy; colorbar;
%%% pause

% compute spectrogram of channel 1 triggered to events E
E1=targon(find(targets==direction)); 
data1=dlfp(:,1); 
[S,t,f,Serr]=mtspecgramtrigc(data1,E1,wintrig,movingwin,params);
figure; imagesc(t,f,10*log10(S)'); axis xy; colorbar;
%%% pause

%
% Analysis - point process stored as times
%

% dsp contains 2 channels of spikes 

data=extractdatapt(dsp,[20 30]); % extract spikes occurring between 20 and 30 seconds and compute their spectrum
[S,f,R,Serr]=mtspectrumpt(data,params);
figure; plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]);
%%% pause

%
% Compute the derivative of the spectrum
%
phi=[0 pi/2];
[dS,f]=mtdspectrumpt(data,phi,params);
figure; plot(f,dS);
%%% pause

%
% Compute the derivative of the time-frequency spectrum
%
data=extractdatapt(dsp,[20 30]);
data1=data(1); data2=data(2);fscorr=[];t=[];
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypt(data1,data2,params);
figure; plot(f,C);
%%% pause

%
% Compute event triggered average spectrum for one of the directions
%
E1=targon(find(targets==direction));
data=dsp(1); 
[S,f,R,Serr]=mtspectrumtrigpt(data,E1,wintrig,params);
figure;plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:))); line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]); 
%%% pause

%
% Compute the matrix of coherencies
%
data=extractdatapt(dsp,[20 30]);
[C,phi,S12,f,zerosp,confC,phierr,Cerr]=cohmatrixpt(data,params,fscorr);

%
% Event triggered spectrogram - first way way
%
data=createdatamatpt(dsp(1),E1,wintrig);
[S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params);
figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar;
%%% pause

%
% Derivative of the time-frequency spectrum
%
data=createdatamatpt(dsp(1),E1,wintrig);
phi=[0 pi/2];
[dS,t,f]=mtdspecgrampt(data,movingwin,phi,params);
figure; subplot(211); imagesc(t,f,squeeze(dS(1,:,:))'); axis xy; colorbar;
subplot(212); imagesc(t,f,squeeze(dS(2,:,:))'); axis xy; colorbar;

%
% Coherogram between the two spike trains
%
data1=createdatamatpt(dsp(1),E1,wintrig);
data2=createdatamatpt(dsp(2),E1,wintrig);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgrampt(data1,data2,movingwin,params,fscorr);
figure;imagesc(t,f,C');axis xy; colorbar
% %%% pause

%
% Event Triggered spectrogram another way
%
data=dsp(1);
[S,t,f,R,Serr]=mtspecgramtrigpt(data,E1,wintrig,movingwin,params);
imagesc(t,f,10*log10(S)'); axis xy; colorbar
%
% Segmented spectrum
%
data=extractdatapt(dsp,[20 30]);
data=data(1);
[S,f,R,varS]=mtspectrumsegpt(data,winseg,params);
plot(f,10*log10(S)); line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]);
%
% Segmented coherency
%
data=extractdatapt(dsp,[20 30]);
data1=data(1);data2=data(2);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpt(data1,data2,winseg,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))
%
% Analysis - hybrid: one continous and one point process stored as times
%
offset=1;
data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpt(data1,data2,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))


data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegcpt(data1,data2,winseg,params,segave,fscorr);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))


data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgramcpt(data1,data2,movingwin,params,fscorr);
figure; subplot(311); imagesc(t,f,C');axis xy; colorbar; subplot(312);imagesc(t,f,10*log10(S1)');axis xy; colorbar; subplot(313); imagesc(t,f,10*log10(S2)');axis xy; colorbar


%  Analysis: Binned spike counts

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds

data=dN;
[S,f,R,Serr]=mtspectrumpb(data,params);
plot(f,10*log10(S),f,10*log10(Serr(1,:)), f,10*log10(Serr(2,:))); %line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]);

data=dN;
phi=[0 pi/2];
[dS,f]=mtdspectrumpb(data,phi,params);
figure; plot(f,dS); 

data=dN;
data1=data(:,1); data2=data(:,2);fscorr=[];t=[];
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypb(data1,data2,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1)); subplot(313); plot(f,10*log10(S2));

E=targon(find(targets==direction)); E=E(find(E>20 & E<450)); [dN,t]=binspikes(dsp,params.Fs,[20 500]);data=dN(:,1);
[S,f,R,Serr]=mtspectrumtrigpb(data,E,wintrig,params);
figure;plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
[C,phi,S12,f,zerosp,confC,phierr,Cerr]=cohmatrixpb(data,params,fscorr);

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
[S,t,f,R,Serr]=mtspecgrampb(data,movingwin,params);
figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar

clear R Serr Cerr S C phierr dN dS data1 data2

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
phi=[0 pi/2];
[dS,t,f]=mtdspecgrampb(data,movingwin,phi,params);

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
data1=data(:,1); data2=data(:,2);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgrampb(data1,data2,movingwin,params,fscorr);

[dN,t]=binspikes(dsp,params.Fs); % extract spikes
dN=dN(:,1);
E=E(1:6);
data=dN;
[S,t,f,R,Serr]=mtspecgramtrigpb(data,E,wintrig,movingwin,params);

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
data=data(:,1);
[S,f,R,varS]=mtspectrumsegpb(data,winseg,params,segave,fscorr);
figure; plot(f,10*log10(S)); line(get(gca,'xlim'),10*log10([R R])); 

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
data1=data(:,1);data2=data(:,2);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpb(data1,data2,winseg,params);

%
% Analysis - hybrid: one continous and one point process stored as counts
%
data1=dlfp(20000:30000,:);data1=data1(:,1:2); [dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data2=dN; data2=data2(1:end,:);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpb(data1,data2,params);

data1=data1(:,1); data2=data2(:,1);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegcpb(data1,data2,winseg,params,segave,fscorr);


data1=dlfp(20000:30000,:); data1=data1(:,1:2);
[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data2=dN; data2=data2(1:end,:);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgramcpb(data1,data2,movingwin,params,fscorr);

