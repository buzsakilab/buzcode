function [ pupildilation ] = GetPupilDilation(basePath)
%This is a detector for extracting pupil diameter from an IR video of the
%eye.
%
%FILE REQUIREMENTS
%   This detector assumes the presence of a few files:
%       basePath/baseName.avi   pupil video
%       basePath/analogin.dat   analog signal of frame pulses from intan 
%   where basePath is a folder of the form: 
%       whateverPath/baseName/
%
%INPUT
%   basePath    (optional) folder in which the recording files live.
%               of the form whateverPath/baseName/
%               Default: current working directory
%
%OUTPUT: saved as basePath/baseName.pupildiameter.behavior.mat
%   pupildilation
%       .t_pulse        time of the pulses recieved by intan
%       .t_interp       frame times interpolated between 1st/last intan pulse
%       .puparea        [0-1] normalized size of the pupil
%       .puparea_pxl    pupil size in pixels
%       .pupcoords      [x y] coordinates of the pupil center
%       .detectorparms  structure of detection paramters
%
%
%Note: you may need gstreamer: https://gstreamer.freedesktop.org/download/
%
%
%USING THE DETECTOR
%This detector some amount of user input to determine best intensity
%thresholds for the pupil and eyelid for blink and unstable frames
%detection. Upon running the function, user will be promped with a series
%of steps, described below with tips for best performance.
%
%Step 0: Make sure your files conform to the buzcode guidelines above
%Step 1: Trace the eye.
%-Trace the eyelid. This serves two purposes. First, excluding dark shadows
%in the corner of the eye that may be mistaken for a pupil. Second, this
%eye outline serves as a bseline for detection of blinking and other
%instabilities. If the eye presented in the first frame is not ideal,
%please press or hold "ESC" until an ideal eye presents itself, as this 
%will be used as the template throughout the rest of the recording.
%Step 2: Trace the pupil.
%-This doesn't have to be perfect, but should avoid including any non-pupil
%pixels.
%Step 3: Adjust pupil/nonpupil threshold
%-Place the threshold in the histogram to get the best separation between 
%the two histograms. Looking for most pupil in the grey image with least
%amount of shadow contamination.
%
%
%DLevenstein 2017
%With many inputs from WMunoz to improve detection quality.
%%

if ~exist('basePath','var')
    basePath = pwd;
end
[baseFolder,baseName] = fileparts(basePath);
%%
WHISKCONTAMINATION = true;

%%

%addpath(fullfile(baseFolder,baseName));

vidName = fullfile(basePath,[baseName,'.avi']);
analogName = fullfile(basePath,['analogin.dat']);

savefile = fullfile(basePath,[baseName,'.pupildiameter.behavior.mat']);
figfolder = fullfile(basePath,'DetectionFigures');
savevid = fullfile(figfolder,[baseName,'.pupilvid.avi']);

SAVEVID = true;
savevidfr = 10;
if SAVEVID
    if ~exist(figfolder,'dir')
        mkdir(figfolder)
    end
    pupdiamVid = VideoWriter(savevid);
    pupdiamVid.FrameRate = 1./(0.015.*savevidfr);
    open(pupdiamVid);
end
%%
%Load in the video
pupilvidobj = VideoReader(vidName);
NumberOfFrames = pupilvidobj.NumberOfFrames;

%Use the first frame to get image dimensions
vidframe=read(pupilvidobj,10);
imagesize=size(vidframe);

%Initiate the pupil area vector, center coordinates
puparea = nan(NumberOfFrames,1);
pupcoords = nan(NumberOfFrames,2);
unstableframes = [];

meanvid = zeros(imagesize(1:2));

%Loop through the frames and get pupil location/area
pupfig = figure;
 colormap('gray')
for ff = 1:NumberOfFrames 
if mod(ff,200)==0
display([num2str(ff),' of ' num2str(NumberOfFrames),' Frames'])
end

    %Read the video and convert to grayscale
	vidframe=read(pupilvidobj,ff); 
    if size(vidframe,3)==3
        vidframe = rgb2gray(vidframe);
    end
    vidframe_orig = vidframe; %Hold on to the original image for later
    
     meanvid = meanvid+single(vidframe_orig)./NumberOfFrames;
    
    %USER: Define the full eye mask: trace the eye
    if ~exist('eyemask','var')        
        subplot(2,2,1)
            imagesc(vidframe_orig);
            title({'Trace the Eye','Hit "ESC" for next frame'})
            set(gca,'ytick',[]);set(gca,'xtick',[]);
        h_eye = imfreehand;
        if isempty(h_eye); continue; end %Go to next frame if user aborts
        eyeline = getPosition(h_eye);
        eyemask = createMask(h_eye);
        noneyeval = 255;
    end
    vidframe(~eyemask) = noneyeval; %Put the noneye pixels to saturation 
    
    %Define the pupil size and intensity threshold - trace the pupil
    if ~exist('pupilmask','var') 
        %Show the Masked frame
        subplot(2,2,3)
            imagesc(vidframe)
            title('Trace the Pupil')
        %Trace the pupil
        h = imfreehand;
        pupilmask = createMask(h);
        
        %Get the pixel values for pupil, eye, and not-pupil eye
        pupilpixels = double(vidframe(pupilmask))./255; 
        irispixels  = double(vidframe(eyemask & ~pupilmask))./255; 
        eyepixels = double(vidframe(eyemask))./255; 
        
        %Blink threshold
        unstablemeanthresh = 0.6; %std
        unstablestdthresh = 0.15;
        eyerange = mean(eyepixels)+unstablemeanthresh.*std(eyepixels).*[-1 1];
        blinkstdthresh = std(eyepixels)+unstablestdthresh.*[-1 1];
        
        intensitybins = linspace(0,1,40);
        pupilhist = hist(pupilpixels,intensitybins);
        irishist = hist(irispixels,intensitybins);
        eyehist = hist(eyepixels,intensitybins);
        
        pupilsizethresh = 15; %Pupil must be larger than this many pixels
        %Starting intensity threshold: 2.5std above mean pupil
        intensitythresh = mean(pupilpixels)+1.*std(pupilpixels); 
        x = 1;
        while ~isempty(x)
            tentativepupil=2.*~im2bw(vidframe,intensitythresh);
            tentativenotpupil=im2bw(vidframe,intensitythresh);
            tentativemap = eyemask+tentativepupil+tentativenotpupil;
            subplot(2,2,3)
                imagesc(tentativemap)
            subplot(2,2,2)
                bar(intensitybins,eyehist)
                hold on
                plot(intensitybins,irishist,'g','linewidth',2)
                plot(intensitybins,pupilhist,'r','linewidth',2)
                plot([1 1].*intensitythresh,get(gca,'ylim'),'r--')
                xlim([0 1])
                legend('Whole Eye','Not Pupil','Pupil',...
                    'location','northwest')
        
            %Show eye with over/under pixels. 
            %allow user to select threshold and show new over/under 
            title({'Click to adjust threshold',...
               'or press RETURN for current threshold'})
            [x,~] = ginput(1);
            hold off
            if isempty(x); break; end
            intensitythresh = x; 
        end
        
        %%
       % set(pupfig,'visible','off')   
    end
     
    
    %Detect blinking very large/small avgerage/std eye pixel intensity
    unstablewindow = 10; %window of frames around detected unstable frames to denote as unstable
  %  integratetime =  %try just using a window in the past?
    UNSTABLE = 0;
    meaneyepixel(ff) = mean(double(vidframe(eyemask))./255);
     stdeyepixel(ff) = std(double(vidframe(eyemask))./255);
    %display(['Mean: ',num2str(meaneyepixel(ff)),'.   Range: ',num2str(eyerange)])
    if meaneyepixel(ff) < eyerange(1) || meaneyepixel(ff) > eyerange(2) ||...
            stdeyepixel(ff) < blinkstdthresh(1) || stdeyepixel(ff) > blinkstdthresh(2)
        UNSTABLE = 1;
        unstableframes = unique([unstableframes,[ff-unstablewindow:ff+unstablewindow]]);
        unstableframes(unstableframes<1 | unstableframes>NumberOfFrames)=[];
        puparea(unstableframes)=nan;
    end
    
    
    %Get the pupil - all pixels below intensity threshold (black)
    pupil=~im2bw(vidframe,intensitythresh);

    if SAVEVID && mod(ff,savevidfr)==0; 
        subplot(4,4,9);imagesc(pupil);
        set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
    %Show the thresholded image

    
    %Remove small objects and holes (i.e. dark/bright spots)
   % pupil=bwmorph(pupil,'close');
    pupil=imclose(pupil,strel('disk',4)); %Number of pixels to try to merge
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,10);imagesc(pupil);
        set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
    
    pupil=bwmorph(pupil,'open');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,11);imagesc(pupil);
        set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
    
    pupil=bwareaopen(pupil,pupilsizethresh);
    pupil=imfill(pupil,'holes');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,12);imagesc(pupil);
        set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
    
     
    % Tagged objects in BW image
    L=bwlabel(pupil);
    % Get areas and tracking rectangle
    out_a=regionprops(L);
    % Count the number of objects
    N_objects=size(out_a,1);
    
    if N_objects < 1 || isempty(out_a) || UNSTABLE || ismember(ff,unstableframes)% Returns if no object in the image
        %Save nans for unstable framess
        X=nan;
        Y=nan;
        puparea(ff) = nan;
        %Stuff so plot doesn't crash
        pam = 1;
        out_a(pam).BoundingBox = [0 0 0 0];
    else

        % Select larger area
        areas=[out_a.Area];
        [area_max pam]=max(areas);
        puparea(ff) = area_max;

        centro=round(out_a(pam).Centroid);
        X=centro(1);
        Y=centro(2);
    end
  
    if SAVEVID && mod(ff,savevidfr)==0;  
    subplot(2,2,1)
    imagesc(vidframe_orig);
    %caxis([0 255])
    hold on
    plot(eyeline(:,1),eyeline(:,2),'b:','linewidth',0.5)
    rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
        'Curvature', [1,1],'LineWidth',0.5)
    plot(X,Y,'g+','markersize',3)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    hold off
    end
    
    pupcoords(ff,1) = X; pupcoords(ff,2) = Y;
    
   if SAVEVID && mod(ff,savevidfr)==0;  
   subplot(4,1,4)
       yrange = [0 max(puparea)];
       windur = 3000; %frames
       earlypoint = max(1,ff-windur);
       plot(earlypoint:ff,(puparea(earlypoint:ff)-yrange(1))./diff(yrange),'k')
       hold on
       plot(ff,puparea(ff),'ro')
       plot(unstableframes,zeros(size(unstableframes)),'r.','markersize',10)
       hold off
       xlim([earlypoint ff])
       %if yrange(1)~=yrange(2); ylim([min(puparea) max(puparea)]); end
       ylim([0 1])
   end
     
   if SAVEVID && mod(ff,savevidfr)==0;
       imgFrame = getframe(gcf);
       writeVideo(pupdiamVid,imgFrame.cdata);
   end
end

%Close the video object
if SAVEVID; close(pupdiamVid); end


%% Load the analogin for the timestamps

timepulses = LoadBinary(analogName,'nChannels',1,'precision','uint16');
%LoadBinary('uint16')

sf_pulse = 1./20000; %Sampling Frequency of the .abf file
t_pulse = [1:length(timepulses)]'.*sf_pulse;

pulsethreshold =1e4;  %Adjust this later to set based on input.
%Find where pulse channel crosses in the upward direction
pulseonsets = find(diff(timepulses>pulsethreshold)==1); 
pulset = t_pulse(pulseonsets);


minpulsedur = 0.003; %Remove double/noise crossings
shortpulses=diff(pulset)<(minpulsedur);
pulset(shortpulses) = [];

interpulse = diff(pulset);
sf_eff = 1./mean(interpulse);

%Check that frame duration is constant up to tolerance (no skipped frames)
tol = 0.001;

if range(interpulse)>tol
    warning('Frame rate is not constant...')
end

% longpulses=diff(pulset)>(2.*expectedpulserate);
% pulset(longpulses) = [];

%hist(diff(pulset))


%% Check that the number of identified pulses = the number of frames
if length(pulset)~=length(puparea); 
    display('WARNING: NUMBER OF FRAMES DON"T MATCH!!!    v sad.');%keyboard; 
end



%% Estimate and check sampling frequency of the camera



%% Smooth and normalize the pupil diameter trace

%0-1 Normalize, smooth the pupil area trace 
%(make window in units of seconds or figure out the best # frames)
smoothwin_s = 0.5;  %s
%smoothwin = 10;  %frames
smoothwin_frames = round(smoothwin_s.*sf_eff); %Calculate window in frames 
puparea_raw = puparea; %In units of pixels
puparea = smooth(puparea,smoothwin_frames,'moving'); %,'rloess'); %try rloess method... needs percentage of total points span
%Long unstable epochs should be kept nan.....
maxarea = max(puparea);
puparea = puparea./maxarea;

%% Linearly interpolate bad frames of duration less than some window...



%% Align the timestamps
trange = pulset([1 end]);
numframes = length(puparea);
t_interp = linspace(trange(1),trange(2),numframes)';

%%
timestamps = pulset(1:NumberOfFrames);
data = puparea;

%% Detection Quality Control Figure
figure

    subplot(4,1,3)
        plot(timestamps,puparea,'k')
        hold on
        plot(timestamps(unstableframes),zeros(size(unstableframes)),'r.','markersize',10)
        hold on
        set(gca,'xticklabel',[])
        ylabel({'Pupil Diameter','(max^-^1)'})
        xlim(t_pulse([1 end]))
    subplot(4,1,4)
        plot(t_pulse,timepulses,'k')
        hold on
        plot(pulset,zeros(size(pulset)),'r+')
        xlabel('t (s)');ylabel('Pulses')
        set(gca,'ytick',[])
        xlim(t_pulse([1 end]))
        
        


    subplot(6,2,2)
        plot(t_pulse,timepulses,'k')
        hold on
        plot(pulset,zeros(size(pulset)),'r+')
        xlim(pulset(1)+[-0.05 0.1])
        xlabel('t (s)');ylabel('Pulses')
        set(gca,'ytick',[])
        title('First Pulse')
        
    subplot(6,2,4)
        hist(interpulse)
        xlabel('Frame Duration')
    subplot(6,2,6)
        hist(puparea,linspace(0,1,20))
        xlim([0 1])
        xlabel('Pupil Diameter')
        
        
    subplot(2,2,1)
    imagesc(meanvid)
    hold on
        plot(pupcoords(:,1),pupcoords(:,2),'.')
        title(baseName)
        
NiceSave('PupilDetection',figfolder,baseName)



%% Behavior struct for saving

pupildilation.samplingRate = sf_eff;
pupildilation.timestamps = timestamps;
pupildilation.data = data;

pupildilation.t_pulse = pulset; 
pupildilation.puparea_pxl = puparea_raw;

pupildilation.eyepixelmean = meaneyepixel;
pupildilation.eyepixelstd = stdeyepixel;
pupildilation.unstableframes = unstableframes;

pupildilation.pupilxy = pupcoords;

pupildilation.detectorname = 'GetPupilDilation';
pupildilation.detectiondate = today('datetime');
pupildilation.detectorparms.mask = eyemask;
pupildilation.detectorparms.intensitythresh = intensitythresh;
pupildilation.detectorparms.unstablewindow = unstablewindow;
pupildilation.detectorparms.blinkstdthresh = blinkstdthresh;
pupildilation.detectorparms.eyerange = eyerange;


save(savefile,'pupildilation')


end

