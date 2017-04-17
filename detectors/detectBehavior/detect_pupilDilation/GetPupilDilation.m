function [ pupildilation ] = GetPupilDilation( baseName,basePath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Note: you may need gstreamer: https://gstreamer.freedesktop.org/download/

%baseName = 'VideoFrame_Test_50Hz_20000_170328_173154';
%recordingsfolder = 'C:\Users\rudylabadmin\Desktop\layers\Recordings';
%recordingsfolder = '/mnt/proraidDL/Database/WMProbeData';

%%

if nargin==0
    [basePath,baseName] = fileparts(pwd);
end
%%


addpath(fullfile(basePath,baseName));

vidName = fullfile(basePath,baseName,[baseName,'.avi']);
abfname = fullfile(basePath,baseName,[baseName,'.abf']);
analogName = fullfile(basePath,baseName,['analogin.dat']);

savefile = fullfile(basePath,baseName,[baseName,'.pupildiameter.behavior.mat']);
savevid = fullfile(basePath,baseName,'DetectionFigures',[baseName,'.pupilvid.avi']);

SAVEVID = true;
savevidfr = 10;
if SAVEVID
    pupdiamVid = VideoWriter(savevid);
    pupdiamVid.FrameRate = 1./(0.015.*savevidfr);
    open(pupdiamVid);
end
%%
%Load in the video (has to be in not .avi for now...)
pupilvidobj = VideoReader(vidName);
NumberOfFrames = pupilvidobj.NumberOfFrames;

%Use the first frame to get image dimensions
vidframe=read(pupilvidobj,10);
imagesize=size(vidframe);

%Initiate the pupil area vector, center coordinates
puparea = nan(NumberOfFrames,1);
pupcoords = nan(NumberOfFrames,2);

%Loop through the frames and get pupil location/area
pupfig = figure;
 colormap('gray')
for ff = 1:NumberOfFrames 
%ff = 1;
if mod(ff,200)==0
display([num2str(ff),' of ' num2str(NumberOfFrames),' Frames'])
end

     vidframe=read(pupilvidobj,ff);
     
     %Convert to greyscale
    if size(vidframe,3)==3
        vidframe = rgb2gray(vidframe);
    end
    vidframe_orig = vidframe;


    
    %Define the mask: trace the eye
    if ff==1
        
        subplot(2,2,1)
        imagesc(vidframe_orig);
        title('Trace the Eye')
        h = imfreehand;
        mask = ~createMask(h);
        %keyboard
%         mask = false(size(vidframe));
%         ymax = 70; ymin = 40; xmax = 60; xmin = 35;
%         mask(xmax:end,:) = true;
%         mask(1:xmin,:) = true;
%         mask(:,1:ymin) = true;
%         mask(:,ymax:end) = true;
%         
%         %
%         maskval =  max(vidframe(:));
        maskval = 255;
    end
    vidframe(mask) = maskval;
    

    %Define the pupil size and intensity threshold - trace the pupil
    if ff==1
        
            %Show the Masked frame
            subplot(2,2,3)
         imagesc(vidframe)
        title('Trace the Pupil')
        
        h = imfreehand;
        pupilmask = createMask(h);
        
        %Get the pixel values for pupil, eye, and not-pupil eye
        pupilpixels = double(vidframe(pupilmask))./255; 
        irispixels  = double(vidframe(~mask & ~ pupilmask))./255; 
        eyepixels = double(vidframe(~mask))./255; 
        
        intensitybins = linspace(0,1,40);
        pupilhist = hist(pupilpixels,intensitybins);
        irishist = hist(irispixels,intensitybins);
        eyehist = hist(eyepixels,intensitybins);
        
       % keyboard
        pupilsizethresh = 40; %Pupil must be larger than this many pixels.
        %Pupil must be darker than this intensity: 2.5std above mean pupil
        intensitythresh = mean(pupilpixels)+1.*std(pupilpixels); 
        %%
        x = 1;
        while ~isempty(x)
        tentativepupil=2.*~im2bw(vidframe,intensitythresh);
        tentativenotpupil=im2bw(vidframe,intensitythresh);
        
        tentativemap = ~mask+tentativepupil+tentativenotpupil;
        subplot(2,2,3)
        imagesc(tentativemap)
        
        subplot(2,2,2)
        bar(intensitybins,eyehist)
        hold on
        plot(intensitybins,irishist,'g','linewidth',2)
        plot(intensitybins,pupilhist,'r','linewidth',2)
        plot([1 1].*intensitythresh,get(gca,'ylim'),'r--')
        xlim([0 1])
        legend('Whole Eye','Not Pupil','Pupil','location','northwest')
        
       %Show eye with over/under pixels. allow user to select threshold and show new over/under 
       title({'Click to adjust threshold', 'or press RETURN for current threshold'})
        [x,~] = ginput(1);
        hold off
        if isempty(x); break; end
        intensitythresh = x; 
        end
          %%

    

        %%
        set(pupfig,'visible','off')   
    end
     
    %Get the pupil - all pixels below intensity threshold (black)
    pupil=~im2bw(vidframe,intensitythresh);

    if SAVEVID && mod(ff,savevidfr)==0; 
        subplot(4,4,9);imagesc(pupil); end
    %Show the thresholded image

    
    %Remove small objects and holes (i.e. dark/bright spots)
    pupil=bwmorph(pupil,'close');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,10);imagesc(pupil); end
    
    pupil=bwmorph(pupil,'open');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,11);imagesc(pupil); end
    
    pupil=bwareaopen(pupil,pupilsizethresh);
    pupil=imfill(pupil,'holes');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,12);imagesc(pupil); end

    % Tagged objects in BW image
    L=bwlabel(pupil);
    % Get areas and tracking rectangle
    out_a=regionprops(L);
    % Count the number of objects
    N=size(out_a,1);
    
    if N < 1 || isempty(out_a) % Returns if no object in the image
        continue
    end

    % Select larger area
    areas=[out_a.Area];
    [area_max pam]=max(areas);
    puparea(ff) = area_max;
    
    centro=round(out_a(pam).Centroid);
    X=centro(1);
    Y=centro(2);
  
    if SAVEVID && mod(ff,savevidfr)==0;  
    subplot(2,2,1)
    imagesc(vidframe_orig);
    hold on
    rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
        'Curvature', [1,1],'LineWidth',1)
    plot(X,Y,'g+')
    hold off
    end
    
    pupcoords(ff,1) = X; pupcoords(ff,2) = Y;
    
   if SAVEVID && mod(ff,savevidfr)==0;  
   subplot(4,1,4)
   plot(1:ff,puparea(1:ff),'k')
   end
     
   if SAVEVID && mod(ff,savevidfr)==0;
       imgFrame = getframe(gcf);
       writeVideo(pupdiamVid,imgFrame.cdata);
   end
end

%Close the video object
if SAVEVID; close(pupdiamVid); end


%0-1 Normalize, smooth the pupil area trace 
%(make window in units of seconds or figure out the best # frames)
smoothwin = 10;  %frames
puparea_pxl = puparea;
puparea = smooth(puparea,smoothwin);
puparea = (puparea-min(puparea))./(max(puparea)-min(puparea));

%% Load the analogin for the timestamps

timepulses = readmulti(analogName,1);
%LoadBinary('uint16')

sf_pulse = 1./20000; %Sampling Frequency of the .abf file
t_pulse = [1:length(timepulses)]'.*sf_pulse;

pulsethreshold =1e4;  %Adjust this later to set based on input.
pulseonsets = find(diff(timepulses<pulsethreshold)==1);
pulset = t_pulse(pulseonsets);

minpulsedur = 0.003; %Remove double/noise crossings

shortpulses=diff(pulset)<(minpulsedur);
pulset(shortpulses) = [];

interpulse = diff(pulset);

% longpulses=diff(pulset)>(2.*expectedpulserate);
% pulset(longpulses) = [];

%hist(diff(pulset))

%% Clampex Pulses
% timepulses = abfload(abfname);
% timepulses = timepulses(:,1);
% 
% sf_abf = 1./20000; %Sampling Frequency of the .abf file
% t_pulse = [1:length(timepulses)]'.*sf_abf;
% 
% pulsethreshold =0.5;  %Adjust this later to set based on input.
% pulseonsets = find(diff(timepulses<pulsethreshold)==1);
% pulset = t_pulse(pulseonsets);


%pulset(1) = []; %remove the first trigger... make this more rigorous later 

%% Check that the number of identified pulses = the number of frames
if length(pulset)~=length(puparea); 
    display('WARNING: NUMBER OF FRAMES DON"T MATCH!!!    v sad.');%keyboard; 
end

%%
figure
subplot(2,1,1)
plot(t_pulse,timepulses,'k')
hold on
plot(pulset,zeros(size(pulset)),'r+')

subplot(2,1,2)
hist(interpulse)

%%
trange = pulset([1 end]);
numframes = length(puparea);
t_interp = linspace(trange(1),trange(2),numframes)';

%%





%% Behavior struct

pupildilation.t_pulse = pulset;
pupildilation.t_interp = t_interp;
pupildilation.puparea = puparea;
pupildilation.puparea_pxl = puparea_pxl;
pupildilation.pupcoords = pupcoords;
pupildilation.detectorname = 'GetPupilDilation';
pupildilation.detectiondate = today('datetime');
pupildilation.detectorparms.mask = mask;
pupildilation.detectorparms.intensitythresh = intensitythresh;


save(savefile,'pupildilation')


end

