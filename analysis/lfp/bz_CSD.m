function [ CSD ] = bz_CSD (lfp, channels, win, spat_sm, temp_sm, doDetrend, plotCSD, plotLFP)

% [ CSD ] = bz_CSD (lfp, twin, spat_sm, temp_sm,doPlot)
%
% INPUTS:
%       LFP         LFP structure in buzcode format:
%                   a structure with fields lfp.data, lfp.timestamps
%       channels    vector with channels to inlcude. If empty take all (default)
%       win         time interval to compute CSD. If empty take all (default)
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       plotCSD      true/false. Default true.

% OUTPUT:
%       CSD         LFP structure in buzcode format
%                    a structure with fields lfp.data, lfp.timestamps

% Antonio FR, 7/18

%% Parse inputs

p = inputParser;
addParameter(p,'channels',1:size(lfp.data,2),@isvector);
addParameter(p,'win',[1 size(lfp.data,1)],@isnumeric);
addParameter(p,'spat_sm',11,@isnumeric);
addParameter(p,'temp_sm',11,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);

parse(p);

%% Compute CSD

lfp_frag = lfp.data(p.Results.win(1):p.Results.win(2),p.Results.channels);

% detrend
if p.Results.doDetrend
   lfp_frag = detrend(lfp_frag')';
end
    
% spatial smoothing
if p.Results.spat_sm > 0
   for ch = 1:size(lfp_frag,2) 
       lfp_frag(:,ch) = smooth(lfp_frag(:,ch),p.Results.spat_sm,'sgolay');
   end
end

% temporal smoothing
if p.Results.temp_sm > 0
   for t = 1:size(lfp_frag,1) 
       lfp_frag(t,:) = smooth(lfp_frag(t,:),p.Results.temp_sm,'sgolay');
   end
end

% calculate CSD 

CSD = diff(lfp_frag,2,2);

%% Plot

if plotCSD
   cmax = max(max(CSD)); 
   
   figure;
   contourf(lfp.timestamps,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
   colormap jet; caxis([-cmax cmax]);
   set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title(CSD); 
   
end

if plotLFP
    
    figure;
    subplot(1,2,1);
    contourf(lfp.timestamps,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD'); 
   
    subplot(1,2,2);
    plotLFPprofile(lfp_frag',scale)
    set(gca,'YDir','reverse','YTickLabel',[]);xlabel('time (s)');ylabel('channel');title('LFP');   
    
end


%%
function plotLFPprofile( lfp,scale,varargin )
%   plotLFPprofile( lfp,scale,varargin )
%   plot vertical profile of LFP from a single shank
%   lfp = Nchannels x Nsamples
%   scale = multiplying factor for plot LFP
%
%   Azahara and Antonio, 2015

    nargin = length(varargin);

    for i =1:nargin
    if strcmp(varargin{i},'DoDetrend')
           lfp = detrend(lfp);
    end
    end

    %figure;   
    for i=1:size(lfp,1)
        offset = 500*(i-1);
        sh_tmp = scale*lfp(i,:) + offset;
        plot(sh_tmp,'k','LineWidth',2); hold on;
        clear sh_tmp
    end
    
    for i =1:nargin
        if strcmp(varargin{i},'NoLabels')
           set(gca,'box','off');axis off;
        end
    end

end

