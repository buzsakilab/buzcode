function [ CSD ] = bz_CSD (lfp, channels, win, spat_sm, temp_sm,doPlot)

% [ CSD ] = bz_CSD (lfp, twin, spat_sm, temp_sm,doPlot)
%
% INPUTS:
%       LFP         LFP structure in buzcode format:
%                   a structure with fields lfp.data, lfp.timestamps
%       channels    vector with channels to inlcude. If empty take all (default)
%       win         time interval to compute CSD. If empty take all (default)
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       doPlot      true/false. Default true.

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
addParameter(p,'doPlot',true,@islogical);
parse(p);

%% Compute CSD

lfp_frag = lfp.data(p.Results.win(1):p.Results.win(2),p.Results.channels);

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

if doPLot
   cmax = max(max(CSD)); 
   
   figure;
   contourf(lfp.timestamps,1:size(CSD,2),CSD',40,'LineColor','none');
   colormap jet; caxis([-cmax cmax]);
   set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title(CSD); 
   
end



end




