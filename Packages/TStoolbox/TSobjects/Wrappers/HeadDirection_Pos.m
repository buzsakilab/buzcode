function [ang,GoodRanges,ep] = HeadDirection_Wrapper(fbasename,varargin)


% loads head-direction from a position file file (ending in .whl)
%
% USAGE
%     [ang,GoodRanges] = HeadDirectionWhl(fbasename)
%     
% INPUT:
%     fbasename: session file basename
%	
% OUTPUT
%     ang: a tsd object of angular values
%     GoodRanges: a intervalSet object where LEDs were successfully detected


% Adrien Peyrache 2011

if ~isempty(varargin)
    w = varargin{1};
    if ~isa(w,'struct')
        error('argument must be a struct')
    end
    whl = w.whl;
    t = w.t;
    GoodRanges = w.GR;
else
    [whl,t,GoodRanges] = LoadPosition_Pos(fbasename);
end

[folder dumy dumy]  = fileparts(fbasename);
if ~isempty(folder)
    folder =  [folder filesep];
end

if exist([folder 'Analysis/GeneralInfo.mat'],'file')
    warning off
    load([folder 'Analysis/GeneralInfo.mat'],'angOffset');
    warning on
end
if ~exist('angOffset','var')
    warning('Angle not calibrated!')
    angOffset = 0;
end

t = t*10000;
dx = whl(:,1)-whl(:,3);
dy = whl(:,2)-whl(:,4);

ep = ~(isnan(dx) | isnan(dy));

tg = t;
tg(~ep) = -1;
tg(end) = -1;
GoodRanges = thresholdIntervals(tsd(t,tg),0,'Direction','Above');

ang = atan2(dy(ep),dx(ep))-angOffset;
ang = mod(ang,2*pi);
ang = tsd(t(ep),ang);


