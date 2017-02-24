function [ang,GoodRanges,ep] = HeadDirection_Wrapper(fbasename,varargin)


% loads head-direction from a position file file (ending in .whl)
%
% USAGE
%     [ang,GoodRanges] = HeadDirectionWhl(fbasename,options)
%     
% INPUT:
%     fbasename: session file basename
%     whlStruct (optional): a structure containing position info (see
%                           LoadPosition_wrapper)
%     angOffset (optional): boolean for angular offset correction (see
%                           CalibrateLEDorientation) (default = 1)
% OUTPUT
%     ang: a tsd object of angular values
%     GoodRanges: a intervalSet object where LEDs were successfully detected


% Adrien Peyrache 2011

whl = [];
angCorrection = 1;

if ~isempty(varargin)
    w = varargin{1};
    
    if ~isa(w,'struct')
        if isa(w,'numeric') | isa(w,'logical')
            angCorrection = w ~= 0;
        else
            error('argument must be a struct or numeric/logical')
        end
    else
        whl = w.whl;
        t = w.t;
        GoodRanges = w.GR;
    end
    
end

if isempty(whl)
    [whl,t,GoodRanges] = LoadPosition(fbasename);
end

[folder dumy dumy]  = fileparts(fbasename);
if ~isempty(folder)
    folder =  [folder filesep];
end

angOffset = 0;

if angCorrection
    if exist([folder 'Analysis/GeneralInfo.mat'],'file')
        warning off
        load([folder 'Analysis/GeneralInfo.mat'],'angOffset');
        warning on
    end
end

%t = t*10000;
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


