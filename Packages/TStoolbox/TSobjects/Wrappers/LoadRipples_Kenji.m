function [ripTsd ripEp] = LoadRipples_Kenji(fbasename,varargin)

% Load ripple times from Kenji dataset
%
% USAGE:
%     [ripTsd ripEp] = LoadRipples_Kenji(fbasename);
%
% INPUT:
%    fbasename: session file basename
%
% OUTPUT:
%    ripTsd: a tsd of ripple times
%    ripEp: a intervalSet of ripple events (beginnings and ends)

% Adrien Peyrache 2013

%parameters
ripFolder = '/mnt/DataAdrien12/KenjiData/Ripple/';
timeFactor = 10000/1250;


load([ripFolder fbasename '.Ripple.mat']);

if isempty(varargin)
    epoch = 'SWS';
else
    epoch = varargin{1};
end

epoch = ['ind_' epoch];
epoch = find(ismember(Ripple_Stat.content,epoch));
epoch = find(Ripple_Stat.Time{3}(:,epoch));

%take negative_peak_time_near_power_peak, 6th index
rip = Ripple_Stat.Time{3}(epoch,6);
%take max sd, 8th column
sd = Ripple_Stat.Time{3}(epoch,8);

ripTsd = tsd(rip*timeFactor,sd);
ripEp = intervalSet(Ripple_Stat.Time{3}(epoch,1)*timeFactor,Ripple_Stat.Time{3}(epoch,2)*timeFactor);
