function [ripTsd ripEp] = LoadRipples(fbasename);

% Load ripple times from a .sts.RIPPLES file
%
% USAGE:
%     [ripTsd ripEp] = LoadRipples(fbasename);
%
% INPUT:
%    fbasename: session file basename
%
% OUTPUT:
%    ripTsd: a tsd of ripple times
%    ripEp: a intervalSet of ripple events (beginnings and ends)

% Adrien Peyrache 2012

fname = [fbasename '.sts.RIPPLES'];
if exist(fname,'file')
    rip = load([fbasename '.sts.RIPPLES']);
    ripTsd = tsd(rip(:,2)*10000,rip(:,4:5));
    ripEp = intervalSet(rip(:,1)*10000,rip(:,3)*10000);
else
    ripTsd = tsd([],[]);
    ripEp = intervalSet([],[]);
end
