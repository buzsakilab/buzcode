function [ripTsd ripEp] = LoadRipples(fbasename);

% USAGE
%     [ripTsd ripEp] = LoadRipples(fbasename);
% wrapper to load ripples from fbasename.sts.RIPPLES file
%     
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