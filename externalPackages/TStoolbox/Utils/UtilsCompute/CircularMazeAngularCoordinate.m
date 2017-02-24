function phi = CircularMazeAngularCoordinate(X, Y, center, zero_ang)
% phi = CircularMazeAngularCoordinate(X, Y)
%
% gives an angular coordinate for the Warp circular maze used for rats
% 7027, 6985
% INPUTS
% X, Y: tsds containing rat's cartesian coords.
% OUTPUT
% phi: tsd containing the angular coord


if nargin < 3  
  %center = [400 239]; % determined for (7027) 2001-7-22 dataset (best tracking
                    % data so far)
		    center = [325 288] % for 6985 -- 2001-5-1		    



		    %zero_ang = +1.6; %this is the angular coordinate of the point in between
                 %the two reward location. It will be taken as zero in
                 %the output tsd
		 zero_ang = 6.13; % for 7166 trinagle maze 2002-4-22
end

t = Range(X, 'ts');
pos = ([ Data(X), Data(Y)] - repmat(center, length(Data(X)), 1));


p = atan2(pos(:,1), pos(:,2));
%p = p+pi;
tp = find(p < zero_ang);

p(tp) = p(tp) + (2 * pi);

p = p - zero_ang;




phi = tsd(t, p);