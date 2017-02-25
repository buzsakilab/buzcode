function orientation  = FlyOrient(subset_frame, threshold)

%FLYORIENT
%   Usage:
%       orientation = FlyOrient(subset_frame, threshold)
%
% This function takes in a subset frame around the fly (calculated by
% FindFly) and discards the 3D data by placing points where the pixel intensity
% is larger than the user chosen threshold.  Then, Principal Components
% Analysis (by way fot the pca1 function) is performed on the resulting 
% scatter plot to find the direction of maximum variance --- this direction
% is taken to be the fly's (ambiguous) orientation.  
% orientation is a vector consisting of two angles (complements) that comprise 
% the body axis. The first element is an angle in the upper half plane; the
% second element is an angle in the lower half plane.

% Written by Dan Valente
% 11 October 2006

%Normalize frame data by pixel of maximum intensity
subset_frame = subset_frame/max(max(subset_frame));

% Put dots where fly is and do PCA on reduced data set
[rows, cols] = find(subset_frame >= threshold);
rows = length(subset_frame(:,1))-rows+1;
x = [cols';rows'];
[xnew, PC, V, data] = pca1(x);

% Find orientation vectors (two, mirrored across diagonal), and group into
% upper half and lower half planes.
a1 = PC(1,1);
b1 = PC(2,1);
a2 = -PC(1,1);
b2 = -PC(2,1);
if (b1 >= 0 );
    orientUHP = atan2(b1,a1);
    orientLHP = atan2(b2,a2);
elseif (b2 >=0);
    orientUHP = atan2(b2,a2);
    orientLHP = atan2(b1,a1);
else
end

% The vector we will return
orientation = [orientUHP orientLHP];

return;
