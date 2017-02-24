function [x, y, bodyline, sqr] = FindFly(chunk, sqrsize)

%FINDFLY
% 
% Usage:
%   [x, y, bodyline, sqr] = FindFly(chunk, sqrsize)
%
% This function takes in a single image matrix (chunk) and finds the fly. 
% First, it finds the brightest pixel in the image, then grabs a square 
% of size 2*sqrsize pixels around the fly.  The center of mass of this 
% square is calculated (mass = pixel level) and that is returned as the 
% fly position. The variable 'sqr' can be returned, which contains the 
% limits defining a bounding box around the fly.
%
% This function also returns the body axis orientation vector. The
% body axis is calcuated using FlyOrient.  Fly Orient returns the bodyline
% vector, which is a 2-element vector containing two angles, one in the upper
% half plane, and one in the lower half plane. These angles (which are 
% complementary) define the body axis. The noise threshold required for
% FlyOrient is currently fixed at 0.2.  If you wish to change this
% threshold, you must go into the FindFly.m file and hard-code a
% different value.


% Written by Dan Valente
% 28 September 2006

noise_thresh = 0.2;    % for the FlyOrient function

brightest_pixel_level = max(max(chunk));
[brightpix_row brightpix_col] = find(chunk >= (brightest_pixel_level));
     
%just in case more spots have the same brightness, only take one...
width = length(chunk(1,:));
height = length(chunk(:,1));

row_pos = brightpix_row(1);
col_pos = brightpix_col(1);
        
y = row_pos;   
x = col_pos;
 
        
%take subset of pixels around fly.  need to take enough to ensure
%entire fly is captured. 
row_lower_limit = row_pos-sqrsize;
row_upper_limit = row_pos+sqrsize;
if (row_lower_limit <= 0)
    row_lower_limit = 1;
end
if (row_upper_limit >= height )
    row_upper_limit = height;
end
        
col_lower_limit = col_pos-sqrsize;
col_upper_limit = col_pos+sqrsize;
if (col_lower_limit <= 0)
    col_lower_limit = 1;
end
if (col_upper_limit >= width )
    col_upper_limit = width;
end

sqr = [row_lower_limit row_upper_limit col_lower_limit col_upper_limit];

%Now grab center of mass of fly (i.e. CM of subset image pixel intensities)  
temp_mat = chunk(row_lower_limit:row_upper_limit,col_lower_limit:col_upper_limit);

x2 = [1:length(temp_mat(1,:))]';
y2 = [1:length(temp_mat(:,1))]';
total = sum(sum(temp_mat));
x = sum(temp_mat*x2)/total+col_lower_limit-1;
y = sum(temp_mat'*y2)/total+row_lower_limit-1;


% Find Fly orientation from this frame (in radians)
bodyline = FlyOrient(temp_mat, noise_thresh);  

return;