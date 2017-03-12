function [dat] = bz_scrubTracking(dat)
% USAGE  
%  [dat] = bz_scrubTracking(dat)
% 
% INPUTS
%   dat - struct created by calling dat = importdata(file.csv) on a csv
%         file exported from a .tak file
%
% OUTPUTS
%   dat - the same struct 1) outliers removed
%                         2) median filtering
%                         3) interpolation over missing frames
%
%
% 
% this functions takes the unedited *.csv exported from Motive, smoothes and
% resorts columns to match the position tracking formatting standards
% 
% David Tingley, 2017
% 
%TODO
% check that units are in meters...


% find extreme outliers outside of tracking volume
dat.data((abs(dat.data(:,7))>mean(dat.data(:,7))+std(dat.data(:,7))*5),7)=nan;
dat.data((abs(dat.data(:,8))>mean(dat.data(:,8))+std(dat.data(:,8))*5),8)=nan;
dat.data((abs(dat.data(:,9))>mean(dat.data(:,9))+std(dat.data(:,9))*5),9)=nan;

% find poorly tracked frames using the error per marker column
dat.data(find(dat.data(:,10)>mean(dat.data(:,10))+std(dat.data(:,10))*10),:)=nan;   

% find dropped frames or large translational shifts and call them nan's here
f = find(diff(dat.data(:,8))==0 | abs(diff(dat.data(:,8)))> mean(std(abs(diff(dat.data(:,8)))))+ std(abs(diff(dat.data(:,8)))) * 10);
dat.data(f+1,8)=nan;
f = find(diff(dat.data(:,9))==0 | abs(diff(dat.data(:,9)))> mean(std(abs(diff(dat.data(:,9)))))+std(abs(diff(dat.data(:,9)))) * 10);
dat.data(f+1,9)=nan;
f = find(diff(dat.data(:,7))==0 | abs(diff(dat.data(:,7)))> mean(std(abs(diff(dat.data(:,7)))))+std(abs(diff(dat.data(:,7)))) * 10);
dat.data(f+1,7)=nan;

% median filtering
dat.data(:,7) = medfilt1(dat.data(:,7),12);
dat.data(:,8) = medfilt1(dat.data(:,8),12);
dat.data(:,9) = medfilt1(dat.data(:,9),12);

% pattern matching interpolation
dat.data(:,7) = fillmissing(dat.data(:,7),'pchip');
dat.data(:,8) = fillmissing(dat.data(:,8),'pchip');
dat.data(:,9) = fillmissing(dat.data(:,9),'pchip');

%  now change column order to new standards
% columns = [7 9 8 3 5 4 6 10 1 2]; % optitrack is Y-oriented instead of Z
% dat.data = dat.data(:,columns);

end