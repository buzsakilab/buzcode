function varargout = ConfineTspOutput(fbasename,outputfilename,varargin)

% USAGE:
%     ConfineTspOutput(fbasename,varargin)
% 
% Keep only LED points detected within the user-specified polygon (ie home
% cage).
% >Assumes you will either enter inputs to specify shape of polygon, or that
% .mpg file is in same directory as the .tsp
%
% INPUT:
% 'fbasename': the file base names ('fbasename.tsp', etc.)
% 'outputfilename': name of file to which data will be written (so as not
% to overwrite original data)
%  varargin: may contain an n by 2 matrix with x and y coordinates
%  respectively in each column specifying the points of a polygon within 
%  which detected points must lie.  If they are outside this polygon, they
%  will be set to -1 -1.
%   


% Brendon Watson 2013

if strcmp(fbasename(end-12:end),'_original.tsp')
    tspname = fbasename;
    fbasename = fbasename(1:end-13);
elseif strcmp(fbasename(end-3:end),'.tsp')
    tspname = fbasename;
    fbasename = fbasename(1:end-4);
else
    tspname = [fbasename,'.tsp'];
end

if ~isempty(varargin)
    polycoords = varargin{1};
    x = polycoords(:,1);
    y = polycoords(:,2);
else
    mpgname = [fbasename,'.mpg'];
%     obj = VideoReader(mpgname);
%     frame = read(obj,1);
    mpgpath = fullfile(cd,mpgname);%HAVE to use full path with mmread (below)
    [frame,audio]=mmread(mpgpath,1);
    frame = frame(1).frames.cdata;
    h = figure;
    imagesc(frame);
    [x,y] = ginput;
    disp('Hit Return when finished selecting points')
%     polycoords = [x y];
    close(h)
end

%read tsp
tspfile = [tspname];
tspdata = load(tspfile);

%gather x,y data for each timepoint for each color
color1x = tspdata(:,2);
color1y = tspdata(:,3);
color2x = tspdata(:,4);
color2y = tspdata(:,5);
color3x = tspdata(:,6);
color3y = tspdata(:,7);

%determine which points are in the area selected (area of animal)
color1inout = inpolygon(color1x,color1y,x,y);
color2inout = inpolygon(color2x,color2y,x,y);
color3inout = inpolygon(color3x,color3y,x,y);

% if outside this set to -1,-1
tspdata(~color1inout,2)=-1;
tspdata(~color1inout,3)=-1;
tspdata(~color2inout,4)=-1;
tspdata(~color2inout,5)=-1;
tspdata(~color3inout,6)=-1;
tspdata(~color3inout,7)=-1;

figure;
hold on
plot(tspdata(:,6),tspdata(:,7),'r.')
plot(tspdata(:,4),tspdata(:,5),'g.')
plot(tspdata(:,2),tspdata(:,3),'b.')
title(fbasename)

% save to output file
if ~strcmp(outputfilename(end-3:end),'.tsp')
    outputfilename = [outputfileename,'.tsp'];
end

dlmwrite(outputfilename,tspdata,'delimiter',' ','newline','pc','precision',12)

if nargout == 1
    varargout{1} = [x y];
else
    varargout = [];
end