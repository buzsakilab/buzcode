function ResampleBinary(inputName,outputName,nChannels,up,down)

%ResampleBinary - Resample binary data file.
%
% Resample binary data file, e.g. create LFP file from raw data file.
%
%  USAGE
%
%    ResampleBinary(inputName,outputName,nChannels,up,down)
%
%    inputName      binary input file
%    outputName     binary output file
%    nChannels      number of channels in the file
%    up             upsampling integer factor
%    down           downsampling integer factor
%
%  NOTE
%
%    The actual resampling ratio is up/down.
%
%    Here is a list of typical values for Spike2 recording systems:
%
%    FROM           TO       UP     DOWN
%    =====================================
%    20000          1025     1      16
%    19531.25       20000    128    125
%    19531.25       1025     128    125*16
%    19841.29...    20000    126    125
%    19841.29...    1025     126    125*16
%    20284          20000    5000   5071
%    20284          1025     5000   5071*16

% Copyright (C) 2004-2011 by Michaël Zugaro, 2004 by Hajime Hirase
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin ~= 5,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ResampleBinary">ResampleBinary</a>'' for details).');
end

%  if nargin == 4,
%  	mode = up;
%  	if ~isa(mode,'char'),
%  		error(['Incorrect mode (type ''help <a href="matlab:help ResampleBinary">ResampleBinary</a>'' for details).']);
%  	end
%  	switch(lower(mode)),
%  		case '19to20',
%  			up = 128;
%  			down = 125;
%  		case '19tolfp',
%  			up = 128;
%  			down = 125*16;
%  		case '20tolfp',
%  			up = 1;
%  			down = 16;
%  		otherwise,
%  			error(['Incorrect mode ''' mode '''(type ''help <a href="matlab:help ResampleBinary">ResampleBinary</a>'' for details).']);
%  	end
%  end

% Open input file and output file
inputFile = fopen(inputName,'r');
outputFile = fopen(outputName,'w');

%
bufferSize = 2^16  - mod(2^16,down); % 16 or 12?
% Number of overlapping points per channel in the resampled signal
% (chosen so that both resampledOverlap and originalOverlap are integers)
resampledOverlap = 8*up;
% Number of overlapping points per channel in the original signal
originalOverlap = resampledOverlap * down/up;

% Read first buffer
[overlapBuffer,count] = fread(inputFile,[nChannels,originalOverlap],'int16');
overlapBuffer = fliplr(overlapBuffer);
frewind(inputFile);
[dataSegment,count] = fread(inputFile,[nChannels,bufferSize],'int16');
dataSegment2 = [overlapBuffer,dataSegment]';
resampled = resample(dataSegment2,up,down);
count2 = fwrite(outputFile,resampled(resampledOverlap+1:size(resampled,1)-resampledOverlap/2,:)','int16');
overlapBuffer = dataSegment2(size(dataSegment2,1)-(originalOverlap-1):size(dataSegment2,1),:);

% Read subsequent buffers
while ~feof(inputFile),
  [dataSegment,count] = fread(inputFile,[nChannels,bufferSize],'int16');
  dataSegment2 = [overlapBuffer;dataSegment'];
  resampled = resample(dataSegment2,up,down);
  count2 = fwrite(outputFile,resampled((resampledOverlap/2+1):size(resampled,1)-resampledOverlap/2,:)','int16');
  overlapBuffer = dataSegment2(size(dataSegment2,1)-(originalOverlap-1):size(dataSegment2,1),:);
end

% Add the last unprocessed portion
resampled = resample(overlapBuffer,up,down);
count2 = fwrite(outputFile,resampled((resampledOverlap/2+1):end,:)','int16');

fclose(outputFile);

