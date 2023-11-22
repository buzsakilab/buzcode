function colors = GrayColors(numcolors)
% Generates RGB Triplets based on red-blue continuum:
% 
% INPUT
% -numcolors: number of colors you want out
% OUTPUT
% - colors a 3 column matrix signifying Red Green and Blue values with each
% column, number of rows = numcolors
% 
% Brendon Watson 2015

r = linspace(.7,0,numcolors)';
g = linspace(.7,0,numcolors)';
b = linspace(.7,0,numcolors)';

colors = [r g b];
