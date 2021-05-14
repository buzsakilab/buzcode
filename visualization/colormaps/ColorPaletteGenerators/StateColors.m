function colors = StateColors
% Generates RGB Triplets based on red-blue continuum:
% 
% INPUT
% -numcolors: number of colors you want out
% OUTPUT
% - colors a 3 column matrix signifying Red Green and Blue values with each
% column, number of rows = numcolors
% 
% Brendon Watson 2015

colors(1,:) = [0 0 0];
colors(2,:) = [0 0 1];
colors(3,:) = [1 0 0];
