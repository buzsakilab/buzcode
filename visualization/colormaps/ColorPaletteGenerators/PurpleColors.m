function colors = PurpleColors(numcolors)
% Generates RGB Triplets based on red-blue continuum:
% 
% INPUT
% -numcolors: number of colors you want out
% OUTPUT
% - colors a 3 column matrix signifying Red Green and Blue values with each
% column, number of rows = numcolors
% 
% Brendon Watson 2015

r = linspace(0,1,numcolors)';
g = zeros(numcolors,1);
b = linspace(1,0,numcolors)';

colors = [r g b];
