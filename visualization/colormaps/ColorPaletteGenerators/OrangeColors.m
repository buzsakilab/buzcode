function colors = OrangeColors(numcolors)
% Generates RGB Triplets based on red-blue continuum:
% 
% INPUT
% -numcolors: number of colors you want out
% OUTPUT
% - colors a 3 column matrix signifying Red Green and Blue values with each
% column, number of rows = numcolors
% 
% Brendon Watson 2015

r = ones(numcolors,1);
g = linspace(1,0,numcolors)';
b = zeros(numcolors,1);

colors = [r g b];
