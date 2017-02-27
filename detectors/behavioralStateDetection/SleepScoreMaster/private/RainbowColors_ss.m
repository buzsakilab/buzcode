function colors = RainbowColors_ss(numcolors)
% Generates RGB Triplets based on the idealized colors below:
% 0. red    [1 0 0]
% 1. orange [1 .5 0]
% 2. yellow [1 1 0]
% 2.5. green  [0 1 0]
% 3. cyan   [0 1 1]
% 4. blue   [0 0 1]
% 5. purple [.5 0 .5]
% 
% INPUT
% -numcolors: number of colors you want out
% OUTPUT
% - colors a 3 column matrix signifying Red Green and Blue values with each
% column, number of rows = numcolors
% 
% Brendon Watson 2015

cidxs = linspace(0,1,numcolors);
r = redfunc(cidxs);
g = greenfunc(cidxs);
b = bluefunc(cidxs);

colors = [r g b];

1;

function y = redfunc(x)
for a = 1:length(x)
    if x(a)<2/5
        y(a) = 1;
    elseif x(a)>=2/5 && x(a)<2.5/5
        y(a) = 1-5*(x(a)-2/5);
    elseif x(a)>=2.5/5 && x(a)<4/5
        y(a) = 0;
    elseif x(a)>=4/5
        y(a) = 0+.5*5*(x(a)-4/5);
    end
end
y = y';

function y = greenfunc(x)
for a = 1:length(x)
    if x(a)<2/5
        y(a) = .5*5*(x(a));
    elseif x(a)>=2/5 && x(a)<3/5
        y(a) = 1;
    elseif x(a)>=3/5 && x(a)<4/5
        y(a) = 1-5*(x(a)-3/5);
    elseif x(a)>=4/5
        y(a) = 0;
    end
end
y = y';

function y = bluefunc(x)
for a = 1:length(x)
    if x(a)<3/5
        y(a) = 0;
    elseif x(a)>=2.5/5 && x(a)<3/5
        y(a) = 5*(x(a)-2.5/5);
    elseif x(a)>=3/5 && x(a)<4/5
        y(a) = 1;
    elseif x(a)>=4/5
        y(a) = 1-.5*5*(x(a)-4/5);
    end
end
y = y';
