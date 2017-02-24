function cleanx = CleanData(x, rnge, choice, epsilon);

%CLEANDATA
% Usage:
%       cleanx = CleanData(x, rnge, choice, epsilon);
%
% This function will take in a data set x, a specific portion of that data,
% rnge, and look for any points below (choice = 'below') or above (
% (choice = 'above') the threshold epsilon.  It will then get rid of those
% points and create linearly extrapolated points in their place.  It
% returns the cleaned data set cleanx.
%
% Note: This function is to be used ONLY when the user knows for a fact that
% spurious points exist in the data set.  It is not very automatic, but
% assures that only KNOWN false points will be removed from the set. In 
% order to use this, one must look at a plot of the data to determine the 
% portion to look at and the threshold to set.

% Written by Dan Valente
% 10 November 2006

cleanx = x;
dirtyx = x(rnge);

if (strcmp(choice, 'below'))
    goodpoints = find(dirtyx > epsilon);
elseif (strcmp(choice, 'above'))
    goodpoints = find(dirtyx < epsilon);
end

goodpoints = goodpoints + rnge(1) - 1;
for i = 1:(length(goodpoints)-1)
    if ((goodpoints(i+1)-goodpoints(i)) ~= 1)
        xf = x(goodpoints(i+1));
        xs = x(goodpoints(i));
        m = (xf-xs)/(goodpoints(i+1)-goodpoints(i));
        count = goodpoints(i);
        while count < goodpoints(i+1)
            cleanx(count) = m*(count-goodpoints(i))+xs; 
            count = count+1;
        end
    end
end

return;