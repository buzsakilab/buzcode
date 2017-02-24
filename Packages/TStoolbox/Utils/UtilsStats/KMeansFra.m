

%KMEANS K-means clustering
%
%       [c s] = KMEANS(x,k)
%       [c s] = KMEANS(x,k, x0)
%       [c s] = KMEANS(x,k, x0)
%
%       K-means clustering for one dimensional data, x.  k is the number of
%       clusters, and x0 if given is the inital centroid for the clusters.
%
%       On return c[i] is the centroid of the i'th cluster.  s is a vector
%       of length equal to x, whose value indicates which cluster the
%       corresponding element of x belongs to.
%
% REF: Tou and Gonzalez, Pattern Recognition Principles, pp 94
%
%
%       Copyright (c) Peter Corke, 1999  Machine Vision Toolbox for Matlab

%               pic 7/92
%
function [c,s] = kmeans(x, K, z)
        deb = 0;

        if nargin == 2,
                %
                % if no initial clusters are given spread the
                % centers evenly over the range min to max.
                %
                z = [0:K-1]/(K-1)*(max(x)-min(x)) + min(x);
        end
        if length(z) ~= K,
                error('initial cluster length should be k')
        end

        %
        % step 1
        %
        zp = z;
        s = zeros(size(x));
        n = length(x);

        iterating = 1;
        k = 1;
        iter = 0;
        while iterating,
                iter = iter + 1;

                %
                % step 2
                %
                for l=1:n,      % for each point l=1..n
                        [y,ind] = min(abs(x(l) - z));
                        s(l) = ind;     % assign index of closest set
                end

                %
                % step 3
                %
                for j=1:K
                        zp(j) = mean( x(s==j) );
                end

                %
                % step 4
                %
                nm = norm(z - zp);
                if deb>0,
                        nm
                end
                if nm == 0,
                        iterating = 0;
                end
                z = zp;
                if deb>0,
                        plot(z);
                        pause(.1);
                end
        end
        if deb>0,
                disp('iterations ');
                disp(iter);
        end
        c = z;

