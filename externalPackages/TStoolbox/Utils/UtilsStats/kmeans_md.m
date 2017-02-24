

%KMEANS_MD K-means clustering
%
%       [c s] = KMEANS_MD(x,k)
%       [c s] = KMEANS_MD(x,k, x0)
%       [c s] = KMEANS_MD(x,k, x0)
%
%       K-means clustering for multi dimensional data, x (each row is a
%       sample) .  k is the number of
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
function [c,s] = kmeans_md(x, K, z)
        deb = 1;

	[n_samples, D] = size(x)
	
        if nargin == 2,
                %
                % if no initial clusters are given spread the
                % centers evenly over the range min to max.
                %
		z = zeros(K,D);
		for ii = 1:D
                z(:,ii) = rand(K,1)*(max(x(:,ii))-min(x(:,ii))) + ...
		    min(x(:,ii));
		end
        end
        if size(z, 1) ~= K,
                error('initial cluster length should be k')
        end

        %
        % step 1
        %
        zp = z;
        s = zeros(n_samples, 1);
        n = n_samples;

        iterating = 1;
        k = 1;
        iter = 0;
        while iterating,
                iter = iter + 1;

                %
                % step 2
                %
                for l=1:n,      % for each point l=1..n
		  for ii = 1:K
		    dist(ii) = norm(x(l,:)-zp(ii,:), 2);
		  end
		  
                        [y,ind] = min(dist);
                        s(l) = ind;     % assign index of closest set
                end

                %
                % step 3
                %
                for j=1:K
		  if any(s==j)
                        zp(j,:) = mean( x(find(s==j), :), 2 );
		  end
		  
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
                        plot(z(:,1), z(:,2));
                        pause(.1);
                end
        end
        if deb>0,
                disp('iterations ');
                disp(iter);
        end
        c = z;

