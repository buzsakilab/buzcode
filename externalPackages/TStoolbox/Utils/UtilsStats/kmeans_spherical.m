function [c, s, quality] = kmeans_spherical(X, K)
%   [c, s] = kmeans_spherical(X, K)
%
%   spherical kmeans, as suggested for example in Dhillon, 2000. 
%   INPUTS: 
%   X:  a n_samples X D matrix where each row is a data point vector.  
%   K:  the number of clusters  
%   OUTPUTS:
%   C: a K x D matrix of the K centroids 
%   S: a 1 x n_samples array of cluster belonging
%   n_iter: the number of K-means iterations performed  
%   the computationally expensive part of the algorithm is implemented in
%   the mex file kmeans_spherical_impl.c
  
% fpbatta 2003
  
  
  
  n_samples = size(X, 1);
  D = size(X, 2);
  xn = sqrt(sum(X.*X, 2));
  
  XN = X ./ repmat(xn, 1, D);
  
  mean_X = mean(XN, 1);
  mean_X = mean_X / sum(mean_X .* mean_X);
  
  
  C_init = repmat(mean_X, K, 1);
  
  eps_init = 0.002;
  
  % generate initial conditions for centroids perturbing the mean vecotr
  % of the data set C_init. the perturbations are normalized so that they
  % have L2 norm  eps_norm
  
  norm_factor = eps_init / sqrt(D/3);
  c1 = norm_factor * (2 * rand(K, D) -1);
  C_init = C_init + c1;

  cn = sqrt(sum(C_init.*C_init, 2));
  C_init = C_init ./ repmat(cn, 1, D);
%  keyboard  

  rp = randperm(n_samples);
  
  np = floor(n_samples/K);
  qv1 = XN(rp,:);
  
  for i = 1:K
    C_init(i,:) = sum(qv1((1+(i-1)*np):(i*np),:),1);
  end
  
  cn = sqrt(sum(C_init.*C_init, 2));
  C_init = C_init ./ repmat(cn, 1, D);
  

  [c, s, quality] = kmeans_refine_impl(XN, K, C_init, 1000);
  
  %alternate code that stops at each iteration, allowing inspection
% $$$   for i = 1:50
% $$$   
% $$$     [c, s, n_iter] = kmeans_refine_impl(XN, K, C_init, 2);
% $$$   
% $$$     q = mean(sum(XN .* c(s,:), 2), 1)
% $$$     Q(i) = q;
% $$$     keyboard
% $$$     
% $$$     C_init = c;
% $$$   end
% $$$   
  
  
  
  
  
  
  
  
  
  