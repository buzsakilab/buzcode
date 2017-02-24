function em = errormatrix(fet,clu)

cluIx = unique(clu);
nClu = length(cluIx);
nDim = size(fet,2);

mu = zeros(nClu,nDim);
sigma = zeros(nDim,nDim,nClu);
p = zeros(1,nClu);
for ii=1:nClu
    f = fet(clu==cluIx(ii),:);
    % Mean fet vector
    mu(ii,:) = mean(f);
    % Covariance. We add 1e-5 to the diagonal to ensure that it will be
    % positive definite
    sigma(:,:,ii) = cov(f)+1e-5*eye(nDim,nDim);
    p(ii) = sum(clu==cluIx(ii));
end
obj = gmdistribution(mu,sigma,p);
%keyboard

em = zeros(nClu,nClu);

for ii=1:nClu
    p = posterior(obj,fet(clu==cluIx(ii),:));
    %keyboard
    em(ii,:) = mean(p);
end
