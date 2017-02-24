function [pc,rpc,ev,rev] = wavePCA(cv)
%
% [pc,rpc,ev,rev] = wavePCA(cv)
% 
% Principal Component Analysis of standardized waveforms from a
% given (unstandardized) waveform covariance matrix cv(nSamp,nSamp).
% If input is a cell array of covariance matrices the outputs are
% corresponding cell arrays. 
%
%
% INPUT: 
%      cv ... nSamp x nSamp wavefrom covariance matrix (unnormalized)
%
% OUTPUT:
%      pc ... column oriented principal components (Eigenvectors)
%      rpc ... column oriented Eigenvectors weighted with their relative amplitudes
%      ev ... eigenvalues of SVD (= std deviation of data projected onto pc)
%      rev ... relative eigenvalues so that their sum = 1
%
% PL 1999%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

if isnumeric(cv)
   sd = sqrt(diag(cv));    
   cvn = cv./(sd*sd');    % standardized convariance matrix: diag=(1,1,1,...,1)
   [u,ev,pc] = svd(cvn);
   ev  = diag(ev);
   rev = ev/sum(ev);
   rpc = pc.*rev(:,ones(length(rev),1))';
elseif iscell(cv)
   nCells=length(cv);
   for iC = 1:nCells
      sd = sqrt(diag(cv{iC}));    
      cvn = cv{iC}./(sd*sd');    % standardized convariance matrix: diag=(1,1,1,...,1)
      [u,ev{iC},pc{iC}] = svd(cvn);
      ev{iC}  = diag(ev{iC});
      rev{iC} = ev{iC}/sum(ev{iC});
      rpc{iC} = pc{iC}.*rev{iC}(:,ones(length(rev{iC}),1))';
   end%for
else
   error(' Input must be a covariance matrix or a cell array of covariance matrices');
end%if

   