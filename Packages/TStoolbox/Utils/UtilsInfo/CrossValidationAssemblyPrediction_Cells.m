function [L,weights] = CrossValidationAssemblyPrediction_Cells(Starget,Qassembly,Qtarget,ep)

% 
% [L,weights] = CrossValidationAssemblyPrediction_Cells(Starget,Qassembly,Qtarget,ep)
% 
% computs the loglikelihood prediction of the spike train Starget from the
% binned spike train matrix of its peers (Qassembly). The prediction is built
% upon a Generalized Linear Model of the multivariate data given by
% Qassembly. The algorithm searches for a set of weights that maximize the
% prediction of the binned spike train Qtarget (the training set) and then
% computes the likelihood (during the test set). Everything is the same as
% in Harris et al., Nature, 2004.
% 
% INPUTS:
%     Starget: a ts object of  spike times of the target cell
%     Qasssembly: a binned spike train matrix of the predictive cells (could be smoothed)
%     Qtarget: the binned spike train of the target cell, same bins and same smoothing
%     ep:an intervalset objet of the epoch during which the prediction is built
%     
% OUTPUTS:
%     L: information rate (bit/s)
%     weights: the optimal weights

%Adrien Peyrache 2014

% Parameters:
nbEp = 10; %Number of cross-validations

epV = regIntervals(ep,nbEp);
fprintf('Launching Cross-validated Peer info\n')

dQ = Data(Qtarget);
nbC = size(dQ,2);

L = NaN(nbEp,nbC);

weights = zeros(nbC,nbEp,size(Data(Qassembly),2));

for ii=1:nbEp
    fprintf('.')
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=1:nbEp-2
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    for c=1:nbC
        
        q = tsd(Range(Qtarget),dQ(:,c));
        
        if any(dQ(:,c)>0)
            
            warning off
            if ii==1
                w = ComputePeerPrediction(q,Qassembly,trainingEp);
                weights(c,ii,:) = w;
            else
                %Use previous weights as initial conditions to improve speed
                w0 = squeeze(mean(weights(c,1:ii-1,:),2));
                w = ComputePeerPrediction(q,Qassembly,trainingEp,w0);
            end
            warning on
            weights(c,ii,:) = w;
             
            %Now compute prediction
            Qtest = Restrict(Qassembly,testEp);
            qt = full(Data(Restrict(q,testEp)));
            rg = Range(Qtest);
            dt = median(diff(rg))/10000;
            pop = full(Data(Qtest));
            r = rate(Starget{c},trainingEp);
            
            Lf = SpkTrainLogLikelihood(qt,dt*modifiedExp(pop*w));
            L0 = SpkTrainLogLikelihood(qt,dt*repmat(r,[size(pop,1) 1]));
            
%             disp(Lf-L0)
            L(ii,c) = Lf-L0;
            
        end
    end
end

fprintf('done\n')

L = nansum(L);
L = L/tot_length(ep,'s');
L = L(:);

weights = squeeze(mean(weights,2));

