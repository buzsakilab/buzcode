% function [L weights] = CrossValidationPeerCheck(S,Q,ep)
function [L,weights] = CrossValidationPeerCheck_Cells(S,Q,ep)
nbEp = 10;

epV = regIntervals(ep,nbEp);
fprintf('Launching Cross-validated Peer info\n')

dQ = Data(Q);
nbC = size(dQ,2);
L = NaN(nbEp,nbC);

weights = zeros(nbC,nbEp,nbC-1);
% compTime = zeros(nbC,nbEp);

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
        popIx = ~ismember((1:nbC),c);
        Qs = tsd(Range(Q),dQ(:,popIx));
        q = tsd(Range(Q),dQ(:,c));
        
        if any(dQ(:,c)>0)
            
            warning off
            if ii==1
                w = ComputePeerPrediction(q,Qs,trainingEp);
                weights(c,ii,:) = w;
            else
                w0 = squeeze(mean(weights(c,1:ii-1,:),2));
                w = ComputePeerPrediction(q,Qs,trainingEp,w0);
            end
            warning on
            weights(c,ii,:) = w;
%             compTime(c,ii) = t;
            
            Qtest = Restrict(Qs,testEp);
            qt = Data(Restrict(q,testEp));
            rg = Range(Qtest);
            dt = median(diff(rg))/10000;
            pop = Data(Qtest);
            r = rate(S{c},trainingEp);
            
            Lf = SpkTrainLogLikelihood(qt,dt*modifiedExp(pop*w));
            L0 = SpkTrainLogLikelihood(qt,dt*repmat(r,[size(pop,1) 1]));
            L(ii,c) = Lf-L0;
            
        end
        
    end

end

fprintf('done\n')

L = nansum(L);
L = L/tot_length(ep,'s');
L = L(:);

weights = squeeze(mean(weights,2));

