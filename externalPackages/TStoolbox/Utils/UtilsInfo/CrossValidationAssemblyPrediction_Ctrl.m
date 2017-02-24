% function [L weights] = CrossValidationPeerCheck(S,Q,ep)
function [L,weights] = CrossValidationAssemblyPrediction_Cells(Starget,Qtarget,Sref,Qref,angBin,ang,ep,lambda,bins)
nbEp = 10;

epV = regIntervals(ep,nbEp);
fprintf('Launching Cross-validated Peer info\n')

nbC = length(Starget);
nbCref = length(Sref);

dQ = Data(Qtarget);

L = NaN(nbEp,nbC);

weights = zeros(nbC,nbEp,size(Data(Qref),2));

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
     
    qTraining = Restrict(Qtarget,trainingEp);
    dAng = Data(Restrict(ang,qTraining));
    
    % Here we assume that bins are one degree...
    xx = floor((length(angBin)-1)*contrast(dAng))+1;
    
    qRef = zeros(length(xx),nbCref);
    for c=1:nbCref
        h = HeadDirectionField(Sref{c},ang,trainingEp);
        h = bins*h(:);
        qRef(:,c) = poissrnd(h(xx),length(xx),1);
    end
    qRef = gaussFilt(qRef,lambda,0);
    qRef = tsd(Range(qTraining),qRef);
    
    for c=1:nbC
        
        q = tsd(Range(Qtarget),dQ(:,c));
        
        if any(dQ(:,c)>0)
            
            warning off
            if ii==1
                w = ComputePeerPrediction(q,qRef,trainingEp);
                weights(c,ii,:) = w;
            else
                w0 = squeeze(mean(weights(c,1:ii-1,:),2));
                w = ComputePeerPrediction(q,qRef,trainingEp,w0);
            end
            warning on
            weights(c,ii,:) = w;
             
            Qtest = Restrict(Qref,testEp);
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

