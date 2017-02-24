function L = CrossValidationPeerCheck(S,Q,ep)

nbEp = 10;

epV = regIntervals(ep,nbEp);
L = zeros(nbEp,1);
weights = zeros(size(Data(Q),2),nbEp);

for ii=1:nbEp
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=1:nbEp-2
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    warning off
    w = ComputePeerPrediction(S,Q,trainingEp);
    weights(:,ii) = w;
    warning on
    
    Qtest = Restrict(Q,testEp);
    rg = Range(Qtest);
    pop = Data(Qtest);
    r = rate(S,trainingEp);
    
    L(ii) = SpkTrainValuation(Restrict(S,testEp),modifiedExp(pop*w),rg,r);
    
end

L = sum(L)/tot_length(ep,'s');