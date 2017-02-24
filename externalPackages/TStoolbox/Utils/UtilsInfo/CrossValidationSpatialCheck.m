function L = CrossValidationSpatialCheck(S,lambda,x1,x2,X,Y,ep)

nbEp = 10;

epV = regIntervals(ep,nbEp);
L = zeros(nbEp,1);

for ii=1:nbEp
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=2:nbEp-1
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    [pf,occH] = PlaceFields_NoSmoothing(S,X,Y,trainingEp,x1,x2);
    pfSm = gaussFilt(pf,lambda);
    pfSm(occH==0) = NaN;
    
    
    
    L(ii) = SpkTrainSpatialValuation(pfSm,S,x1,x2,X,Y,testEp,trainingEp);
    
end

L = sum(L)/tot_length(ep,'s');
