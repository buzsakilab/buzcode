function L = CrossValidationHDCheck(S,ang,ep)

nbEp = 10;

epV = regIntervals(ep,nbEp);
L = zeros(nbEp,1);

for ii=1:nbEp
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=1:nbEp-2
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    [h,B] = HeadDirectionField(S,ang,trainingEp);

    L(ii) = SpkTrainHDValuation(h,S,B,ang,testEp,trainingEp);
    
end

L = sum(L)/tot_length(ep,'s');