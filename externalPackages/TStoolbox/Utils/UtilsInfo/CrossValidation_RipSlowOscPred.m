function L = CrossValidationHDCheck(rip,dwEp,allEp)

nbEp = 10;

epV = regIntervals(allEp,nbEp);
L = zeros(nbEp,1);
Q = MakeQfromS(tsdArray({rip}),50);


for ii=1:nbEp
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=1:nbEp-2
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    
    upSt = allEp - dwEp;
    
    ripRdw = rate(rip,intersect(dwEp,trainingEp));
    ripRup = rate(rip,intersect(upSt,trainingEp));
    
    q = Restrict(Q,testEp);
    t = Range(Q);
    r = zeros(length(t),1);
    [~,ix] = Restrict(ts(t),intersect(dwEp,testEp));
    r(ix) = ripRdw;
    [~,ix] = Restrict(ts(t),intersect(upSt,testEp));
    r(ix) = ripRup;
            
    L(ii) = SpkTrainValuation(Restrict(rip,testEp),r,t,rate(rip,trainingEp));
    
end

L = sum(L)/tot_length(allEp,'s');