function L = CrossValidationSpatialCheck_Cells(S,ang,ep)

nbEp = 10;
% dispBar = 1;

epV = regIntervals(ep,nbEp);
L = NaN(nbEp,length(S));
fprintf('Launching Cross-validated HD info\n')

for ii=1:nbEp
    fprintf('.')

    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=2:nbEp-1
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    for c=1:length(S)
        s = Restrict(S{c},trainingEp);
        if length(s)
            [h,B] = HeadDirectionField(s,ang,trainingEp);
            L(ii,c) = SpkTrainHDValuation(h,S{c},B,ang,testEp,trainingEp);
        end
    end
    
end
% close(wb)
fprintf('done\n')

L = nansum(L);
L = L/tot_length(ep,'s');
L = L(:);