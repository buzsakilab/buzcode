function L = CrossValidationGLM_Cells(Q,Pred,ep)

nbEp = 10;
dispBar = 0;
binSize = 0.5; %in cms

epV = regIntervals(ep,nbEp);
L = NaN(nbEp,length(S),length(lambda));
fprintf('Launching Cross-validated spatial info\n')

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
            if c==1;
                [pf,occH,bx,by] = PlaceFields_NoSmoothing(s,X,Y,trainingEp,binSize);
            else
                pf = PlaceFields_NoSmoothing(s,X,Y,trainingEp,bx,by,occH);
            end
            
            for ll=1:length(lambda)
                pfSm = gaussFilt(pf,lambda(ll)/binSize);
                try
                L(ii,c,ll) = SpkTrainSpatialValuation(pfSm,S{c},bx,by,X,Y,testEp,trainingEp);
                catch
                    keyboard
                end
                
            end
        end
    end
    
end
fprintf('done\n')
L = nansum(L);
if length(size(L))>2
    L = squeeze(L);
else
    L = L(:);
end
    
L = L/tot_length(ep,'s');