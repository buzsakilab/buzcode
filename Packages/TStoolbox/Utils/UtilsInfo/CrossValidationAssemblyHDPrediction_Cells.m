function L = CrossValidationPeerHDCheck(Starget,Qassembly,Qtarget,ang,ep)

nbEp = 10;

fprintf('Launching Cross-validated Peer info\n')

epV = regIntervals(ep,nbEp);
dQ = Data(Qtarget);
nbC = size(dQ,2);
L = NaN(nbEp,nbC);

% weights = zeros(nbC,nbEp,nbC-1);

angq = Restrict(ang,Qtarget);

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
            [h,B] = HeadDirectionField(Starget{c},ang,trainingEp);
            warning off
                w = ComputePeerHDPrediction(q,Qassembly,ang,h,B,trainingEp);
%                 weights(c,ii,:) = w;
            warning on

            Qtest = Restrict(Qassembly,testEp);
            qt = Data(Restrict(q,testEp));
            rg = Range(Qtest);
            dt = median(diff(rg))/10000;
            pop = Data(Qtest); %in Hz
            r = rate(Starget{c},trainingEp);

            angq = Restrict(ang,Qtest);
            dAng= Data(angq);
            xx = floor((length(B)-1)*contrast(dAng))+1;
            h = h(:);
            f = h(xx);
            
            Lf = SpkTrainLogLikelihood(qt,dt*f.*modifiedExp(pop*w));
            L0 = SpkTrainLogLikelihood(qt,dt*repmat(r,[size(pop,1) 1]));
            L(ii,c) = Lf-L0;

                if 0
                    disp(L(ii,c))
                figure(1),clf
                plot(rg,f)
                hold on
                plot(rg,modifiedExp(pop*w),'r')
                plot(rg,f.*modifiedExp(pop*w),'g')

                keyboard
                end

        end
    end
end

fprintf('done\n')

L = nansum(L);
L = L/tot_length(ep,'s');
L = L(:);
