function [L,weights] = CrossValidationPeerHDCheck_Cells(S,Q,ang,ep)

nbEp = 10;

epV = regIntervals(ep,nbEp);
L = NaN(nbEp,length(S));
fprintf('Launching Cross-validated Peer info\n')

dQ = Data(Q);
nbC = length(S);
weights = zeros(nbC,nbEp,nbC-1);

angq = Restrict(ang,Q);

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
        s = S{c};
        popIx = ~ismember((1:nbC),c);
        Qs = tsd(Range(Q),dQ(:,popIx));
        q = tsd(Range(Q),dQ(:,c));
        
        if length(s)
            [h,B] = HeadDirectionField(s,ang,trainingEp);
            warning off
            if ii==1
                w = ComputePeerHDPrediction(q,Qs,ang,h,B,trainingEp);
            else
                w0 = squeeze(mean(weights(c,1:ii-1,:),2));    
                w = ComputePeerHDPrediction(q,Qs,ang,h,B,trainingEp,w0);
            end
            warning on
            weights(c,ii,:) = w;
            
            Qtest = Restrict(Qs,testEp);
            qt = Data(Restrict(q,testEp));
            rg = Range(Qtest);
            dt = median(diff(rg))/10000;
            pop = Data(Qtest); %in Hz
            r = rate(s,trainingEp);

            angq = Restrict(ang,Qtest);
            dAng= Data(angq);
            xx = floor((length(B)-1)*contrast(dAng))+1;
            h = h(:);
            f = h(xx);
            
            Lf = SpkTrainLogLikelihood(qt,dt*f.*modifiedExp(pop*w));
            %Lff = SpkTrainLogLikelihood(qt,dt*f);
            L0 = SpkTrainLogLikelihood(qt,dt*repmat(r,[size(pop,1) 1]));
            
            L(ii,c) = Lf-L0;

                if 0
                    disp(L(ii,c))
                figure(2),clf
                plot(rg,f)
                hold on
                plot(rg,modifiedExp(pop*w),'r')
                plot(rg,f.*modifiedExp(pop*w),'g')
                plot(rg,qt/dt,'c')
                keyboard
                end

        end
    end
end


fprintf('done\n')

L = nansum(L);
L = L/tot_length(ep,'s');
L = L(:);
weights = squeeze(mean(weights,2));
