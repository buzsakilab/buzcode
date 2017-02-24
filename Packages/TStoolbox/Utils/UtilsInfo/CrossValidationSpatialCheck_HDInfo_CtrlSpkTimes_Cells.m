function L = CrossValidationSpatialCheck_HDInfo_CtrlSpkTimes_Cells(S,lambda,X,Y,ep,varargin)

nbEp = 10;
dispBar = 1;
binSize = 0.5; %in cms
nbRand = 10;

exDisp = 0;
if ~isempty(varargin)
    exDisp = 1;
end

epV = regIntervals(ep,nbEp);
L = zeros(nbEp,length(S),nbRand);
      
S = Restrict(S,ep);
X = Restrict(X,ep);
Y = Restrict(Y,ep);

fprintf('Launching Cross-validated spatial info for jittered spike times\n')
for ii=1:1%nbEp
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=2:nbEp-1
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    fprintf('..')
     
    for c=1:length(S)
    %fprintf('\b%d',c)
   
        l = zeros(nbRand,1);
        for rd=1:nbRand
            Sct = [];
            st = Start(trainingEp);
            for e=1:length(Start(trainingEp))
                isi = diff(Range(Restrict(S{c},subset(trainingEp,e))));
                randp = randperm(length(isi));
                s = st(e) + cumsum(isi(randp));
                if isempty(Sct)
                    Sct = s;
                else
                    Sct = [Sct;s];
                end
            end
            Sct = ts(sort(Sct));

            if c==1 & rd == 1
                [pf,occH,bx,by] = PlaceFields_NoSmoothing(Sct,X,Y,trainingEp,binSize);
            else
                pf = PlaceFields_NoSmoothing(Sct,X,Y,trainingEp,bx,by,occH);
            end

            pfSm = gaussFilt(pf/nbRand,lambda/binSize);
            L(ii,c,rd) = SpkTrainSpatialValuation(pfSm,S{c},bx,by,X,Y,testEp,trainingEp);
        end
        
        if exDisp
            pf = PlaceFields_NoSmoothing(S{c},X,Y,trainingEp,bx,by,occH);
            x = Data(Restrict(X,testEp));
            y = Data(Restrict(Y,testEp));
            xs = Data(Restrict(X,Restrict(S{c},testEp)));
            ys = Data(Restrict(Y,Restrict(S{c},testEp)));
            
            figure(1),clf
            subplot(1,3,1)
                imagesc(bx,by,gaussFilt(pf,5/0.5));
                hold on
                plot(x,y,'w')
                plot(xs,ys,'ks','MarkerFaceColor','k','MarkerSize',4)
            subplot(1,3,2)
                imagesc(bx,by,pfSm);
                hold on
                plot(x,y,'w')
                plot(xs,ys,'ks','MarkerFaceColor','k','MarkerSize',4)
            subplot(1,3,3)
                polar(B,angHisto{c})
            keyboard
            
            if 0
                epswrite([figDir 'SpatialHDInfo_Example_c49'],'size','screen')
                
            end
            
        end
        
        
        
    end
    
end
fprintf('done\n')
% close(h)

if length(size(L))>2
    L = squeeze(sum(L));
    L = nanmean(L,2);
else
    L = sum(L);
    L = L(:);
end
L = L/tot_length(ep,'s');
% L = L(:);