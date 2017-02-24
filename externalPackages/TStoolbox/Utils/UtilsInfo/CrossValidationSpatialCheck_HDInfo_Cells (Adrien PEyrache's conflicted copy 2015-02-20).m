function L = CrossValidationSpatialCheck_HDInfo_Cells(S,lambda,X,Y,ep,ang,angHisto,B,varargin)

nbEp = 10;
dispBar = 1;
binSize = 0.5; %in cms
nbRand = 10;

exDisp = 0;
if ~isempty(varargin)
    exDisp = 1;
end

epV = regIntervals(ep,nbEp);
L = NaN(nbEp,length(S),nbRand);
      
aep = Restrict(ang,ep);
da = Data(aep);
da = abs(repmat(da,[1 length(B)-1])-repmat(B(1:end-1)',[size(da,1) 1]));
da(da>pi) = 2*pi-da(da>pi);
[dummy,minIx] = min(da');
minIx = tsd(Range(aep),minIx(:));

st = Start(ep);
en = End(ep);
time = [st(1):10:en(end)]'; %ms resolution
time = Range(Restrict(ts(time),ep));
minIx = Data(Restrict(minIx,time));
fprintf('Launching Cross-validated spatial info w.r.t. HD signal\n')
for ii=1:nbEp
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=2:nbEp-1
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    fprintf('.')
     
    for c=1:length(S)
        for rd=1:nbRand
           h = angHisto(:,c);
            Sct = SpkTrainGenerator(h(minIx),0.001);
            Sct = ts(time(Sct));

            if c==1 & rd==1
                [pf,occH,bx,by] = PlaceFields_NoSmoothing(Sct,X,Y,trainingEp,binSize);
            else
                pf = PlaceFields_NoSmoothing(Sct,X,Y,trainingEp,bx,by,occH);
            end
            pfSm = gaussFilt(pf,lambda/binSize);
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

L = nansum(L);
if length(size(L))>2
    L = squeeze(L);
    L = nanmean(L,2);
else
    L = L(:);
end
L = L/tot_length(ep,'s');
