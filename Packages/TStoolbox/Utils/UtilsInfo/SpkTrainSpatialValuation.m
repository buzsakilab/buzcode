function L = SpkTrainSpatialValuation(pf,S,x1,x2,X,Y,testEp,trainingEp)


X = Restrict(X,testEp);
xs = Data(X);
ys = Data(Restrict(Y,testEp));

t = Range(X);
f = zeros(length(t),1);
r = rate(S,trainingEp);
S = Restrict(S,testEp);

xx = floor((length(x1)-1)*contrast(xs))+1;
yy = floor((length(x2)-1)*contrast(ys))+1;

ix = length(x2)*(xx-1)+yy;
pfl = pf(:);

f = pfl(ix);
f(isnan(f))=0;
%f = tsd(t,f);

L = SpkTrainValuation(S,f,t,r);

if 0
    figure(10),clf
    imagesc(x1,x2,pf),axis xy
    hold on,plot(xs,ys,'w')
    xst = Data(Restrict(X,S));
    yst = Data(Restrict(Y,S));
     hold on,plot(xst,yst,'r*')
     
     figure(11),clf
    imagesc(x1,x2,pf2),axis xy
    hold on,plot(Data(Restrict(X,trainingEp)),Data(Restrict(Y,trainingEp)),'w')
    xst = Data(Restrict(X,Restrict(S,trainingEp)));
    yst = Data(Restrict(Y,Restrict(S,trainingEp)));
     hold on,plot(xst,yst,'r*')
     
     
end






