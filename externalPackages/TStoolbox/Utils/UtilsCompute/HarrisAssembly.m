function w = HarrisAssembly(Q,cellIx,sigma,dt)

N = size(Q,2);
peerIx = [1:N];
peerIx(ismember(peerIx,cellIx))=[];

t = [-200:200]';
gw = 1/sqrt(2*pi*sigma^2)*exp(-t.^2/(2*sigma^2));

Qs = convn(Q,gw,'same');
sa = Qs(:,peerIx);
n = Q(:,cellIx);

l = @(x) sum(-HarrisExp(sa*x)*dt+n.*log(HarrisExp(sa*x)*dt))-sum(x.^2)/4;

%  w0 = zeros(length(peerIx),1);
w0 = randn(length(peerIx),1);
%  keyboard
w = fminunc(l,w0);