function weights = ComputePeerHDPrediction(q,Q,ang,pf,x1,ep,varargin)

% weights = ComputePeerHDPrediction(S,Q,ep)
% S: ts of the cell
% Q: binned spike train matrix (tsd object) of all other cells
% ang: tsd of angles 
% ep: intervalset of interest

q = Data(Restrict(q,ep));
Q = Restrict(Q,ep);
rg = Range(Q);
dt = median(diff(rg))/10000;
pop = Data(Q); %In Hz

ang = Restrict(ang,Q);
dAng= Data(ang);
xx = floor((length(x1)-1)*contrast(dAng))+1;
pf = pf(:);
pf = pf;%/max(pf);
f = pf(xx);

if ~isempty(varargin)
    weights0 = varargin{1};
else
    weights0 = randn(size(pop,2),1);
end

L = @(x)-SpkTrainLogLikelihood(q,dt*f.*modifiedExp(pop*x)) + 0.25*sum(x.^2);
options = optimoptions('fminunc','display','off');

problem.objective = L;
problem.options = options;
problem.solver = 'fminunc';
problem.x0 = weights0;

weights = fminunc(problem);


if 1
 
figure(1),clf
plot(rg,f)
hold on
plot(rg,modifiedExp(pop*weights),'r')
plot(rg,f.*modifiedExp(pop*weights),'g')

end

