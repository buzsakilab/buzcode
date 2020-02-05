function b = getCubicBSplineBasis(x,nknots,isCirc)
% Calculates set of splines with n evenly spaces knots to describe
% fluctuations in x
% Author: Abed Ghanbari

% input
% x = timestamps of neuron normalized by range
% nknots = number of knots required to evens sample time range
% isCirc = wrap time around ?


%output = beta weights for the nknots

% nknots = 3;
% x = linspace(0,1,256);
% isCirc=false;
% x = linspace(0,2*pi,256);
% isCirc=true;

if ~isCirc & nknots>1
    nknots=nknots-1; % for consistency across circ/non-circ
    weights = ones(nknots+1,1);
    knots = linspace(-2/nknots,1+2/nknots,nknots+5); %define know positions
    s = fastBSpline(knots,weights); 
    b = s.getBasis(x);
else
    knots = linspace(-2*2*pi/nknots,2*pi+2*2*pi/nknots,nknots+5);
    weights = ones(nknots+1,1);
    s = fastBSpline(knots,weights);
    b = zeros(length(x),nknots+1);
    for k=-4:4
        xk = x+k*2*pi;
        b = b+s.getBasis(xk);
    end
    b = b(:,1:end-1);
end
% plot(x,b)
% length(knots)-length(weights)-1
% axis tight
% imagesc(b)

if nknots>1
    b = [b(:,1)*0+1 b];
end