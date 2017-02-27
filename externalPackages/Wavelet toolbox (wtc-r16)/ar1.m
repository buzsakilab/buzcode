function [g,a,mu2]=ar1(x)
% AR1 - Allen and Smith AR(1) model estimation.
% Syntax: [g,a,mu2]=ar1(x);
%
% Input:  x - time series (univariate).
%
% Output: g - estimate of the lag-one autocorrelation.
%         a - estimate of the noise variance.[unbiased_var(x)~=a^2/(1-g^2)]
%         mu2 - estimated square on the mean.
%
% AR1 uses the algorithm described by Allen and Smith 1995, except
% that Matlab's 'fzero' is used rather than Newton-Raphson.
%
% Fzero in general can be rather picky - although
% I haven't had any problem with its implementation
% here, I recommend occasionally checking the output
% against the simple estimators in AR1NV.
%
% Alternative AR(1) estimatators: ar1cov, ar1nv, arburg, aryule
%
% Original version written by Eric Breitenberger.
%
% Updated,optimized&stabilized by Aslak Grinsted 2003-2007
%


x=x(:);
N=length(x);
m=mean(x);
x=x-m;

% Lag zero and one covariance estimates: 
c0=x'*x/N;
c1=x(1:N-1)'*x(2:N)/(N-1);

%how the expression below was derived:
%mu2=(-1/N+(2/N^2)*((N-g^N)/(1-g)-g*(1-g^(N-1))/(1-g)));
%solve('c0t=c0+c0t*(mu2)','c0t')
%c0t=(c0*N^2*(-1+g)/(-N^2+g*N^2+N+N*g-2*g));
%solve('g=(c1+c0t*mu2)/(c0+c0t*mu2)','g')
%g=1/2/c0/N^2*(c1*N+c0*N^2+2*c0-2*c1+c1*N^2-c0*N-(4*c0^2+4*c1^2-4*c1^2*N-6*c0^2*N^3+5*c0^2*N^2+c1^2*N^4+c0^2*N^4-8*c0*c1-4*c0^2*N-3*c1^2*N^2+2*c1^2*N^3-2*c0*c1*N^2-2*c0*N^4*c1+8*c0*c1*N+4*c0*N^3*c1)^(1/2));
%other solution:
%g(2)= 1/2/c0/N^2*(c1*N+c0*N^2+2*c0-2*c1+c1*N^2-c0*N+(4*c0^2+4*c1^2-4*c1^2*N-6*c0^2*N^3+5*c0^2*N^2+c1^2*N^4+c0^2*N^4-8*c0*c1-4*c0^2*N-3*c1^2*N^2+2*c1^2*N^3-2*c0*c1*N^2-2*c0*N^4*c1+8*c0*c1*N+4*c0*N^3*c1)^(1/2));

%A.Grinsteds substitutions: 
B=-c1*N-c0*N^2-2*c0+2*c1-c1*N^2+c0*N;
A=c0*N^2;
C=N*(c0+c1*N-c1);
D=B^2-4*A*C;
if D>0
    g=(-B-sqrt(D))/(2*A);
else
    warning('REDNOISE:unboundAr1','Can not place an upperbound on the unbiased AR1.\n\t\t -Series too short or too large trend.');
    g=nan; 
end


if nargout>1
    mu2=-1/N+(2/N^2)*((N-g^N)/(1-g)-g*(1-g^(N-1))/(1-g)); %allen&smith96(footnote4)
    c0t=c0/(1-mu2);
    a=sqrt((1-g^2)*c0t);
end



%--------------------OLD STUFF---------------------------
%     function g=oldar1(c0,c1,N)
%         g0=c1/c0; % Initial estimate for gamma
% 
% 
%         o=optimset('fzero');
%         o.TolX=1e-5;
%         o.Display='off';
% 
% 
%         % Find g by getting zero of 'gammest':
%         g=fzero(@gammest,g0,o); %-updated to accomodate newer versions of matlab
%         if isnan(g) %optimization failed use homemade zero finder... (aslak 2005)
%             ginout=[g0*.999;g0];
%             ginout(:,2)=[gammest(ginout(1,1)); gammest(ginout(2,1));];
%             ii=0;
%             while (ii<10)&(abs(ginout(end,2))>1e-6)
%                 ii=ii+1;
%                 dx=ginout(end,1)-ginout(end-1,1);
%                 dy=ginout(end,2)-ginout(end-1,2);
%                 xguess=ginout(end,1)-ginout(end,2)*dx/dy; %linear extrap
%                 ginout(end+1,:)=[xguess gammest(xguess)];
%             end
%             [m,mi]=nanmin(abs(ginout(:,2)));
%             g=ginout(mi,1);
%             if (m>1e-6)
%                 g=nan;
%             end
%         end
% 
%         function gout=gammest(gin)
%             % GAMMEST - used by AR1 to compute
%             % a function for minimization by fzero.
%             %
%             % Written by Eric Breitenberger.      Version 1/21/96
%             % Please send comments and suggestions to eric@gi.alaska.edu
%             %
% 
%             %gk=1:N-1; gk=gin.^gk; %old slow gk
%             %gk=cumprod(repmat(gin,1,N-1)); %optimization - aslak 2007
%             %mu2=(1/N)+(2/N^2)*sum((N-1:-1:1).*gk);
%             mu2=-1/N+(2/N^2)*((N-gin^N)/(1-gin)-gin*(1-gin^(N-1))/(1-gin));
%             %^^ aslaks07 optimization: explicit sum from allen&smith-jclim96(footnote4)
%             gout=(1-g0)*mu2+g0-gin;
%             if gout>1, gout=nan; end
%         end
%     end