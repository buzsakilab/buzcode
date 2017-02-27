function [meantheta,R,sigma,confangle,kappa]=anglemean(theta);
% calculates the mean of angles
%
%    [meantheta,anglestrength,sigma,confangle,kappa]=anglemean(theta); 
%    
% anglestrength: can be thought of as the inverse variance. [varies between 0 and one]
% sigma: circular standard deviation
% confangle: a 95% confidence angle (confidence of the mean value)
% kappa: an estimate of kappa in the Von Mises distribution
%
% check: http://www.cosy.sbg.ac.at/~reini/huber_dutra_freitas_igarss_01.pdf
%
% Aslak Grinsted 2002
theta=mod(theta(:),2*pi);
n=length(theta);
S=sum(sin(theta));
C=sum(cos(theta));
meantheta=atan2(S,C);

if nargout<2
    return
end
%if ((S>C)&(C>0))
%    meantheta=atan(S/C);
%elseif (C<0)
%    meantheta=atan(S/C)+pi;
%else
%    meantheta=atan(S/C)+2*pi;
%end
%meantheta=mod(meantheta,2*pi);
Rsum=sqrt(S^2+C^2);
R=Rsum/n;

if (R<.53)
    kappa=2*R+R^3+5*R^5/6;
elseif (R<.85)
    kappa=-0.4+1.39*R+0.43/(1-R);
else
    kappa=1/(R^3-4*R^2+3*R);
end



%
%
% circular standard deviation:
sigma=sqrt(-2*log(R));


if nargout<4
    return
end


%conflim=.95;  
%a=(length(theta)-Rsum)/chi2pdf()

%this is true if the 
chi2=3.841; % = chi2inv(.95,1)
if ((R<.9)&(R>sqrt(chi2/(2*n))))
    confangle=acos(sqrt(2*n*(2*Rsum^2-n*chi2)/(4*n-chi2))/Rsum);
elseif (R>.9)
    confangle=acos(sqrt(n^2-(n^2-Rsum^2)*exp(chi2/n))/Rsum);
else %R is really really small ... 
    confangle=pi/2;
    warning('Confidence angle not well determined.')
    %this is not good, but not important because the confidence is so low anyway...
end
