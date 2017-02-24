function f=doubleexponential_ap(params,x)
% f=doubleexponential_ap(params,x)
% f=params(1)*exp(params(2)*x)+params(3) + params(4)*exp(params(5)*x);

%f=params(1)*exp(params(2)*x) + params(3)*exp(params(4)*(x-params(5)));
f = params(1)*(x+0.01).^-1-params(2)*(x+0.01).^-2;