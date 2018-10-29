function y = sigmoid(x,c,a)
% sigmoid evaluates a simple sigmoid function along x: 
% 
%         ______1________
%     y =        -a(x-c)
%          1 + e^
% 
%% Syntax 
% 
% y = sigmoid(x)
% y = sigmoid(x,c)
% y = sigmoid(x,c,a)
% 
%% Description
% 
% y = sigmoid(x) generates a sigmoid function along x. 
% 
% y = sigmoid(x,c) makes a sigmoid that scaled from zero to one, where
% c corresponds to the x value where y = 0.5. If c is not specified, a
% default value of c = 0 is assumed. 
% 
% y = sigmoid(x,c,a) specifies a, the rate of change. If a is close to
% zero, the sigmoid function will be gradual. If a is large, the sigmoid
% function will have a steep or sharp transition. If a is negative, the 
% sigmoid will go from 1 to zero. A default value of a = 1 is assumed if 
% a is not declared. 
% 
%% Example 1 
% A simple sigmoid: 
% 
% x = -10:.01:10;
% plot(x,sigmoid(x))
% 
%% Example 2: 
% Make a sigmoid function along x = 1 to 100, such that y(x=60) = 0.5: 
% 
% x = 1:100; 
% y = sigmoid(x,60);
% plot(x,y,'b','linewidth',2) 
% 
% Now do the same thing as above, but make the transition more gradual: 
% 
% y2 = sigmoid(x,60,0.1); 
% hold on
% plot(x,y2,'r','linewidth',2) 
% legend('default a = 1','a = 1/10') 
% legend boxoff 
% 
%% Author Info:
% Chad Greene, May 28, 2015.  
% http://www.chadagreene.com

%% Parse Inputs: 

narginchk(1,3) 
if nargin<3
    a = 1; 
else
    assert(isscalar(a)==1,'a must be a scalar.') 
end

if nargin<2
    c = 0; 
else
    assert(isscalar(c)==1,'c must be a scalar.') 
end

%% Perform mathematics: 

y = 1./(1 + exp(-a.*(x-c)));
