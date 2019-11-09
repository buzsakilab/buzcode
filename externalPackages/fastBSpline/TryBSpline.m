%% Load this function in Cell Mode and read the code for example usage
clear
%Fit a noisy measurement with a spline
x = (0:.5:10)';
y = sin(x*pi*.41-.9)+randn(size(x))*.2;
%Triple knots at the end
knots = [0,0,0,0:10,10,10,10];

%Least squares fit for the weights of a spline
s = fastBSpline.lsqspline(knots,3,x,y);

clf;
rg = -2:.005:12;
subplot(2,1,1);
plot(x,y,'o',rg,s.evalAt(rg),rg,sin(rg*pi*.41-.9));
legend('measured','approximated','actual');
subplot(2,1,2);
plot(rg,s.getBasis(rg)*diag(s.weights));


%%
%Fit a noisy measurement with a smoothness-penalized spline (p-spline)
x = (0.5:.5:9.5)';
y = sin(x*pi*.41-.9)+randn(size(x))*.2;
%double knots at the end
knots = [-1,-1,0:.5:10,11,11];

%Notice there are as many knots as observations
%Because there are so many knots, this is now an exact interpolant
s1 = fastBSpline.lsqspline(knots,3,x,y);

%Fit penalized on the smoothness of the spline
s2 = fastBSpline.pspline(knots,3,x,y,.7);

clf;
rg = -2:.005:12;
plot(x,y,'o',rg,s1.evalAt(rg),rg,s2.evalAt(rg));
legend('measured','interpolant','smoothed');

%%
%Define knots and weights directly and evaluate derivatives and so forth
knots = [-1,-1,0:10,11,11];
weights = sin((0:10)*pi*.41-.9);

%Create spline object
s = fastBSpline(knots,weights);

%first spatial derivative
sp = s.dx;

%second spatial derivative
spp = sp.dx;

%evaluate within this range
rg = -2:.02:12;
subplot(3,1,1);
plot(rg,s.evalAt(rg),rg,sp.evalAt(rg),rg,spp.evalAt(rg))
legend('Spline','Derivative','2nd Derivative');
title('Spline and its derivatives');

subplot(3,1,2);
plot(rg,s.getBasis(rg));
title('Basis functions underlying the spline');

subplot(3,1,3);
plot(rg,bsxfun(@times,s.getBasis(rg),s.weights'));
title('Weighted basis functions');

%%
%Now compare m-file vs C speed
%Only tested in Linux 64-bit but should work as is on other platforms
CompileMexFiles

%%
%Try without MEX acceleration
%For a cubic spline
knots = -2:102;
weights = randn(100,1);
x = sort(rand(1e5,1)*110)-5;
y = sort(rand(3e4,1)*110)-5;
z = sort(rand(1e5,1)*110)-5;
s = fastBSpline(knots,weights);
s.usemex = false;
fprintf('Without Mex acceleration:\n');
tic;
g1 = s.evalAt(x);
toc;
tic;
B1 = s.getBasis(y);
toc;
tic;
d1 = s.Btimesy(x,z);
toc;


%Try with MEX acceleration
s.usemex = true;
fprintf('With Mex acceleration:\n');
tic;
g2 = s.evalAt(x);
toc;
tic;
B2 = s.getBasis(y);
toc;
tic;
d2 = s.Btimesy(x,z);
toc;

fprintf('Relative numerical difference b/w m and C version:\n');
max(abs(g1-g2)/max(abs(g1)))
max(max(abs(B1-B2)))
max(abs(d1-d2)/max(abs(d1)))

%%
%For a linear spline
knots = -1:101;
weights = randn(100,1);
x = sort(rand(3e5,1)*110)-5;
y = sort(rand(1e5,1)*110)-5;
z = sort(rand(3e5,1)*110)-5;
s = fastBSpline(knots,weights);
s.usemex = false;
fprintf('Without Mex acceleration:\n');
tic;
g1 = s.evalAt(x);
toc;
tic;
B1 = s.getBasis(y);
toc;
tic;
d1 = s.Btimesy(x,z);
toc;


%Try with MEX acceleration
s.usemex = true;
fprintf('With Mex acceleration:\n');
tic;
g2 = s.evalAt(x);
toc;
tic;
B2 = s.getBasis(y);
toc;
tic;
d2 = s.Btimesy(x,z);
toc;

fprintf('Relative numerical difference b/w m and C version:\n');
max(abs(g1-g2)/max(abs(g1)))
max(max(abs(B1-B2)))
max(abs(d1-d2)/max(abs(d1)))