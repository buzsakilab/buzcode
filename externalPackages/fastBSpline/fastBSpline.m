classdef fastBSpline
    %fastBSpline - A fast, lightweight class that implements 
    %non-uniform B splines of any order
    %
    %Matlab's spline functions are very general. This generality comes at
    %the price of speed. For large-scale applications, including model
    %fitting where some components of the model are defined in terms of
    %splines, such as generalized additive models, a faster solution is 
    %desirable.
    %
    %The fastBSpline class implements a lightweight set of B-spline
    %features, including evaluation, differentiation, and parameter fitting. 
    %The hard work is done by C code, resulting in up to 10x acceleration
    %for evaluating splines and up to 50x acceleration when evaluating 
    %of spline derivatives. 
    %
    %Nevertheless, fastBSplines are manipulated using an intuitive, high-
    %level object-oriented interface, thus allowing C-level performance 
    %without the messiness. Use CompileMexFiles to compile the required 
    %files. If mex files are not available, evaluation will be done in .m 
    %code, so you may still use the code if you can't use a compiler for your
    %platform.
    %
    %B splines are defined in terms of basis functions:
    %
    % y(x) = sum_i B_i(x,knots)*weights_i
    %
    %B (the basis) is defined in terms of knots, a non-decreasing sequence 
    %of values. Each basis function is a piecewise polynomial of order
    %length(knots)-length(weights)-1. The most commonly used B-spline is
    %the cubic B-spline. In that case there are 4 more knots than there
    %are weights. Another commonly used B-spline is the linear B-spline,
    %whose basis function are shaped like tents, and whose application
    %results in piecewise linear interpolation. 
    %
    %The class offers two static functions to fit the weights of a spline: 
    %lsqspline and pspline. It includes facilities for computing the basis 
    %B and the derivatives of the spline at all points.
    %
    %Constructor:
    %
    %sp = fastBSpline(knots,weights);
    %
    %Example use:
    %
    %%Fit a noisy measurement with a smoothness-penalized spline (p-spline)
    % x = (0:.5:10)';
    % y = sin(x*pi*.41-.9)+randn(size(x))*.2;
    % knots = [0,0,0,0:.5:10,10,10,10]; 
    %%Notice there are as many knots as observations
    % 
    %%Because there are so many knots, this is an exact interpolant
    % sp1 = fastBSpline.lsqspline(knots,3,x,y);
    %%Fit penalized on the smoothness of the spline
    % sp2 = fastBSpline.pspline(knots,3,x,y,.7);
    % 
    % clf;
    % rg = -2:.005:12;
    % plot(x,y,'o',rg,sp1.evalAt(rg),rg,sp2.evalAt(rg));
    % legend('measured','interpolant','smoothed');
    % 
    %fastBSpline properties:
    % outOfRange - Determines how the spline is extrapolated outside the
    %              range of the knots
    % knots      - The knots of the spline (read only)
    % weights    - The weights of the spline (read only)
    %
    %fastBSpline Methods:
    %  fastBSpline - Construct a B spline from weights at the knots
    %  lsqspline - Construct a least-squares spline from noisy measurements
    %  pspline   - Construct a smoothness-penalized spline from noisy
    %              measurements
    %  evalAt    - Evaluate a spline at the given points
    %  getBasis  - Get the values of the underlying basis at the given points
    %  Btimesy   - Evaluate the product getBasis(x)'*y
    %  dx        - Returns another fastBSpline object which computes the derivative of
    %              the original spline
    %
    %Disclaimer: fastBSpline is not meant to replace Matlab's spline functions; 
    %            it does not include any code from the Mathworks
    properties
        %knots - The vector of knots.
        knots = [];
        
        %weights - The vector of weights
        weights = [];
        
        %If outOfRange = fastBSpline.ZERO, outside the knots the spline = 0
        %   outOfRange = fastBSpline.CONSTANT, for x > knots(end), sp(x) =
        %   x(knots(end)), and similarly for x < knots(1)
        outOfRange = fastBSpline.ZERO;
        
        %usemex - Whether use mex acceleration is used
        usemex = false;
    end
    
    methods        
        function this = fastBSpline(knots,weights)
        %sp = fastBSpline(knots,weights)
        %
        %Class constructor for fastBSpline
        %The order of the spline is = length(knots)-length(weights) - 1
            this.knots = knots;
            this.weights = weights;     
            %Check that the required C files are compiled
            if exist('evalBSpline','file') == 3 && exist('evalBin','file') == 3
                this.usemex = true;
            else
                warning('fastBSpline:nomex','Mex files are not available, this will be a bit slower. Try running CompileMexFiles.');
            end
        end
        
        function s = dx(this)
        %spp = sp.dx
        %Given sp a spline, spp is another fastBSpline object such that 
        %spp(x) is the derivative of sp evaluated at x
            this.weights = [0;this.weights;0];
            this.knots = this.knots([1,1:end,end]');
            wp = protect(this.order*diff(this.weights) ./ ...
                 (this.knots(this.order + (2:length(this.weights)))-this.knots((2:length(this.weights)))));
            s = fastBSpline(this.knots(2:end-1),wp);
        end
        
        function Sx = evalAt(this,x)
        %Sx = sp.evalAt(x) - Evaluate the spline at x
        %x can be a vector or N-D matrix. size(Sx) == size(x)
        %This function is potentially accelerated through mex
            sz = size(x);
            x = x(:);
            
            order = this.order;
            if this.outOfRange == fastBSpline.CONSTANT
                x = max(min(x,this.knots(end)-1e-10),this.knots(1)+1e-10);
            end
            
            if this.usemex
                [~,thebins] = histc(x,this.knots);
                firstknot = max(min(thebins-order,length(this.weights)),1);
                 lastknot = min(thebins+1,length(this.weights)+1);

                Sx = evalBSpline(x,firstknot,lastknot,this.knots,this.weights,order);
            else
                %Fallback
                Sx = zeros(prod(sz),1);
                for ii = 1:length(this.knots)-order-1
                    v = x > this.knots(ii) & x <= this.knots(ii+order+1);
                    Sx(v) = Sx(v) + this.weights(ii)*bin(this.knots,ii,order,x(v));
                end
            end
            
            Sx = reshape(Sx,sz);
        end
        
        function Sx = getBasis(this,x)
        %B = sp.getBasis(x) - Gets the basis B sampled at points x
        %size(B) == [numel(x),length(sp.weights)]
            x = x(:);
            if this.outOfRange == fastBSpline.CONSTANT
                x = max(min(x,this.knots(end)-1e-10),this.knots(1)+1e-10);
            end
            order = this.order;
            if this.usemex
                [~,thebins] = histc(x,this.knots);
                firstknot = max(min(thebins-order,length(this.weights)),1);
                 lastknot = min(thebins+1,length(this.weights)+1);

                Sx = evalBin(x,firstknot,lastknot,this.knots,this.weights,order);
            else
                %Fallback
                Sx = zeros(numel(x),length(this.weights));
                for ii = 1:length(this.knots)-this.order-1
                    Sx(:,ii) = bin(this.knots,ii,this.order,x);
                end
            end
        end
        
        function d = Btimesy(this,x,y)
        %sp.Btimesy(x,y) - Computes the product sp.getBasis(x)'*y.
        %
        %This product is frequently needed in fitting models; for 
        %example, if E is an error function, r is the output of the 
        %spline evaluated at x, and y = dE/dr, then the derivative 
        %dE/dw is given by sp.Btimesy(x,y) by the chain rule. 
        %
        %This product can be computed efficiently in C without actually 
        %forming the B matrix in memory. 
            
            %Sanity check
            x = x(:);
            y = y(:);
            if length(x) ~= length(y)
                error('length(x) and length(y) must be equal');
            end
            if this.usemex
                [~,thebins] = histc(x,this.knots);
                firstknot = max(min(thebins-this.order,length(this.weights)),1);
                 lastknot = min(thebins+1,length(this.weights)+1);

                d = evalBinTimesY(x,firstknot,lastknot,this.knots,this.weights,this.order,y);
            else
                %Fallback
                d = this.getBasis(x)'*y;
            end
        end
        
        function x = order(this)
        %x = sp.order - Returns the order of the spline sp
            x = length(this.knots)-length(this.weights)-1;
        end
        
        function [this] = set.knots(this,k)
        %Set knots
            this.knots = k(:);
            if ~issorted(this.knots)
                error('Knots must be non-decreasing');
            end
            this.checkNum();
        end
        
        function [this] = set.weights(this,w)
        %Set weights
            this.weights = w(:);
            this.checkNum();
        end        
    end
    
    methods(Access=protected)
        function [] = checkNum(this)
            %Check whether number of weights is consistent with number of knots
            if length(this.knots) <= length(this.weights)
                error('length(knots) must be > than length(weights)');
            end
        end
    end
    
    methods(Static)
        function sp = lsqspline(knots,order,xi,yi)
        %BSpline.lsqspline(knots,order,xi,yi)
        %Fit the weights of the spline with given knots and order based 
        %on a least-squares fit of the data yi corresponding to xi. This 
        %function is static.
            xi = xi(:);
            yi = yi(:);
            sp = fastBSpline(knots,ones(length(knots)-order-1,1));
            sp.outOfRange = fastBSpline.CONSTANT;
            B = sp.getBasis(xi);
            w = B\yi;
            sp = fastBSpline(knots,w);
            sp.outOfRange = fastBSpline.CONSTANT;
        end
        
        function sp = pspline(knots,order,xi,yi,lambda)
        %BSpline.pspline(knots,order,xi,yi,lambda)
        %Like lsqspline but with a smoothness penalty of strength lambda 
        %on the weights of the spline. This function is static.
            xi = xi(:);
            yi = yi(:);
            sp = fastBSpline(knots,ones(length(knots)-order-1,1));
            sp.outOfRange = fastBSpline.CONSTANT;
            B = sp.getBasis(xi);
            e = ones(size(B,2),1);
            D = spdiags([-e 2*e -e], -1:1, size(B,2), size(B,2));
            w = [B;lambda*D]\[yi;zeros(size(B,2),1)];
            sp = fastBSpline(knots,w);
            sp.outOfRange = fastBSpline.CONSTANT;
        end
    end
        
    properties(Constant,Hidden)
        ZERO = 0;
        CONSTANT = 1;
    end
end

function x = protect(x)
    %Define 0/0 == 0
    x(isnan(x) | abs(x) == Inf) = 0;
end

function y = bin(knots,i,n,x)
    %bin - evaluate basis functions at given points. This is a fallback
    %when C files are not available
    delta = 1e-10;
    if n == 0
        y = x<knots(i+1) & x>= knots(i);
    elseif n == 1
        %Removes one degree of recursion
        y = (x-knots(i))/(delta + knots(i+1)-knots(i)).*(x<knots(i+1) & x>= knots(i)) + ...
            (knots(i+2) - x)/(delta + knots(i+2)-knots(i+1)).*(x<knots(i+2) & x>= knots(i+1));
    elseif n == 2

        m = (x<knots(i+2) & x>= knots(i+1));
        y = (x-knots(i))/(knots(i+2)-knots(i)+delta).*( ... 
                (x-knots(i))/(delta + knots(i+1)-knots(i)).*(x<knots(i+1) & x>= knots(i)) + ...
                (knots(i+2) - x)/(delta + knots(i+2)-knots(i+1)).*m) + ...
                ...
            (knots(i+3) - x)/(knots(i+3)-knots(i+1)+delta).*(...
                (x-knots(i+1))/(delta + knots(i+2)-knots(i+1)).*m + ...
                (knots(i+3) - x)/(delta + knots(i+3)-knots(i+2)).*(x<knots(i+3) & x>= knots(i+2)));
    else
        y = (x-knots(i))/(knots(i+n)-knots(i)+delta).*bin(knots,i,n-1,x) + ...
            (knots(i+n+1) - x)/(knots(i+n+1)-knots(i+1)+delta).*bin(knots,i+1,n-1,x);
    end
end
