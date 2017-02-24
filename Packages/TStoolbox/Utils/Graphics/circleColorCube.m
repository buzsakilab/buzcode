function circleColor(M, x, y, cmap)

    if nargin < 5
        Mcol = M;
    end
    
    if length(x) ~= size(M, 2) | length(y) ~= size(M,1)
        error('x and y should have the same length as dimensions of M');
    end
    
    

    n = 100;
    l = linspace(0, 2*pi, n);
    th = linspace(0, 2*pi, n);
    r = 1;
    rh = ones(1, n)*r;
    [xx, yy] = pol2cart(th, rh);

    Mmax = max(M(:));
    Mmax = 1;
    M = ceil(64*M / Mmax);
    
    max_r = min(min(diff(x)), min(diff(y))) / (64*2);
    
    %cc = resample(hsv, length(x), 64);
    zz = linspace(-1,1,length(y));
    for i = 1:length(y)
        cc = brighten(hsv, zz(i));
        cq(i,:,:) = resample(cc, length(x),64);
    end
    
    cq(cq > 1) = 1;
    cq(cq < 0) = 0;
    
    
    
    
    max_r = 1/128;
    colormap(cmap)
    for i = 1:length(y)
        for j = 1:length(x)
            [xx, yy] = pol2cart(th, rh * M(i,j) * max_r);
            xx = xx + x(j);
            yy = yy + y(i);
            patch(xx, yy, repmat(cq(i,j,:), [1 n 1]) );
            hold on 
        end
    end
    