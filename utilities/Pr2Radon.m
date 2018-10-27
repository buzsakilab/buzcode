function [slope integral] = Pr2Radon(varargin)
% computes max projection line, using radon transform, on Pr matrix

if nargin == 1
    Pr = varargin{1};
    plotting = 0;
elseif nargin == 2
    Pr = varargin{1};
    plotting = varargin{2};
end

Pr(isnan(Pr)) = 0;
theta = 0:.5:180;%0:.5:180;
[R xp] = radon(Pr,theta);
y(1) = floor((size(Pr,1)+1)./2);
x(1) = floor((size(Pr,2)+1)./2);


[Y,I] = max(R);
% [a b]  = max(makeLength(Y,length(Y)*10));
% [pks locs] = findpeaks(Y);
% [pks locs] = max(Y);
locs = 1:length(Y);
% search all peaks for max integral
for pk = 1:length(locs)
%     [a b]  = max(Y);
%     b = locs(pk);
    angle = theta(locs(pk));
    % offset = x(1) - (I(b)) + 1;
    % xpp = makeLength(xp,length(xp)*10);
    % II = makeLength(I,length(I)*10);
%     offset = xp(I(b));
%     [a b] = max(R(:,locs(pk)));
    offset = xp(I(locs(pk)));
    if offset == 0
        offset = .01;
    end

    y(2) = y(1) + offset*sin(deg2rad(-angle));
    x(2) = x(1) + offset*cos(deg2rad(-angle));
    coeffs = fPolyFit(x, y, 1);
    xx = 1:size(Pr,2);
    yy = (-1/coeffs(1))*(xx - x(1)) + y(1) - offset;

    coeffs = fPolyFit(xx, yy, 1);
    slope(pk) = coeffs(1);

    if abs(slope(pk)) < 2 * size(Pr,1) ./ size(Pr,2) & abs(slope(pk)) > 1.5 % rise/run limit to calc integral (must be in the frame)
        for i=1:length(xx)
%             inds = Restrict(xx(i)-5:xx(i)+5,[1 size(Pr,2)]);
            if yy(i) > .5 & yy(i) < size(Pr,1) + .5
                curve(i) = Pr(round(yy(i)),xx(i));
%                 curve_fit(i) = sum(Pr(round(yy(i)),inds));
            else
                curve(i) = nan;
%                 curve_fit(i) = nan;
            end
        end
        integral(pk) = nanmean(curve); clear curve;
%         integral(pk) = nanmean(curve_fit); clear curve_fit
    else
        integral(pk) = NaN; clear curve;
        slope(pk) = NaN;
    end
    
end
integral = double(integral); % weird typecasting fix
[integral idx] = max(integral);
slope = slope(idx);

if plotting
    offset = xp(I(locs(idx)));
    if offset == 0
        offset = .01;
    end
    angle = theta(locs(idx));
    y(2) = y(1) + offset*sin(deg2rad(-angle));
    x(2) = x(1) + offset*cos(deg2rad(-angle));
    coeffs = fPolyFit(x, y, 1);
    xx = 1:size(Pr,2);
    yy = (-1/coeffs(1))*(xx - x(1)) + y(1) - offset;
    
    
%     figure(10)
    cla
    imagesc(Pr)
    hold on
    line([x(1) x(2)],[y(1) y(2)])
    plot(xx,yy,'r')
%     plot(xx,yy-5,'w')
%     plot(xx,yy+5,'w')
    title(integral)
%     pause
end

% iter = 0;
% while integral < 1 & ~isnan(integral) % this may be a bad idea
%     [a b] = min(abs((theta-angle)));
%     R(find(xp==offset),b) = 0;
%     y(1) = round(size(Pr,1)./2);
%     x(1) = round(size(Pr,2)./2);
% 
% 
%     [Y,I] = max(R);
%     [a b]  = max(Y);
%     angle = theta(b);
%     % offset = x(1) - (I(b)) + 1;
%     offset = xp(I(b));
% 
% 
%     y(2) = y(1) + offset*sin(deg2rad(-angle));
%     x(2) = x(1) + offset*cos(deg2rad(-angle));
%     coeffs = polyfit(x, y, 1);
%     xx = 1:size(Pr,2);
%     yy = (-1/coeffs(1))*(xx - x(1)) + y(1) - offset;
% 
%     coeffs = polyfit(xx, yy, 1);
%     slope = coeffs(1);
% 
%     if abs(slope) < size(Pr,1) ./ size(Pr,2)% rise/run limit to calc integral (must be in the frame)
%         for i=1:length(xx)
%             if yy(i) > .5 & yy(i) < size(Pr,1) + .5
%             curve(i) = Pr(round(yy(i)),xx(i));
%             end
%         end
%         integral = nansum(curve);
%     else
%         integral = NaN;
%         slope = NaN;
%     end
%     integral = double(integral); % weird typecasting fix
%     iter = 1 + iter
%     if iter > 10
%        disp('and youre out of tries...')
%        integral = NaN;
%        slope = NaN;
%     end
% end

