function [pdf, sigma] = circ_ksdensity(T, ti, method, sigma)
%CIRC_KSDENSITY Probability density function of circular data.
%
% -------------
% INPUT
% -------------
% T - circular data values, n-by-1 vector
% ti - vector location in angles of even-spaced points where 
%     the pdf is to be computed (1:360 for example)
% method - {'msni'} | 'msn' | 'niqr' | 'mad' | 'std'
%     - type of method for infering the smoothing kernel width
%         from sample data statistics
%     - the methods differe in the statistics used for measuring 
%         the data spread and the choice of a multiplicative constant.
%     std - standard deviation of sample data
%           Use: good for normal & inappropriate for multimodal distributions
%           Ref: Silverman:45 (3.28)
%     niqr - normalized interquartile range
%           Use: for long-tailed distributions, but not multimodal
%           Ref: Silverman:47 (3.29), Wand:60
%           Note: Implemented in Matlab's ksdensity
%     mad - median absolute deviation
%           Use: for long-tailed distributions with possible outliers
%           Ref: Bowman:31
%           Note: Implemented in Matlab's ksdensity
%                 for data with censored values
%     msn - minimum of standard deviation and interquartile range
%           Use: compromise for skewed, multimodal distributions
%           Ref: Silverman:47 (3.30)
%     msni - msi improved
%           Use: as above, with less oversmoothing
%           Ref: Silverman:48 (3.31), Wand:62
% sigma - scalar giving the width of the smoothing kernel,
%     i.e. the sigma of the normal distribution used as kernel.
%     This argument is optional.
% 
% -------------
% OUTPUT
% -------------
% pdf - probaibility density function
% sigma - standard deviation of the smoothing kernel
% 
% -------------
% NOTES
% -------------
% 1. The kernel used in this function is a normal distribution 
% with an automatically computed optimal sigma as presented in:
% - Silverman B. W. (1998), Density Estimation for Statistics
% and Data Analysis, Chapman & Hall / CRC, Boca Raton (FL), 47-8.
% - Bowman Adrian W. & Adelchi Azzalini (1997) - Applied Smoothing 
% Techniques for Data Analysis, Oxford University Press, 31.
% - Wand M. P. & M. C. Jones (1995) - Kernel Smoothing, Chapman & 
% Hall, London, 60-3.
%
% 2. The difference with Matlab's ksdensity function
% is that this function is adaped to circular data, such as wind
% orientation. Using Matlab's function will give biased values 
% at the extremities of the pdf for circular data.
% 
% 3. The 'mad' method using the CircStat toolbox usually doesn't 
% accept datasets larger than about 1000 values, so in such
% cases a random selection of values is done on the input dataset.
% 
% -------------
% REQUIREMENTS & ACKNOWLEDGMENTS
% -------------
% The functions uses the CircStat toolbox for circular data
% by Philipp Berens. Download from:
% http://www.mathworks.com/matlabcentral/fileexchange/10676
%
% If this is however not available, it will use statistics
% for non-circular data to compute the smoothing kernel. 
% With the possible exception of msii method the results
% will then be in general slightly oversmoothed.
%
% -------------
% EXAMPLE
% -------------
% % Data: local orientations in degrees along the contour 
% % of a handwriting sample of the author; 
% % binned circular axial data between 0°-179°
% or_frq = [49,48,38,36,46,45,49,33,39,40,37,43,42,37,50,...
%     43,46,60,41,38,34,46,59,53,38,57,37,41,58,56,51,58,56,...
%     58,85,57,66,80,64,69,54,64,68,64,84,83,80,83,80,75,81,...
%     76,74,86,95,83,81,75,56,65,68,77,80,78,93,95,86,87,96,...
%     54,80,48,78,77,74,64,53,68,67,51,56,55,61,52,48,43,39,...
%     33,27,30,32,31,29,32,30,36,39,36,41,36,31,38,33,25,36,...
%     38,33,32,23,33,21,29,34,30,22,25,27,33,36,23,27,37,24,...
%     31,27,24,32,26,24,26,23,42,23,21,25,26,33,32,38,24,31,...
%     23,22,33,23,38,30,27,33,29,28,31,42,27,37,35,38,37,48,...
%     29,37,32,35,27,42,34,33,32,39,50,39,26,38,35,32,28,45,...
%     48,45,42;];
% 
% % convert frequency counts to measured data
% or_raw = [];
% for k = 1:180
%     or_raw = [or_raw; ones(or_frq(k),1)*(k-1)];
% end
% 
% % compute pdf
% [pdf.std,sigma.std] = circ_ksdensity(or_raw,0:179,'std');
% idx = floor(random('Uniform',1,length(or_raw),1000,1));
% [pdf.mad,sigma.mad] = circ_ksdensity(or_raw(idx),0:179,'mad');
% [pdf.niqr,sigma.niqr] = circ_ksdensity(or_raw,0:179,'niqr');
% [pdf.msn,sigma.msn] = circ_ksdensity(or_raw,0:179,'msn');
% [pdf.msni,sigma.msni] = circ_ksdensity(or_raw,0:179,'msni');
% [pdf.msni04,sigma.msni04] = circ_ksdensity(or_raw,0:179,'msni',...
%     sigma.msni/4);
% [pdf.ksdensity,ti,sigma.ksdensity] = ksdensity(or_raw',0:179);
% 
% % display results
% figure('name','Circular data pdf-s');
% 
% subplot(3,3,1)
% plot(0:179,or_frq)
% yd = (max(or_frq)-min(or_frq))/10;
% axis([-1 180 min(or_frq)-yd max(or_frq)+yd])
% title('frequency count')
% xlabel('Angle (degrees)')
% ylabel('Frequency')
% 
% % same scaling for all plots to allow comparison
% p = [pdf.std;pdf.mad;pdf.niqr;pdf.msn;pdf.msni; ...
%     pdf.msni04;(pdf.ksdensity)'];
% yul = max(p);
% yll = min(p);
% yul = yul + (yul - yll)/10;
% yll = yll - (yul - yll)/10;
% 
% subplot(3,3,2)
% plot(0:179,pdf.msni04)
% axis([-1 180 yll yul])
% title('msni/4')
% xlabel(['s = ' num2str(sigma.msni04)])
% ylabel('Probability')
% 
% subplot(3,3,3)
% plot(0:179,pdf.msni)
% axis([-1 180 yll yul])
% title('msni')
% xlabel(['s = ' num2str(sigma.msni)])
% 
% subplot(3,3,4)
% plot(0:179,pdf.std)
% axis([-1 180 yll yul])
% title('std')
% xlabel(['s = ' num2str(sigma.std)])
% 
% subplot(3,3,5)
% plot(0:179,pdf.niqr)
% axis([-1 180 yll yul])
% title('niqr')
% xlabel(['s = ' num2str(sigma.niqr)])
% 
% subplot(3,3,6)
% plot(0:179,pdf.msn)
% axis([-1 180 yll yul])
% title('msn')
% xlabel(['s = ' num2str(sigma.msn)])
% 
% subplot(3,3,7)
% plot(0:179,pdf.mad)
% axis([-1 180 yll yul])
% title('mad')
% xlabel(['s = ' num2str(sigma.mad)])
% 
% subplot(3,3,9)
% plot(0:179,pdf.ksdensity)
% axis([-1 180 yll yul])
% title('ksdensity')
% xlabel(['s = ' num2str(sigma.ksdensity)])
%
% -------------
% CREDITS
% -------------
% Vlad Atanasiu
% atanasiu@alum.mit.edu | http://www.waqwaq.info/atanasiu/
% 2010.12.20
%     - creation

% compute frequency data
Tf = hist(T,ti);
Tf = Tf';
tn = length(ti);

% padd data
t2 = floor(tn/2);
Tp = [Tf(t2+1:tn); Tf; Tf(1:t2)];

% compute kernel width
if nargin < 4
    if nargin < 3
        method = 'msni';
    end
    % checks for CircStat toolbox availability
    if exist('circ_std','file') > 0
        T = circ_axial(circ_ang2rad(T));
        k = (4/3)^0.2;
        switch method
            case 'msni'
                k = ((280*pi^0.5)/729)^0.2;
                s = min(circ_std(T), circ_iqr(T)/1.349);
            case 'msn'
                s = min(circ_std(T), circ_iqr(T)/1.349);
             case 'niqr'
                s = iqr(T)/1.349;
            case 'mad'
                if length(T) > 1000
                    idx = floor(random('Uniform',...
                        1,length(T),1000,1));
                    T = T(idx);
                end
                s = circ_median(abs(T - circ_median(T)))/0.6745;
            otherwise
                s = circ_std(T);
        end
        sigma = k*s*length(T)^-0.2;
        sigma = sigma*180/pi;
    else
        k = (4/3)^0.2;
        switch method
            case 'msni'
                k = ((280*pi^0.5)/729)^0.2;
                s = min(std(T), iqr(T)/1.349);
            case 'msn'
                s = min(std(T), iqr(T)/1.349);
             case 'niqr'
                s = iqr(T)/1.349;
            case 'mad'
                s = median(abs(T - median(T)))/0.6745;
            otherwise
                s = std(T);
        end
        sigma = k*s*length(T)^-0.2;
    end
end

% make Gaussian kernel
ti = 1:length(Tp);
mu = length(ti)/2 + 1;
G = (exp( -((ti - mu).^2)./ (2*sigma^2) ))';

% convolve kernel with data
pdf = ifft(fft(Tp).*fft(G));

% remove padding
pdf = pdf(tn-t2+1:tn-t2+tn);

% normalize values
pdf = pdf./sum(pdf);

