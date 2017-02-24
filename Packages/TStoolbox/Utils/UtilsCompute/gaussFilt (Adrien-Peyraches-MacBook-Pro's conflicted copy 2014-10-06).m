function X = gaussFilt(X,N1,varargin)

% Xsm = gaussFilt(X,N1,N2)
% 
% 1 or 2-d gaussian filtering with border correction
% Input:
%     X: vector (or matrix) to smooth
%     N1: sd of gausswin for the first dimension. If N1==0, no smoothing
%     N2 (optionnal): in case of X is a matrix, sd of smoothing on the 2nd dimension. If omitted, N2=N1
%     normalized (optionnal): if 'normalized' is true (default), the
%     gaussian windows are normalized to preserve the sum over all values
%     in the vector
% Output:
%     Xsm: smoothed matrix
%     
% Adrien Peyrache 2011

if length(size(X))>2
    error('Too many dimensions!')
end

N2=0;
normBool = 1;
if ~isempty(varargin)
        N2 = varargin{1};
        if length(varargin)==2
            normBool = varargin{2};
        end
    else
        N2=N1;
 end
 
 if min(size(X))>1
    
    l1 = size(X,1);
    l2 = size(X,2);
    
   
    
    Xf = flipud(X);
    X0 = [[fliplr(Xf) Xf fliplr(Xf)]];
    X0 = [X0;[fliplr(X) X fliplr(X)]];
    X0 = [X0;[fliplr(Xf) Xf fliplr(Xf)]];
    X = X0;
    clear X0 Xf
   
    if N1==N2 & N1>0
        mxWinSize = min(max(size(X)),round(5*N1));

        h = fspecial('gaussian', mxWinSize, N1);
        X = imfilter(X,h);

    elseif N1~=0 | N2~=0
        
        N1 = round(10*N1);
        N2 = round(10*N2);
        if N1~=0 | N2~=0
            if N1>0 & N2==0
                gw = gausswin(2*N1,5); %alpha = 0.05
            elseif N2>0 & N1==0
                gw = gausswin(2*N2,5)';
            elseif N1>0 & N2>0
                gw = gausswin(2*N1,5)*gausswin(2*N2,5)';
            end
            if normBool
                gw = gw/sum(gw(:));
            end
            X = convn(X,gw,'same');
        end

    end
     X = X(l1+1:2*l1,l2+1:2*l2);

else
    
    X = X(:);
   if N1~=0

        l1 = size(X,1);
        X = [flipud(X);X;flipud(X)];
        N1 = round(10*N1);
        gw = gausswin(N1,5); %alpha = 0.05
        if normBool
            disp('Norm')
            gw = gw/sum(gw);
        end
        X = convn(X,gw,'same');
        X = X(l1+1:2*l1);
   end

end
