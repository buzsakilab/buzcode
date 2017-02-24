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

transpX = 0;
if N1==0 && N2>0
    X = X';
    N1 = N2;
    N2=0;
    transpX = 1;
end

if min(size(X))>1
    
    l1 = size(X,1);
    l2 = size(X,2);
    Xf = flipud(X);
    
    if N2>0
        X0 = [[fliplr(Xf) Xf fliplr(Xf)]];
        X0 = [X0;[fliplr(X) X fliplr(X)]];
        X0 = [X0;[fliplr(Xf) Xf fliplr(Xf)]];
        X = X0;
        clear X0
    else
        X = [Xf;X;Xf];
    end
    clear Xf
    
    if N1~=0 | N2~=0
        
        %s.d. is N/(2*alpha) and we want s.d. = original N.
        %so here, with alpha=0.05, we multiply N by 10
        
        N1 = round(10*N1);
        N2 = round(10*N2);
        if N1~=0 | N2~=0
            if N1>0 & N2==0
                gw = gausswin(N1,5); %alpha = 0.05
            elseif N1>0 & N2>0
                gw = gausswin(N1,5)*gausswin(N2,5)';
            end
            if normBool
                gw = gw/sum(gw(:));
            end
            X = convn(X,gw,'same');
        end

    end
    if N2>0
        X = X(l1+1:2*l1,l2+1:2*l2);
    else
        X = X(l1+1:2*l1,:);
    end    

else
   if size(X,1)==1 
        X = X';
        transpX = 1;
   end
   if N1~=0
        l1 = size(X,1);
        X = [flipud(X);X;flipud(X)];
        
        %same as above
        N1 = round(10*N1);
        gw = gausswin(N1,5); %alpha = 0.05
        if normBool
            gw = gw/sum(gw);
        end
        X = convn(X,gw,'same');
        
        X = X(l1+1:2*l1);

   end
    
end

if transpX
    X = X';
end
