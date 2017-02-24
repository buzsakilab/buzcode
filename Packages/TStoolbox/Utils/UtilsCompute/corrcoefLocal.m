function C = corrcoefLocal(q,T,varargin)

% This function computes correlation coefficient getting rid of the global
% fluctuation. Inspired by Renart et al., Science, 2010
% Input:
%    q: matrix of binned spike trains
%    T: smoothing window. The spk trains will be convolved with a mexican
%    hat of T positive width and 4T negative width
% Optional:
% C = corrcoefLocal(q,T,'show')
% will show a wait bar.

showWb=0;
if nargin>2
    while length(varargin)
        if strcmp(varargin{1},'show')
            showWb=1;
            varargin = varargin{2:end};
            if length(varargin)>1
                varargin = varargin{2:end};
            else
                vararing = {};
            end
        elseif strcmp(varargin{1},'epoch')
            ep = varargin{2};
            if length(varargin)>2
                varargin = varargin{3:end};
            else
                vararing = {};
            end
        end
        
    end
end


t = [-20*T:0.2:20*T];
k = normpdf(t,0,T)-normpdf(t,0,T*sqrt(17)); %sqrt(T^2+J^2), J = 4T, so...

for ii=1:size(q,2)
   if showWb
        h = waitbar(ii/size(q,2));
   end
   q(:,ii) = conv(q(:,ii),k,'same');
    
end
if showWb
    close(h)
end
C = q'*q;
d = diag(C);
C = C./sqrt(d*d');