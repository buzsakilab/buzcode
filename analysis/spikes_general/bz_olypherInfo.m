function [track_info,pos_info_val] = bz_olypherInfo(data,round_to,smoothing)
% USAGE
%        [track_info,pos_info_val] = bz_olypherInfo(data,round_to,smoothing)
%
% INPUTS
%         data - matrix (M x N x D) where M is the number of cells to analyze,
%                N is the number of trials for each cell, and D is the number
%                of time bins
%         round_to - integer value that data is discritized to. A value of
%                    2 means all data will be rounded to nearest 2 (i.e.
%                    2,4,8...)
%         smoothing - 0 no smoothing, else smooth with N bins
% OUTPUTS
%         track_info - matrix (M x D-2) of information scores across all 
%                      trials or behavior windows (N)
%         pos_info_val - matrix (M x N x D) of all information values 
%                        that are calculated  
%
% this function calculates the information carried in the firing rate of single
% neurons per spatial/temporal bin
%
%  Written by David Tingley
%  UCSD Cognitive Neuroscience
%  1/15/12

%% TODO
% - convert to varargin with input parser
% - add 'exclude' input to remove 0's from info calculation
% - 
          
M = size(data,1);
N = size(data,2); 
D = size(data,3);  

if M == 0
    M = 1
end
if N == 0
    N = 1
end
if D == 0
    D = 1
end
pos_info_val = zeros(M,N,D);
a = N*D;

%% Rounding 
if smoothing ~= 0
    for i = 1 : M
        for k = 1:N
            data(i,k,:) = smooth(squeeze(data(i,k,:)),smoothing).*smoothing;
        end
    end
end
data = round(data./round_to)*round_to;

%% Info Analysis
     
for i = 1 : M
      for x = 1 : D   
         for k = 1:N

                q = data(i,k,x);
 
                pKx = ((length(find(data(i,:,x) == q)))/N);

                pK = ((length(find(data(i,:,:) == q)))/a);

                if pK == 0 || pKx == 0 || pKx < pK
                    pos_info_val(i,k,x) = pos_info_val(i,k,x);
                else
                    pos_info_val(i,k,x) = pos_info_val(i,k,x) + (pKx*log2(pKx/pK));
                end
            
         end        
     end

end

for i = 1:M
track_info(i,:) = sum(pos_info_val(i,:,2:end-1),2);
end





