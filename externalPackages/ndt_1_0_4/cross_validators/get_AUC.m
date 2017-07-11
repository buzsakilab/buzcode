function [ROC_AUC_val roc_curve] = get_AUC(target_present_vals,  target_absent_vals)

%  This helper function calculate the area under an ROC curve based on a vector of values
%   when a target was present and a vector of values when a target was absent.  
%
% Input params:
%   target_present_vals:  column vector of values for the target present class
%   target_absent_vals:   column vector of values for the target absent class
%
% Outputs:
%   ROC_AUC_val:  The area under the ROC curve
%   roc_curve:  The ROC curve (defined here as the TP rate at each FP example)
%
% Notes:
%
%  1.  We are assuming that the target absent class has smaller values than the target present class.
%  
%  2.  The areas under the ROC curve is calculated from the the first negative
%       example  - i.e., the point where the FP rate is zero is not included.  This is done deliberately 
%       because the point at 0 FP can be ill defined, i.e., there can be multiple TP rates for a 0 FP rate.
%       This ambiguity at 0 FP could lead to slightly different ROC area values if more negative examples
%       were included, thus if the ROC area was calculated including a FP rate of 0
%       the property that the ROC area is invariant to the number of negative/positive test points
%       could potentially be violated (although the distortion would in practice probably be very small).

%==========================================================================

%     This code is part of the Neural Decoding Toolbox.
%     Copyright (C) 2011 by Ethan Meyers (emeyers@mit.edu)
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%==========================================================================   





% if row vectors are passed as arguments instead of column vectors, transpose them
if size(target_present_vals, 1) == 1
    target_present_vals = target_present_vals';
end
        
if size(target_absent_vals, 1) == 1
    target_absent_vals = target_absent_vals';
end


ROC_AUC_val = 0;


for iOrder = 1:2    % This is looped through twice to deal with tie values 

    
    if iOrder == 1
        all_vals = [target_present_vals; target_absent_vals];
        [tilda_junk, sorted_inds] = sort(all_vals, 'descend');
        all_indicators = [ones(size(target_present_vals)); zeros(size(target_absent_vals))];
    else
        all_vals = [target_absent_vals; target_present_vals];
        [tilda_junk, sorted_inds] = sort(all_vals, 'descend');
        all_indicators = [zeros(size(target_absent_vals)); ones(size(target_present_vals))];
    end


    sorted_indicators = all_indicators(sorted_inds);

    roc = cumsum(sorted_indicators);
    roc = roc((sorted_indicators == 0));

    roc = roc./length(target_present_vals);   % roc curve

    roc_curve_diff_orders(iOrder, :) = roc;   % save the roc curve


    
    % To calculate the area under the ROC I need to integrate over all FP_rate values.  Since each interval in the FP_rate range is the same, I can just calculate the mean TP_value between successive FP_points 
    %   and then afterward multiply by 1/num_target_absent points. If the TP_rate values at FP_rate(1) = A, FP_rate(2) = B, ... FP_rate(26) = Z, then we the AUC = ( (A+B./2) + (B+C./2) + ... + (Y+Z./2)) .* 1/num_negative points.
    %   But (A+B)./2 + (B+C)./2 + ... (Y+Z)./2 =  (A + B + B + C + C + D ... + X + Y + Y + Z)./2 =  (A + Z)./2  + (B + C + .. + Y), so AUC = ((A + Z)./2 + sum(B...Y))./num_target_absent_values.  

    % I am dividing by length(target_absent_vals) - 1, because there is one less interval compared to the number of target absent
    %  values since two target absent values define an interval and because I am not using the point where the FP rate is zero
    %  since this point is not well defined (and including it can lead to slight biases).
    
    %ROC_AUC_val_diff_orders(iOrder) = (sum(roc(2:end-1)) + (roc(1) + roc(end))./2)./ (length(target_absent_vals) - 1); 

    ROC_AUC_val = ROC_AUC_val  + ((sum(roc(2:end-1)) + (roc(1) + roc(end))./2)./ (length(target_absent_vals) - 1))./2;
    
end



 
%ROC_AUC_val = (ROC_AUC_val_diff_orders(1) + ROC_AUC_val_diff_orders(2))./2;   %mean(ROC_AUC_val_diff_orders);   % averaging to deal with ties

roc_curve = (roc_curve_diff_orders(1, :) + roc_curve_diff_orders(2, :))./2;   %mean(roc_curve_diff_orders);




