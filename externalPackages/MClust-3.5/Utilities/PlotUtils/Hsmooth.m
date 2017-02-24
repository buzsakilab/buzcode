function [S, bb] = Hsmooth(H)
%
%  [S, bb] = Hsmooth(H)
%
% smoothes a 2D histogram H (= 2D array)
% with a 6-th order double low pass firls bb (linear phase least-square FIR filter)
% by lipa

% create filter
if ~isempty(which('firls'))
	b = firls(6, [0 .5 .5 1], [1 .5 .5 0]);  
	bb = kron(b',b);    % 2D filter = tensor product of 1D filters
	
	S = filter2(bb,H);  % first pass (introduces a linear phase shift)
	S = filter2(bb',S);  % second pass (compensates phase shift)
else
	warning('MCLUST:ToolboxUnavailable', 'Skipping contour smoothing. Signal processing toolbox unavailable.');
	S = H;
end