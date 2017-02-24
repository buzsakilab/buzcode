function t = thresholdHysteresis(tsa, thr1, thr2, varargin)

%  Hysteresis thresholding 
%  
%  	USAGE:
%  	t = thresholdHysteresis(tsd, thr1, thr2, varargin) 
%  	
%  	returns the crossings of thr2 that are associated to crossing of thr1 
%  	
%  	INPUTS:
%  	tsa  - a tsd object 
%  	thr1 - the "high threshold"
%  	thr2 - the "low threshold"
%  	
%  	OUTPUT:
%  	t - the tsd of threshold crossing
%  	
%  	OPTIONS:
%  	'Crossing' - possible values are 'Rising' (default), and 'Falling'
%  	'Order'    - 'Before' look for crossings of thr2 that coem before thr1
%  		      'After' look for crossings of thr2 that come after thr1

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html



 opt_varargin = varargin;
 
 defined_options = dictArray({ { 'Crossing', { 'Rising' , {'char' } } },
		               { 'Order', { 'Before', { 'char' } } } });
 
 getOpt;
  
 t1 = threshold(tsa, thr1, 'Crossing', Crossing);
 
 t2 = threshold(tsa, thr2, 'Crossing', Crossing); 
 
 switch Order
  case 'Before'
   align = 'prev';
  case 'After'
   align = 'next';
  otherwise
   error('Unrecognized option value');
 end
 
 
 [tt, ix] = Restrict(t2, t1, 'align', align);
 
 ix = unique(ix);
 
 t = subset(t2, ix);
 
 
 
 
 