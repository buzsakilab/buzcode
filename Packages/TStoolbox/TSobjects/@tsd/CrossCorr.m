function Cr = CrossCorr(tsa, tsb, binsize, nbins, varargin)

% Computes cross-correlation between TSD objects
%  
%  USAGE:
%  C = CrossCorr(tsa, tsb, binsize, nbins, options)
%  
%  INPUTS:
%  tsa 	- the triggering time series
%  tsb 	- the triggered time series 
%  binsize	- the size of the bins
%  nbins 	- the total number of bins in the returned cross-correlogram
%  
%  OPTIONS:
%  timeUnits 	  - the units for binsize, defaults to the units of tsa
%  fix_boundaries - if set to non-zero (default), will discard the points in the
%    		    triggering time series that are too close to the extremes of the
%  		    triggered time series, and would cause a spurious decline effect in the
%  		    cross-correlogram
%  errors 	  - if 'none' (default), only the cross-correlogram is returned in
%  		    the Data of the output tsd
%  		    if  'std', the standard deviation is returned as the second column of
%  		    the Data portion of the output
%  		    if 'sem', the standard error of the mean is returned instead

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
  
  defined_options = dictArray({ { 'timeUnits', {tsa.time_unit, {'Units', 'char'} } },...
		                { 'fix_boundaries', {1, {'numeric'} } } , ...
		                { 'errors', {'none', {'char'} } } , ...
		                { 'timeUnitsOut', {time_units('s'), {'Units', ...
		    'char'} } } });
  
  
                            
  
		   
  opt_varargin = varargin;
  getOpt;
  
  binsize = binsize * convert(time_units(timeUnits), tsa.time_unit);
  
  bs = 1 /  convert(tsa.time_unit, time_units(timeUnitsOut));
  
  
  
  if isa(tsb, 'tsd')
    T2 = tsb;
% $$$     if isempty(Range(tsa))
% $$$       nb = nbins;
% $$$       if mod(nb, 2) == 0
% $$$ 	nb = nb+1;
% $$$       end
% $$$       B = linspace((- binsize * floor(nbins / 2)):(- binsize * floor(nbins ...
% $$$ 						  / 2)), nb);
% $$$       B = B';
% $$$       C = NaN * ones(nb, 1);
% $$$       C2 = NaN * ones(nb, 1);
      
      
      
      
    switch errors
     case 'none'
      [C, B] = CrossCorr_c(Range(tsa), Range(T2, tsa.time_unit), binsize, ...
			   nbins, fix_boundaries);
      Cr = tsd(B, bs * C);
     case 'std'
      [C, B, C2] = CrossCorr_c(Range(tsa), Range(T2, tsa.time_unit), binsize, ...
			   nbins, fix_boundaries);
      Cr = tsd(B, bs * [C C2]);
     case 'sem'
      [C, B, C2] = CrossCorr_c(Range(tsa), Range(T2, tsa.time_unit), binsize, ...
			   nbins, fix_boundaries, 1);
      Cr = tsd(B, bs * [C C2]);

     otherwise
      error('Unrecognized option');
    end
  elseif isa(tsb, 'tsdArray')
    Cr = tsdArray;
    for i = 1:length(tsb)
      T2 = tsb{i};
      switch errors
       case 'none'
	[C, B] = CrossCorr_c(Range(tsa), Range(T2, tsa.time_unit), binsize, ...
			     nbins, fix_boundaries);
	Cr{i,1} = tsd(B, bs * C);
       case 'std'
	[C, B, C2] = CrossCorr_c(Range(tsa), Range(T2, tsa.time_unit), binsize, ...
				 nbins, fix_boundaries);
	Cr{i,1} = tsd(B, bs * [C C2]);
       case 'sem'
	[C, B, C2] = CrossCorr_c(Range(tsa), Range(T2, tsa.time_unit), binsize, ...
				 nbins, fix_boundaries, 1);
	Cr{i,1} = tsd(B, bs * [C C2]);

       otherwise
	error('Unrecognized option');
      end
    end
    
  else
    error('invalid value for tsb');
  end
  
      
    
      
  
  
  