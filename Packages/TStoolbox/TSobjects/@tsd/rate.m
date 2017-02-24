function r = rate(tsa, varargin)

% TODO: overlap with intervalRate

% r = rate(tsa, is, TimeUnits) returns average rate of tsd object, measures in
% its timeSpan period, or in an arbitrary interval 
% 
% INPUTS:
% tsa: a tsd object
% an intervalSet , if present rate will be calculated on that range only  
% OUTPUTS:
% r: the rate   
% OPTIONS:
% TimeUnits: a units object or the abbreviation string representing the time
% units for rate calculation, rate will be returned in inversed time
% units (e.g. 's' => Hz) 
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


  iSet = [];
  if ~isempty(varargin)
    if isa(varargin{1}, 'intervalSet')
      iSet = varargin{1};
      if length(varargin) > 2
	varargin = varargin(2:end);
      else
	varargin = {};
      end
    end
  end
  
  opt_varargin = varargin;
  
  defined_options = dictArray( ...
  { { 'TimeUnits', {time_units('s'), {'char', 'units'} } } } );
  
  getOpt;
    
  TimeUnits = time_units(TimeUnits);
    
  
  if ~isempty(iSet)
    tsa = Restrict(tsa, iSet);
  end
  
  t = tsa.t;
  if ~isempty(t)

      if ~isempty(iSet)
        l = tot_length(iSet, TimeUnits);
      elseif length(t)>1
        cnvrt = convert(tsa.time_unit, TimeUnits);
        l = cnvrt * (t(end) - t(1));
      else
        l = eps;  
      end

      if l == 0 
        l = eps;
      end
      r = length(t) / l;
  else
      r = 0;
  end