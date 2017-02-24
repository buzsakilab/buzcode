function [tmd] = triggeredMeanData(tsa, events, tmax, varargin)

%  Event-triggered mean  of the Data component, triggered on events (TODO)
%  
% TODO: Overlap wth the ETAverage function. Need cleaning
%  
%  USAGE:
%   [tmd] = triggeredMeanData(tsa, events, tmax, options)
%   
% INPUTS:
% tsa: a tsd object, assumed to have constant step timestamps 
% events: a ts(d) object containing events to trigger the peth
% tmax: the resulting peth will run from -tmax to tmax, with a step equal
% to the step of tsa.
% OUTPUTS:
% tmd: a tsd containing the time triggered average, in the same time
% units as tsa
% OPTIONS:
% timeUnits: the units for tmax 
% error: if set to 'sdev', will return standard deviation as the second,
%        if set to 'sem', will return standard errors of the mean 
% column in the output tsd
 
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
  
  
  defined_options = dictArray({{  'timeUnits', { 'ts', {'Units', 'char'} } }, ...
		               { 'errors', {'none', {'char'} } }});
  
  opt_varargin = varargin;
  getOpt;
  
  timeUnits = time_units(timeUnits);
  
  tmax = tmax * convert(timeUnits, tsa.time_unit);
  
  t = Range(tsa);
  step = median(diff(t));
  
  
  sz = ceil(tmax/step);

  [R, ev] = Restrict(tsa, events);
  
  if(size(ev, 2) ~= 1)
    error('ev should be column');
  end
  
  
  B = -sz:sz;
  dd = Data(tsa);
  
  if size(dd,2) ~= 1
    error('defined only for 1 column data');
  end
  
  dd = [dd ; NaN];
  
    evtot = ev;
  C= zeros(size(B));
  C2= zeros(size(B));
  norm = zeros(size(B));
  
  dstep = 10000;
  
  for i = 1:dstep:length(evtot)
    ev = evtot(i:min(length(evtot), i+dstep-1));
    tms = repmat(ev, 1, length(B))+repmat(B, length(ev), 1);
    t1 = length(t);
    tms(find(tms < 1)) = t1+1;
    tms(find(tms > t1)) = t1 + 1;
    C = C + nansum(dd(tms));
    C2 = C2 + nansum(dd(tms).^2);
    norm = norm + sum(isfinite(dd(tms)));
  end
    
  lt = norm;
  C = C ./ lt;
  C2 = C2./(lt-1) - (lt./(lt-1)).*C.*C;
  
  C2 = sqrt(C2);

  B = B * step;
  
  switch errors
   case 'sdev'
    
    tmd = tsd(B', [C', C2'], 'TimeUnits', tsa.time_unit);
   case 'sem'
        tmd = tsd(B', [C', (C2./ sqrt(norm))'], 'TimeUnits', tsa.time_unit);
   case 'none'
    tmd = tsd(B', C', 'TimeUnits', tsa.time_unit);
  end
  