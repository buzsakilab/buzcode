function [R, il] = Restrict(tsa, varargin)

%  Restrict TSD to a new timestamp
%  
%  USAGE:
% R = Restrict(tsa, varargin) 
%
%  There are two main ways to invoke this function: in the first an array
%  of timestamps is specified, and the tsd made out of the points in tsa that
%  approximated those timestamps is returned 
%  for these usage call like this
%  R = Restrict(tsa, k, 'OptionName', 'OptionValue'), where k may be a
%  numerical array, or a ts, tsd object
%  The option 'align' specifies the type of approximation to
%  produce. possible values are
%  'prev': returns the points that immediately preceded each of the
%  timestamps
%  'next': returns the points that immediately followed each of the
%  timestamps
%  'closest': returns the points that best approximated the timestamps
%  'equal': returns the points that had the same time as the
%  timestamps. If perfect alignment cannot be achieved for all the points,
%  it returns an error.
%  
%  In the other use, a series of intervals is specified, and the
%  Restrict-ed tsd returns the datapoints included in those intervals. For
%  this use call like:
%  R = Restrict(tsa, t0, t1, 'OptionName', 'OptionValue'), with t0, t1
%  start and end of restrict intervals or
%  R = Restrict(tsa, is, 'OptionName', 'OptionValue')
%  with is intervalSet object containing the restriciting intervals
%  [R, ix] = Restrict(tsa, ...) returns, as second output, the indices of
%  the selected points  
  
  
% the current version is an evolution of A. David Redish's code
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% released under the GPL   
  

% persistent mex_count;

% $$$ if isempty(mex_count)
% $$$   mex_count = 0;
% $$$ end
% $$$ 
  
restrict_type = '';
align_type = 'closest';

while length(varargin) > 0
  v = varargin{1};
  varargin = varargin(2:end);
  
  if isa(v, 'numeric') | isa(v, 'tsd')
    if ~isempty(restrict_type)
      error('Can only specify align or intervals');
    end
    st = v;
    time_type = 0;
    if length(varargin) > 0
        if isa(varargin{1}, 'numeric') | isa(varargin{1}, 'tsd') % then it is interbals
            v = varargin{1};
            varargin = varargin(2:end);
            en = v;
            restrict_type = 'intervals';
        else
            restrict_type = 'align';
        end
    else
        restrict_type = 'align';
    end

    
  elseif isa(v, 'intervalSet') %intervals expressed in intervalSet
    if ~isempty(restrict_type)
      error('Can only specify align or intervals');
    end
    st = v;
    restrict_type = 'iset';
  elseif isa(v, 'char') % read in option
    if length(varargin) > 0
      ov = varargin{1};
      varargin = varargin(2:end);
    else
      error('must specify value for option');
    end
    
    time_type = 0;
    
    switch(v)
        case 'align'
            if ~isa(ov, 'char')
                error('value for option align must be string')
            end
            switch(ov)
                case {'prev', 'next', 'closest', 'equal'} % option will be ignored for int
                    align_type = ov;
                otherwise
                    error('Unrecognized option value');
            end
        case 'time'
            switch(ov)
                case 'original'
                    time_type = 0;
                case 'align'
                    time_type = 1;
            end
            % put more options here
        otherwise
            error('Unrecognized option');
    end
  end
end

if isempty(restrict_type)
  error('Must specify arguments to Restrict!');
end

if strcmp(restrict_type, 'intervals')
  st = intervalSet(st, en);
  restrict_type = 'iset';
end


switch(restrict_type)
 case 'iset'
% $$$   mex_count = mex_count+1;
% $$$   display(mex_count);
% $$$   if mex_count == 456 
% $$$     42;
% $$$     keyboard;
% $$$   end
  

  if ~isempty(Start(st, tsa.time_unit))
    ix = Restrict_idx_iSet(Range(tsa, tsa.time_unit), ...
			   Start(st, tsa.time_unit), ...
			   End(st, tsa.time_unit));
  else 
    ix = zeros(0,1);
  end
  
    
 case 'align'
  if isa(st, 'tsd');
    st = Range(st, 'ts');
  end
  
  if (~isempty(st)) & length(tsa) > 0 
    ix = Restrict_idx_align(Range(tsa, 'ts'), st, align_type);
  else
    ix = zeros(0,1);
  end
  
end

tt = tsa.t(ix);
if strcmp(restrict_type, 'align') & time_type == 1
    tt = st;
end


if isa(tsa, 'ts')
  R = ts(tt);
else
  R = tsd(tt, SelectAlongFirstDimension(tsa.data, ix));
end

if nargout == 2
  il = ix;
end
