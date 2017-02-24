% this is a script, not a function, so it can assecc the function
% workspace. it will process the array opt_varargin looking for option
% name/option value pairs. 
% it expects to find a defined_options dict array, with option names as
% keys. The value is a cell array first element contains the default
% value, the second is a cell array of admissible types for the value
% for example. Empty value mneas don't check type. 
% defined_options = 
% { length: { 120, {'numeric'}}, 
%   mode: { 'none', {'char'}}, 
%  units: { time_units('ts'), {'char', 'units'} },
%   }


% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


getOpt_specified_options = dictArray;

% parse the argument list
while length(opt_varargin(:)) > 0
  getOpt_opt = opt_varargin{1};
  getOpt_value = opt_varargin{2};
  if length(opt_varargin) > 3
    opt_varargin = opt_varargin(3:end);
  else
    opt_varargin = {};
  end

  
  
  if ~has_key(defined_options, getOpt_opt)
    error(['Unrecognized option ' getOpt_opt]);
  end
  
  do = defined_options{getOpt_opt};
  getOpt_opt_types = do{2};
  
  if ~isempty(getOpt_opt_types)
    getOpt_good_type = 0;
    
    for getOpt_tc = getOpt_opt_types
      getOpt_t = getOpt_tc{1};
      if isa(getOpt_value, getOpt_t)
	getOpt_good_type = 1;
	break
      end
    end
    if ~getOpt_good_type
      error([ 'Option value for option ', getOpt_opt, ...
	      ' of unproper type']);
    end
  end
  
  if isa(getOpt_value, 'char')
    getOpt_value = strrep(getOpt_value, '''', '''''');
    eval([getOpt_opt, ' =  ''', getOpt_value ''' ; ']);
  else
    eval([getOpt_opt ' =  getOpt_value ;']);
  end
  getOpt_specified_options{getOpt_opt} = 1;

end


% now set the default arguments for options that were not set

for getOpt_kc = (keys(defined_options))'
  getOpt_k = getOpt_kc{1};
  if ~has_key(getOpt_specified_options, getOpt_k)
    do = defined_options{getOpt_k};
    getOpt_default = do{1};
    eval([getOpt_k ' =  getOpt_default ;']);
  end
end


clear getOpt_*
  
  
  
  
  