function O = dictArray(varargin)
% O = dictArray(inp) constructor for the dictArray class, inspired to
% Python dict class. 
%
% Input could either be another dictArray object (copy constructor) or a
% cell array like the following
% { { Key1, Value1}, 
%   { Key2, Value2}, 
%   Key3,
%   ... 
%   { KeyN, ValueN} } 
% If the element in the cell array is a key, with no value associated,
% the value will be empty
% OPTIONS:
% 'KeysReadOnly': if set to nonzero, the keys will not be modifiable
% after construction
% named dictArray to avoid conflicts with ADR's class dict  
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  

   keys = {};
   values = {};
   KeysReadOnly = 0;
   
  
  if length(varargin) > 0
    
    inp = varargin{1};

     
    
    
    opt_varargin = varargin(2:end);

    KeysReadOnly = 0;
    if ~isempty(opt_varargin)
      if ~isa(opt_varargin{1}, 'char')
	error('option name must be string');
      end
      switch opt_varargin{1}
       case 'KeysReadOnly'
	KeysReadOnly = opt_varargin{2};
       otherwise 
	error('unrecogniezed option');
      end
    end
    
      
    if isa(inp, 'dictArray')
      O = inp;
      O.keysReadOnly = KeysReadOnly;
      return
    end
    
    if ~isa(inp, 'cell')
      error('input must be cell array of keys or key/value pairs');
    end
    
    inp = inp(:);
    
    
    % check validity of input pairs
    
    for i = 1:length(inp)
      pair = inp{i};
      
      if ~( (iscell(pair) & length(pair) == 2) | ...
	    (isa(pair, 'char') ) )
	error('each cell in input must be a key or a key/value pair');
      end
      
      if length(pair) == 2
	if ~isa(pair{1}, 'char')
	  error('keys must be string');
	end
      else
	if ~isa(pair, 'char')
	  error('keys must be string');
	end
      end
      
      
      
    end
    
    keys = cell(length(inp), 1);
    values = cell(length(inp), 1);
    
    for i = 1:length(inp)
      pair = inp{i};
      
      if length(pair) == 2
	keys{i} = pair{1};
	values{i} = pair{2};
      else
	keys{i} = pair;
	values{i} = [];
      end
      
    end
    
    
    
    
  else % constructor was invoked with no argument, make empty dictArray
    
    keys = {};
    values = {};
    
  end
  
  O.keys = keys;
  O.values = values;
  O.keysReadOnly = KeysReadOnly;
  O = class(O, 'dictArray');
  
  
  
  
  
    
      
    
  
  