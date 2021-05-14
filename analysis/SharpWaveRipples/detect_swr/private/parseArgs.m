function params = parseArgs(params, input)
% This function takes in a structure, params, that defines a set of default
% parameter values, named according to its fieldnames. Input is either a
% structure with corresponding field name and user defined values or a set
% of property/value pairs, as formatted by varargin.

% 20-03-2014 JDL2

% parameter names for the calling routine
paramnames              = fieldnames(params);

if isstruct(input{1})
    % The user supplied parameters as a structure
    options   = input{1};
    usernames = fieldnames(options);
    [id, loc] = ismember(paramnames,usernames); %find param name matches
    inds      = find(id);
    % overwrite defaults with user input values
    for ii = inds'
        params.(paramnames{ii}) = options.(usernames{loc(ii)});
    end
elseif ischar(input{1})
    % The user supplied parameters as property/value pairs
    % count arguments
    nArgs = length(input);
    if round(nArgs/2)~=nArgs/2
        error('Input requires a struct of params or propertyName/propertyValue pairs.');
    end
    
    for pair = reshape(input,2,[]) % pair is {propName;propValue}
       inpName = pair{1};
       % Is this one of the relevant parameters?
       if any(ismember(paramnames,pair{1}))
          % overwrite deaults with user input values
          params.(inpName) = pair{2};
       end
       
    end
else
    error('Input requires a struct of params or propertyName/propertyValue pairs.');
end