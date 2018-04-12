function mixtures = generateMixtures(targets,interferers,varargin)
%GENERATEMIXTURES Generate arrays of mixtures from targets and interferers.
% 
%   MIXTURES = IOSR.BSS.GENERATEMIXTURES(TARGETS,INTERFERERS) generates an
%   array of MIXTURE objects using SOURCE objects TARGETS and INTERFERERS.
%   Each mixture object contains one target element and SIZE(INTERFERERS,2)
%   interferers (i.e. as many interferers as there are columns). How the
%   inputs are combined is determined by the nature of the input or the
%   'COMBINE' setting (see below).
% 
%   MIXTURES = IOSR.BSS.GENERATEMIXTURES(...,'PARAMETER',VALUE) allows
%   additional options to be specified. The options are listed below ({}
%   indicate defaults).
% 
%     Independent variables
% 
%       'azimuths'      : {zeros([1 size(interferers,2)+1])} | array
%           Specify the azimuths for the sources. The numeric array should
%           have one more column than INTERFERERS. The first column is the
%           target azimuth; subsequent columns are assigned to the
%           interferer sources such that the Nth azimuth is assigned to the
%           N-1th interferer column (N>1). Specifying this parameter will
%           overwrite the source property.
%       'elevations'    : {zeros([1 size(interferers,2)+1])} | array
%           Specify the elevations for the sources. The specification is
%           the same as for azimuths.
%       'sofa_paths'         : {[]} | str | cellstr
%           Specify the SOFA file(s) for the mixtures as one or more paths
%           to SOFA files that are convolved with sources. The parameter
%           should be a character array or cell array of strings.
%       'tirs'          : {0} | scalar
%           Specify the target-to-interferer ratios for the mixtures.
% 
%       TARGETS and INTERFERERS are also considered independent variables.
% 
%     Settings
% 
%       'cache'         : {false} | true
%           Specify whether mixtures are cached as audio files, rather than
%           being rendered ono-the-fly. The setting invokes
%           mixture.write(). Files are named 'mixture-NNNN.wav' where NNNN
%           is a counter padded with leading zeros. Target and interferer
%           signals are also cached.
%       'combine'       : 'all' | 'rows'
%           Determines how the inputs are combined. If the independent
%           variables have the same number of rows (ignoring empty or
%           scalar variables) then the variables are combined on a
%           row-by-row basis, with scalar variables applied to each row.
%           This is the 'rows' option, which is the default when the above
%           conditions apply; the conditions must be met in order to
%           combine variables in this way. Alternatively, variables are
%           combined in all combinations. This is the 'all' option, which
%           is the default when variables have rows of different lengths;
%           variables may each have any number of rows.
%       'decomposition' : {'stft'} | 'gammatone'
%           The time-frequency decomposition used by the mixture.
%       'folder'        : {'mixture_temp'} | str
%           Specify a folder for storing cached audio files.
%       'fs'            : {16000} | scalar
%           Sampling frequency for the mixtures.
% 
%   See also IOSR.BSS.MIXTURE, IOSR.BSS.SOURCE, SOFALOAD.

%   Copyright 2016 University of Surrey.

    %% get input
    
    IVs = struct(...
        'azimuths',zeros([1 size(interferers,2)+1]),...
        'elevations',zeros([1 size(interferers,2)+1]),...
        'sofa_paths',[],...
        'tirs',0);
    
    settings = struct(...
        'fs',16000,...
        'cache',false,...
        'combine',[],...
        'folder','mixture_temp',...
        'decomposition','stft',...
        'stft',[],...
        'gammatone',[]);
    
    validOptions = [fieldnames(IVs); fieldnames(settings)];
    
    %% check input
    
    % overwrite the default settings
    IVs = overwrite(IVs,validOptions,varargin);
    settings = overwrite(settings,validOptions,varargin);
    
    % decomposition settings
    if isempty(settings.stft)
        settings.stft = struct('win',1024','hop',512);
    end
    if isempty(settings.gammatone)
        settings.gammatone = struct('cfs',iosr.auditory.makeErbCFs(20, settings.fs/2, 64),...
            'frame',round(0.01*settings.fs));
    end
    
    % check IVs
    assert(isa(targets,'iosr.bss.source'), 'iosr:generateMixtures:invalidTargets', '''TARGETS'' must be of type iosr.bss.source')
    assert(isa(interferers,'iosr.bss.source'), 'iosr:generateMixtures:invalidInterferers', '''INTERFERERS'' must be of type iosr.bss.source')
    assert(isnumeric(IVs.azimuths), 'iosr:generateMixtures:invalidAzimuths', '''AZIMUTHS'' must be numeric')
    assert(size(IVs.azimuths,2)==size(interferers,2)+1, 'iosr:generateMixtures:invalidAzimuths', '''AZIMUTHS'' should have one more column than INTERFERERS')
    assert(isnumeric(IVs.elevations), 'iosr:generateMixtures:invalidElevations', '''ELEVATIONS'' must be numeric')
    assert(size(IVs.elevations,2)==size(interferers,2)+1, 'iosr:generateMixtures:invalidElevations', '''ELEVATIONS'' should have one more column than INTERFERERS')
    if ~isempty(IVs.sofa_paths)
        if ischar(IVs.sofa_paths)
            IVs.sofa_paths = cellstr(IVs.sofa_paths);
        elseif ~iscellstr(IVs.sofa_paths)
            error('iosr:generateMixtures:pathInvalid','sofa_paths should be a char array or cell array of strings')
        end
    end
    
    % check settings
    assert(islogical(settings.cache), 'iosr:generateMixtures:invalidCache', '''CACHE'' must be logical')
    assert(ischar(settings.folder), 'iosr:generateMixtures:invalidFolder', '''FOLDER'' must be a char array')
    if ~isempty(settings.combine)
        assert(ischar(settings.combine), 'iosr:generateMixtures:invalidCombine', '''COMBINE'' must be a char array')
    end
    assert(isnumeric(settings.fs) && isscalar(settings.fs), 'iosr:generateMixtures:invalidFs', '''FS'' must be a numeric scalar')
    
    % ensure column vectors
    IVs.sofa_paths = IVs.sofa_paths(:);
    IVs.tirs = IVs.tirs(:);
    targets = targets(:);
    
    % append targets and interferers to IVs
    IVs.targets = targets;
    IVs.interferers = interferers;
    
    %% make mixtures
    
    % check if vars have equal number of rows (ignore scalar of empty)
    [equal,iterations] = check_vars_equal_rows(IVs);
    
    % automatically determine combine method
    if isempty(settings.combine)
        if equal
            settings.combine = 'rows';
        else
            settings.combine = 'all';
        end
    end
    
    var_size = vars_rows(IVs);
    
    % do some things according to combine mode
    switch lower(settings.combine)
        case 'all'
            iterations = prod(max(var_size,1)); % recalculate if rows were equal
            IV_size = [...
                max(numel(IVs.targets),1),...
                max(size(IVs.interferers,1),1),...
                max(numel(IVs.sofa_paths),1),...
                max(numel(IVs.tirs),1),...
                max(size(IVs.azimuths,1),1),...
                max(size(IVs.elevations,1),1)];
        case 'rows'
            assert(equal, 'iosr:generateMixtures:rowsInvalidInput', 'Properties must have an equal number of rows if specifying ''COMBINE'' mode ''ROWS''')
        otherwise
            error('iosr:generateMixtures:unknownCombine','Unknown ''COMBINE'' property')
    end
    
    % create the mixtures
    if settings.cache
        disp('Writing wav files.')
    end
    mixtures(iterations,1) = iosr.bss.mixture; % preallocate
    for m = 1:iterations
        switch lower(settings.combine)
            % get settings for current iteration
            case 'rows'
                % read rows together
                target = copy(getIV('targets',m));
                interferer = copy(getIV('interferers',m));
                sofa_path = getIV('sofa_paths',m);
                tir = getIV('tirs',m);
                azimuths = getIV('azimuths',m);
                elevations = getIV('elevations',m);
            case 'all'
                % work through each variable separately
                [n,p,q,r,s,t] = ind2sub(IV_size,m);
                target = copy(getIV('targets',n));
                interferer = copy(getIV('interferers',p));
                sofa_path = getIV('sofa_paths',q);
                tir = getIV('tirs',r);
                azimuths = getIV('azimuths',s);
                elevations = getIV('elevations',t);
        end
        % apply spatial settings to sources
        target.azimuth = azimuths(1);
        target.elevation = elevations(1);
        for i = 1:length(interferer)
            interferer(i).azimuth = azimuths(i+1);
            interferer(i).elevation = elevations(i+1);
        end
        % calculate mixtures
        mixtures(m,1) = iosr.bss.mixture(...
            target,...
            interferer,...
            'tir',tir,...
            'sofa_path',sofa_path,...
            'fs',settings.fs,...
            'decomposition',settings.decomposition,...
            'gammatone',settings.gammatone,...
            'stft',settings.stft);
        if settings.cache
            mixtures(m,1).write([settings.folder filesep sprintf('mixture-%05d.wav',m)])
        end
    end
    if settings.cache
        disp('Done.')
    end
    
    function val = getIV(field,row)
    %GETIV retrieve independent variable
    
        % function to index into the option field
        index = @(n,F) mod(n-1,size(F,1))+1;
    
        % return data
        if isempty(IVs.(field)) % return empty
            val = [];
        else % return value
            if ~iscellstr(IVs.(field))
                M = index(row,IVs.(field));
                val = IVs.(field)(M,:);
            else
                M = index(row,IVs.(field)(:));
                val = IVs.(field){M};
            end
        end
        
    end
    
end

function opts = overwrite(opts,validOptions,vgin)
%OVERWRITE overwrite the default properties with varargin

    % count arguments
    nArgs = length(vgin);
    if round(nArgs/2)~=nArgs/2
       error('iosr:generateMixtures:nameValuePair','generateMixtures needs propertyName/propertyValue pairs')
    end
    optionNames = fieldnames(opts);
    % overwrite defults
    for pair = reshape(vgin,2,[]) % pair is {propName;propValue}
       IX = strcmpi(pair{1},optionNames); % find match parameter names
       if any(IX)
          % do the overwrite
          opts.(optionNames{IX}) = pair{2};
       else
           if ~any(strcmpi(pair{1},validOptions))
                error('iosr:generateMixtures:unknownOption',['''' pair{1} ''' is not a valid option'])
           end
       end
    end
end

function [equal,num] = check_vars_equal_rows(varargin)
%CHECK_VARS_EQUAL_ROWS check variables have the same number of elements
% 
%   CHECK_VARS_EQUAL_ROWS(A,B,C,...) and CHECK_VARS_EQUAL_ROWS(OPTS) checks
%   the variables A, B, C, ..., or the structure OPTS, to test whether the
%   variables/fields have the same number of elements. Scalar or empty
%   variables/fields are ignored.

    % get number of elements in each variable/field
    nums = max(vars_rows(varargin{:}),1);

    % determine equality
    nums = nums(nums>1);
    
    if ~isempty(nums)
        equal = all(nums==nums(1)); % determine equality
        % return number of elements
        if equal
            num = nums(1);
        else
            num = prod(nums);
        end
    else
        % one row
        equal = true;
        num = 1;
    end

end

function num = vars_rows(varargin)
%VARS_ROWS return the number of rows in each variable
% 
%   VARS_ROWS(A,B,C,...) and VARS_ROWS(OPTS) returns the number of rows in
%   variables A, B, C, ..., or the fields of structure OPTS.

    % convert struct to cell array
    if isstruct(varargin{1}) && length(varargin)==1
        vgin = struct2cell(varargin{1});
    else
        vgin = varargin;
    end

    % return number of elements
    num = cellfun(@(x) size(x,1),vgin);

end
