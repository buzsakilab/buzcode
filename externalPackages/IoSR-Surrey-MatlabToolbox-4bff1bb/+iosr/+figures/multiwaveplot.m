function varargout = multiwaveplot(varargin)
%MULTIWAVEPLOT Stacked line plots from a matrix or vectors
% 
%   Multiwaveplot draws a series of stacked lines (one on top of the other)
%   contained in the rows of an input 2-D matrix. Each line has a
%   designated row on the plot; the first row is plotted at the bottom of
%   the plot.
% 
%   IOSR.FIGURES.MULTIWAVEPLOT(Z) draws a series of lines contained in the
%   rows of the matrix Z.
% 
%   IOSR.FIGURES.MULTIWAVEPLOT(X,Y,Z) plot the lines against the data in X
%   and Y, such that X specifies the common x-data, and Y determines the
%   vertical position of each row. Hence, length(X)==size(Z,2) and
%   length(Y)==size(Z,1).
% 
%   Z may also be a vector, of the same length as X and Y; the data
%   contained in X and Y will be used to construct a matrix for plotting.
% 
%   IOSR.FIGURES.MULTIWAVEPLOT(...,'parameter','value') allows a number of
%   options to be specified. The options are:
% 
%   ({} indicates the default value)
% 
%   'gain'          : {1} | scalar
%       This parameter scales the height of each line. With gain=1
%       (default), the height of the tallest line will be limited so as not
%       to encroach on the adjacent line; all other lines are scaled
%       accordingly. With gain~=1 the height of each line will
%       decrease/increase by a factor gain.
%   'horizonWidth'  : {1} | scalar
%       The plot can be made to narrower (or wider) at the top, to give the
%       impression of a horizon. This may be useful, for example, when y
%       data represent time. With horizonWidth=1 (default), the top and
%       bottom of the plot have the same width. With horizonWidth<1, the
%       plot is narrower at the top by a factor of horizonWidth. With
%       horizonWidth>1, the plot is narrower at the bottom by a factor of
%       1/horizonWidth. Specifying horizonWidth>1 will cause the x-axis to
%       be moved to the top of the plot.
%   'horizonOffset' : {0} | scalar
%       This parameter specifies the vanishing point of the horizon. With
%       horizonOffset=0 (default), the vanishing point is in the centre.
%       With horizonOffset=-1, the vanishing point is on the left of the
%       plot. With horizonOffset=1, the vanishing point is on the right of
%       the plot. Intermediate values specify intermediate vanishing
%       points.
%   'mode' : {'plot'} | 'fill'
%       This parameter plots the data with the specified mode. The default
%       mode depends on gain: for gain<=1, mode defaults to 'plot' and
%       plots lines for each row of Z. For gain>1, mode defaults to 'fill'
%       and instead plots white patch objects for each row of Z, covering
%       the area under each line, such that lines with a lower row index
%       mask those with a higher row index.
%   'reverseY' : {false} | true
%       Determine the order in which data are plotted vertically. With
%       reverseY=false (default), the first row is plotted at the bottom of
%       the plot. With reverseY=true, the first row is plotted at the top
%       of the plot and the axis direction is reversed, but the horizon
%       effect is unchanged. In both cases, lower rows are plotted in front
%       of higher rows.
%
%   H = IOSR.FIGURES.MULTIWAVEPLOT(...) returns a vector of handles to
%   patch (for mode = 'fill') or lineseries (for mode = 'plot') graphics
%   objects, one handle per line.
% 
%   See also FILL, PATCH, PLOT, IMAGESC.

%   Copyright 2016 University of Surrey.

    %% Derive inputs
    
    firstPar = find(cellfun(@(x) ~isnumeric(x),varargin),1,'first');
    if isempty(firstPar)
        firstPar = nargin+1;
    end
    overrides = {};
    
    % get inputs
    switch firstPar-1
        case 1
            Z = varargin{1};
            X = [];
            Y = [];
        case 2
            error('iosr:multiwaveplot:nargin','Wrong number of input arguments.')
        otherwise
            X = varargin{1};
            X = reshape(X,1,length(X));
            Y = varargin{2};
            Y = reshape(Y,1,length(Y));
            Z = varargin{3};
    end
    
    if firstPar < nargin
        overrides = varargin(firstPar:end);
    end
    
    % derive data
    if isvector(Z)
        assert(~isempty(X), 'iosr:multiwaveplot:invalidX', 'If Z is a vector, you must specify X and Y.')
        assert(~isempty(Y), 'iosr:multiwaveplot:invalidY', 'If Z is a vector, you must specify X and Y.')
        % if vector, make matrix
        assert(length(X)==length(Z) && length(Y)==length(Z), 'iosr:multiwaveplot:invalidInputs', 'If Z is a vector, X and Y must be vectors of the same length.')
        x = unique(X)';
        y = unique(Y)';
        wave = NaN(length(y),length(x));
        for n = 1:length(x)
            for m = 1:length(y)
                IX = X==x(n) & Y==y(m);
                wave(m,n) = mean(Z(IX));
            end
        end
        [r,~] = size(wave);
    else
        % if matrix
        if isempty(X)
            X = 1:size(Z,2);
        end
        if isempty(Y)
            Y = 1:size(Z,1);
        end
        assert(ismatrix(Z), 'iosr:multiwaveplot:invalidZ', 'Z must be a 2-D matrix')
        wave = Z;
        [r,c] = size(wave);
        assert(c==length(X), 'iosr:multiwaveplot:invalidX', 'X must be the same length as Z has columns.')
        assert(r==length(Y), 'iosr:multiwaveplot:invalidY', 'Y must be the same length as Z has rows.')
        if isempty(X)
            x = 1:r;
        else
            x = reshape(X,1,length(X));
        end
        if isempty(Y)
            y = 1:c;
        else
            y = Y;
        end
    end
    
    % get options
    
    options = struct(...
        'gain',1,...
        'mode',[],...
        'horizonWidth',1,...
        'horizonOffset',0,...
        'reverseY',false);

    % read parameter/value inputs
    if ~isempty(overrides) % if parameters are specified
        % read the acceptable names
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(overrides);
        if round(nArgs/2)~=nArgs/2
           error('iosr:multiwaveplot:nameValuePair','MULTIWAVEPLOT needs propertyName/propertyValue pairs')
        end
        % overwrite defults
        for pair = reshape(overrides,2,[]) % pair is {propName;propValue}
           IX = strcmpi(pair{1},optionNames); % find match parameter names
           if any(IX)
              % do the overwrite
              options.(optionNames{IX}) = pair{2};
           else
              error('iosr:multiwaveplot:unknownOption','%s is not a recognized parameter name',pair{1})
           end
        end
    end
    
    % set default plot mode
    if isempty(options.mode)
        if options.gain > 1
            options.mode = 'fill';
        else
            options.mode = 'plot';
        end
    end
    
    assert(all(diff(x)>=0), 'iosr:multiwaveplot:invalidX', 'X must be increasing')
    assert(all(diff(y)>=0), 'iosr:multiwaveplot:invalidY', 'Y must be increasing')
    
    %% Plot

    % horizon
    assert(isscalar(options.horizonWidth), 'iosr:multiwaveplot:invalidHorizonWidth', '''horizonWidth'' must be a scalar.')
    assert(isscalar(options.horizonOffset), 'iosr:multiwaveplot:invalidHorizonOffset', '''horizonOffset'' must be a scalar.')
    assert(options.horizonWidth>=0, 'iosr:multiwaveplot:invalidHorizonWidth', '''horizonWidth'' must be greater than or equal to 0.')
    assert(options.horizonOffset>=-1 && options.horizonOffset<=1, 'iosr:multiwaveplot:invalidHorizonOffset', '''horizonOffset'' must be in the range [-1,1].')
    if options.reverseY % correct horizon when axis reversed
        options.horizonWidth = 1/options.horizonWidth;
    end
    % calculate widths
    horizonWidths = linspace(1,options.horizonWidth,length(y));
    horizonWidths = horizonWidths./max(horizonWidths);
    horizonX = linspace(-1,1,length(x))-options.horizonOffset;
    
    % data scaling
    if min(wave(:))>=0
        adjust = 1; % for positive data (e.g. correlograms)
        if options.reverseY
            y0 = y(end);
        else
            y0 = y(1);
        end
    else
        adjust = 0.5; % for wave data varying on zero
        if options.reverseY
            y0 = y(end)+((options.gain*adjust)*(y(end)-y(end-1)));
        else
            y0 = y(1)-((options.gain*adjust)*(y(2)-y(1)));
        end
    end

    wave_max = max(abs(wave(:))); % scale relative to max of wave

    for n = 1:r
        if ~all(wave(n,:)==0) % just in case all values are zero
            wave(n,:) = (adjust/wave_max).*wave(n,:);
        end
    end

    scale = zeros(r,1); % this parameter allows for non-linear y-values, scaling each channel accordingly
    h = zeros(r,1);

    % plot
    ax = newplot;
    hold on;
    
    % plotting order
    if options.reverseY
        set(ax,'ydir','reverse')
        order = 1:r;
        shift = @(z,y) y-z;
    else
        order = r:-1:1;
        shift = @(z,y) z+y;
    end

    for n = order
        % calculate scaling factor
        try
            scale(n) = y(n+1)-y(n);
        catch
            scale(n) = y(n)-y(n-1);
        end
        wave(n,:) = wave(n,:).*scale(n);
        % scale data for horizon effect
        xinterp = interp1(horizonX,x,horizonWidths(n).*horizonX);
        yscale = reshape(shift(options.gain.*wave(n,:),y(n)),1,size(wave,2));
        % plot
        switch options.mode
            case 'plot'
                h(n) = plot(xinterp,yscale,'k');
            case 'fill'
                xhor = [x(1) xinterp(1) xinterp xinterp(end) x(end)];
                yhor = [y(n) y(n) yscale y(n) y(n)];
                xa=[xhor(1) xhor xhor(end) xhor(1)];
                ya=[y0 yhor y0 y0];
                IX = ~(isnan(xa) | isnan(ya));
                h(n) = fill(xa(IX),ya(IX),'w');
            otherwise
                error('iosr:multiwaveplot:unknownMode','Unknown MODE specified: must be ''fill'' or ''plot''')
        end
    end

    % look at every row to search for max, taking y offset into account
    if options.reverseY
        ymin = [y(1)-((options.gain*adjust)*(y(2)-y(1))) y0];
    else
        ymin = [y0 y(end)+((options.gain*adjust)*(y(end)-y(end-1)))];
    end

    % set plot parameters
    hold off
    axis([min(x) max(x) ymin]); % set axis limits
    set(gca,'layer','top') % move axis to top layer
    ytick = get(gca,'ytick');
    set(gca,'ytick',ytick(ytick<=max(y))); % remove ticks beyond y limits
    box on
    if options.horizonWidth>1
        set(gca,'XAxisLocation','top')
    end

    % output
    varargout(1:nargout) = {h};

end
