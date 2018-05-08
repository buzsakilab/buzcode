function [y,x,g] = tab2box(Xin,Yin,Gin)
%TAB2BOX Prepare tabular data for boxPlot function
% 
%   Y = IOSR.STATISTICS.TAB2BOX(XIN,YIN) prepares data, in tabular form,
%   for use in the BOX_PLOT function. Specifically, XIN and YIN are vectors
%   (for example column vectors from a results table). Y is an N-by-P
%   numeric array, where P is the number of unique elements in XIN, and N
%   is the maximum number of occurences of any individual element of XIN.
%   In cases where elements in XIN do not occur an equal number of times,
%   columns in YIN are padded with NaNs.
% 
%   Y = IOSR.STATISTICS.TAB2BOX(XIN,YIN,GIN) returns an
%   N-by-P-by-G-by-I-by-J... numeric array, where G is the number of unique
%   elements in the first column of GIN, I is the number of unique elements
%   in the second column of GIN, etc. GIN is a numeric or cell matrix with
%   as many rows as XIN or YIN. The input facilitates hierarchical grouping
%   of the data in the box plot, with the grouping order determined by the
%   column order of GIN.
% 
%   [Y,X] = IOSR.STATISTICS.TAB2BOX(...) returns the unique values in XIN
%   to X.
% 
%   [Y,X,G] = IOSR.STATISTICS.TAB2BOX(...) returns the unique values in GIN
%   to G. G is a cell vector whereby the Nth element contains a vector of
%   length SIZE(Y,N+2).
% 
%   Example
% 
%       % Prepare data for box_plot grouped by two variables
%   
%       % load data
%       % (requires Statistics or Machine Learning Toolbox)
%       load carbig
% 
%       % arrange data
%       [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);
%   
%       % sort
%       IX = [1 3 2]; % order
%       g = g{1}(IX);
%       y = y(:,:,IX);
% 
%       % plot
%       figure
%       h = iosr.statistics.boxPlot(x,y,...
%           'boxColor','auto','medianColor','k',...
%           'scalewidth',true,'xseparator',true,...
%           'groupLabels',g,'showLegend',true);
%       box on
%       title('MPG by number of cylinders and period')
%       xlabel('Number of cylinders')
%       ylabel('MPG')
% 
%   See also IOSR.STATISTICS.BOXPLOT.

%   Copyright 2016 University of Surrey.

    %% validate input

    % validate Xin
    if ischar(Xin)
        Xin = cellstr(Xin);
    end
    assert(isvector(Xin), 'iosr:tab2box:invalidInput', 'Xin must be a vector')

    % validate Yin
    assert(isvector(Yin), 'iosr:tab2box:invalidInput', 'Yin must be a vector')
    assert(isnumeric(Yin), 'iosr:tab2box:invalidInput', 'Yin must be numeric')

    % validate Gin
    if nargin<3
        Gin = [];
    else
        if ischar(Gin)
            Gin = cellstr(Gin);
        end
        assert(isvector(Gin) || size(Gin,1)==numel(Yin), 'iosr:tab2box:invalidInput', 'Gin must be a vector or a matrix with as many rows as Y has elements.')
        if isvector(Gin)
            Gin = Gin(:);
        end
    end

    %% group

    % unique values
    x = getUniques(Xin);
    if iscell(x)
        x = x{1};
    end
    g = getUniques(Gin);
    if ~iscell(g)
        temp = cell(1,size(g,2));
        for c = 1:size(g,2);
            temp{c} = g(:,c);
        end
        g = temp;
    end

    % preallocate cell array
    gdims = cellfun(@length,g);
    dims = [1,length(x),gdims];
    yc = cell(dims);
    subgidx = cell(1,length(gdims));

    % put data into cells
    for a = 1:length(x)
        IXx = findIX(Xin,x,a);
        for b = 1:prod(gdims)
            [subgidx{:}] = ind2sub(gdims, b); % get group index
            subidx = [{1} {a} subgidx{:}]; % make main index
            IXg = true(size(Xin)); % select everything at first
            for c = 1:length(subgidx) % narrow it down
                IXg = IXg & findIX(Gin(:,c),g{c},subgidx{c});
            end
            % return y as column
            yc{subidx{:}} = reshape(Yin(IXx & IXg),sum(IXx & IXg),1);
        end
    end

    %% convert to numeric array

    try % see if can concat directly
        y = cell2mat(yc);
    catch % else pad with NaNs
        maxSizeY = max(cellfun(@length,yc(:)));
        for a = 1:numel(yc)
            if length(yc{a}) < maxSizeY
                yc{a} = [yc{a}; NaN(maxSizeY-length(yc{a}),1)];
            end
        end
        y = cell2mat(yc);
    end

end

function IX = findIX(L,l,n)
%FINDIX find entry in array from lookup array and index

    if iscellstr(L)
        IX = strcmp(L,l{n});
    elseif iscell(L)
        IX = cell2mat(L)==l(n);
    else
        IX = L==l(n);
    end
end

function out = getUniques(d)
%GETUNIQUES return unique entries from vector

    if isvector(d) % ensure column vector
        d = d(:);
    end
    dims = zeros(1,size(d,2));
    u = cell(1,size(d,2));
    % put unique values into columns
    for c = 1:size(d,2)
        if iscellstr(d(:,c)) || isnumeric(d(:,c))
            u{c} = unique(d(:,c));
        else
            u{c} = unique(cell2mat(d(:,c)));
        end
        dims(c) = length(u{c});
    end
    try % to make a numeric array
        out = cell2mat(u);
    catch % keep as cell array
        out = u;
    end

end
