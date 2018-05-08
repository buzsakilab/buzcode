function h = qqPlot(varargin)
%QQPLOT Quantile-quantile plot with patch option
%
%   IOSR.STATISTICS.QQPLOT(Y) displays a quantile-quantile plot of the
%   sample quantiles of Y versus theoretical quantiles from a normal
%   distribution. If the distribution of Y is normal, the plot will be
%   close to linear.
% 
%   IOSR.STATISTICS.QQPLOT(X,Y) displays a quantile-quantile plot of two
%   samples. If the samples come from the same distribution, the plot will
%   be linear.
% 
%   The inputs X and Y should be numeric and have an equal number of
%   elements; every element is treated as a member of the sample.
% 
%   The plot displays the sample data with the plot symbol 'x'.
%   Superimposed on the plot is a dashed straight line connecting the first
%   and third quartiles.
% 
%   IOSR.STATISTICS.QQPLOT(...,MODE) allows the appearance of the plot to
%   be configured. With MODE='line' (default), the plot appears as
%   described above. With MODE='patch', the data are plotted as a patch
%   object, with the area bound by the x-distribution and the linear fit
%   shaded grey. With mode='both' the two appearances are combined.
%   
%   IOSR.STATISTICS.QQPLOT(...,MODE,METHOD) and
%   IOSR.STATISTICS.QQPLOT(...,[],METHOD) allows the method for calculating
%   the quartiles, used for the fit line, to be specified. The default is
%   'R-8'. Type 'help iosr.statistics.quantile' for more information. The
%   latter form of the function call uses the default mode.
% 
%   H = IOSR.STATISTICS.QQPLOT(...) returns a two- or three-element vector
%   of handles to the plotted object. The nature of the handles depends
%   upon the mode. In all cases, the first handle is to the sample data,
%   the second handle is to the fit line. With MODE='patch' or MODE='both',
%   there is third handle to the patch object.
% 
%   Example
% 
%      % Display Q-Q plots for the rand and randn functions
%      figure
%      subplot(2,1,1)
%      iosr.statistics.qqPlot(rand(20),'patch')
%      subplot(2,1,2)
%      h = iosr.statistics.qqPlot(randn(20),'patch');
%      set(h(3),'FaceColor','r') % change fill color
% 
%   See also IOSR.STATISTICS.QUANTILE, IOSR.STATISTICS.BOXPLOT.

%   Copyright 2016 University of Surrey.

    %% determine X and Y

    IXn = cellfun(@(x) isnumeric(x) & ~isempty(x),varargin);

    switch sum(IXn)
        case 0
            error('iosr:qqPlot:noData','No input data specified')
        case 1
            % compare to normal distrbution
            Y = get_input_sample(varargin,IXn);
            p = (.5:length(Y))/length(Y);
            X = sqrt(2)*erfinv(2*p - 1);
            x_label = 'Standard normal quantiles';
            y_label = 'Sample quantiles';
        case 2
            % compare to input data distribution
            Y = get_input_sample(varargin,find(IXn,1,'last'));
            X = get_input_sample(varargin,find(IXn,1,'first'));
            assert(isequal(size(X),size(Y)), 'iosr:quantile:invalidInput', 'Input data must be the same size')
            x_label = 'X quantiles';
            y_label = 'Y quantiles';
        otherwise
            error('iosr:qqPlot:unkonwnInput','Unknown input specified')
    end

    %% determine mode and method

    % find inputs
    IXc = cellfun(@(x) ischar(x) | isempty(x),varargin);
    switch sum(IXc)
        case 0
            mode = [];
            method = [];
        case 1
            mode = varargin{IXc};
            method = [];
        case 2
            mode = varargin{find(IXc,1,'first')};
            method = varargin{find(IXc,1,'last')};
        otherwise
            error('iosr:qqPlot:unknownString','Unknown string specified')
    end

    % defaults
    if isempty(mode)
        mode = 'line';
    end
    if isempty(method)
        method = 'R-8';
    end

    %% calculate fit to first and third quartile

    % quartiles
    q1x = iosr.statistics.quantile(X,.25,[],method);
    q3x = iosr.statistics.quantile(X,.75,[],method);
    q1y = iosr.statistics.quantile(Y,.25,[],method);
    q3y = iosr.statistics.quantile(Y,.75,[],method);

    % slope
    slope = (q3y-q1y)./(q3x-q1x);
    centerx = (q1x+q3x)/2;
    centery = (q1y+q3y)/2;

    % fit
    maxx = max(X);
    minx = min(X);
    maxy = centery + slope.*(maxx - centerx);
    miny = centery - slope.*(centerx - minx);

    % lines
    X_fit = linspace(minx,maxx,length(X));
    Y_fit = linspace(miny,maxy,length(X));

    %% plot data

    hold on

    if strcmpi(mode,'patch') || strcmpi(mode,'both')
        Hp = patch([X fliplr(X_fit)],[Y fliplr(Y_fit)],[0.5 0.5 0.5]);
        set(Hp,'edgecolor','none')
    else
        Hp = NaN;
    end

    switch lower(mode)
        case 'line'
            linestyle = 'xk';
        case 'patch'
            linestyle = '-k';
        case 'both'
            linestyle = 'xk';
        otherwise
            error('iosr:qqPlot:unknownMode','Unknown mode specified')
    end

    H = plot(X,Y,linestyle,X_fit,Y_fit,'--k');

    hold off

    % axis labels
    xlabel(x_label)
    ylabel(y_label)

    box on; axis tight
    set(gca,'layer','top')

    % return handle
    if nargout>0
        h = H;
        if isobject(Hp) || ishandle(Hp)
            h = [h; Hp];
        end
    end

end

function Z = get_input_sample(z,IX)
%GET_INPUT_SAMPLE get sample, order, and convert to vector
    Z = z{IX};
    Z = sort(Z(:))';
end
