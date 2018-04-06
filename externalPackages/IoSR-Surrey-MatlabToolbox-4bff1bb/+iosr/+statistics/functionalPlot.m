classdef (Abstract, CaseInsensitiveProperties = true) functionalPlot < ...
        iosr.statistics.statsPlot
%FUNCTIONALPLOT Abstract superclass for functional plots
% 
%   As an abstract class, this class cannot be instantiated. It provides no
%   public methods. The class is a super class for:
%       - iosr.statistics.functionalBoxPlot
%       - iosr.statistics.functionalSpreadPlot
% 
%   The class does, however, provide various properties to subclasses.
%   
%   IOSR.STATISTICS.FUNCTIONALSPREADPLOT properties:
%       mainLineColor           - Specifies the color of the central line.
%                                 See 'addPrctilesLineColor' for details of
%                                 how to specify colors. The default is
%                                 'auto'.
%       mainLineStyle           - Specifies the style of the central line.
%                                 The property may be specified as a char
%                                 array (e.g. '-'), or as a length-P cell
%                                 vector of strings (if passing line styles
%                                 for each functional plot). The default is
%                                 '-'.
%       mainLineWidth           - Specifies the width of the central line.
%                                 The property may be specified as a
%                                 length-P numeric vector, or as a cell
%                                 vector (if passing line styles for each
%                                 functional plot). The default is 1.
%       outerLineColor          - Specifies the color of the outer lines.
%                                 See 'mainLineStyle' for details of how to
%                                 specify colors. The default is 'auto'.
%       outerLineStyle          - Specifies the style of the outer lines.
%                                 See 'mainLineStyle' for details of how to
%                                 specify styles. The default is '-'.
%       outerLineWidth          - Specifies the width of the outer line.
%                                 See 'mainLineWidth' for details of how to
%                                 specify widths. The default is 1.
%       outlierLineWidth        - Specifies the width of the outlier marker
%                                 line. See 'mainLineWidth' for details of
%                                 how to specify widths. The default is 1.
%       showOutliers            - Specifies whether to show outliers. The
%                                 property must be a logical value. The
%                                 default is false.
%       spreadAlpha             - Set the alpha of the shaded region(s).
%                                 The default is 0.5.
%       spreadBorderLineColor   - Set the color of the shaded regions'
%                                 / region's border. See
%                                 'addPrctilesLineColor' for details of how
%                                 to specify colors. The default is 'none'.
%       spreadBorderLineWidth   - Set the width of the shaded regions'
%                                 / region's border. The default is 1.
%       spreadColor             - Set the color of the shaded region(s).
%                                 See 'addPrctilesLineColor' for details of
%                                 how to specify colors. The default is
%                                 'auto'.
%       whiskers                - Specify whether the box is connected to
%                                 the outer line via a whisker. The line is
%                                 drawn with the same style as the outer
%                                 line. The default is false.
% 
%   These properties can be referenced using dot notation - e.g. H.BOXCOLOR
%   where H is an instance of the FUNCTIONALSPREADPLOT object - or using
%   the SET and GET methods - e.g. GET(H,'BOXCOLOR'). Both methods are
%   case-insensitive.
% 
%   Note that some handles will be empty unless the associated option is
%   turned on.
% 
%   See also IOSR.STATISTICS.BOXPLOT, IOSR.STATISTICS.QUANTILE,
%   IOSR.STATISTICS.FUNCTIONALBOXPLOT,
%   IOSR.STATISTICS.FUNCTIONALSPREADPLOT, IOSR.STATISTICS.STATSPLOT.

    properties (AbortSet)
        mainLineColor = 'auto'          % The color of the central line.
        mainLineStyle = '-'             % The style of the central line.
        mainLineWidth = 2               % The width of the central line.
        outerLineColor = 'auto'         % The color of the outer lines.
        outerLineStyle = '-'            % The style of the outer lines.
        outerLineWidth = 0.5            % The width of the outer lines.
        outlierLineWidth = 1            % Width of the outlier marker edges.
        showOutliers = true             % Turn outliers on and off.
        spreadAlpha = 0.5               % The alpha of the shaded regions.
        spreadColor = 'auto'            % The color of the shaded regions.
        spreadBorderLineWidth = 1       % The width of the shaded region borders.
        spreadBorderLineColor = 'none'  % The color of the shaded region borders.
        whiskers = false                % Specify whether a line connects the box to the outer lines.
    end
    
    methods
        
        function set.showOutliers(obj, val)
            assert(numel(val)==1 && islogical(val), 'iosr:functionalPlot:invalidShowOutliers', ...
                '''SHOWOUTLIERS'' must be true or false.')
            obj.showOutliers = val;
            obj.draw();
        end
        
        % main lines
        
        function val = get.mainLineColor(obj)
            val = obj.parseColor(obj.mainLineColor);
        end
        
        function set.mainLineColor(obj,val)
            obj.mainLineColor = val;
            obj.draw();
        end
        
        function val = get.mainLineStyle(obj)
            val = obj.parseProps(obj.mainLineStyle, false);
        end
        
        function set.mainLineStyle(obj,val)
            obj.mainLineStyle = val;
            obj.draw();
        end
        
        function val = get.mainLineWidth(obj)
            val = obj.parseProps(obj.mainLineWidth, false);
        end
        
        function set.mainLineWidth(obj,val)
            obj.mainLineWidth = val;
            obj.draw();
        end
        
        % outer lines
        
        function val = get.outerLineColor(obj)
            val = obj.parseColor(obj.outerLineColor);
        end
        
        function set.outerLineColor(obj,val)
            obj.outerLineColor = val;
            obj.draw();
        end
        
        function val = get.outerLineStyle(obj)
            val = obj.parseProps(obj.outerLineStyle, false);
        end
        
        function set.outerLineStyle(obj,val)
            obj.outerLineStyle = val;
            obj.draw();
        end
        
        function val = get.outerLineWidth(obj)
            val = obj.parseProps(obj.outerLineWidth, false);
        end
        
        function set.outerLineWidth(obj,val)
            obj.outerLineWidth = val;
            obj.draw();
        end
        
        % outlier settings

        function val = get.outlierLineWidth(obj)
            val = obj.parseProps(obj.outlierLineWidth, false);
        end
        
        function set.outlierLineWidth(obj,val)
            obj.outlierLineWidth = val;
            obj.draw();
        end
        
        % set spread alpha
        
        function set.spreadAlpha(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:functionalPlot:invalidSpreadAlpha', ...
                '''SPREADALPHA'' must be a numeric scalar')
            obj.spreadAlpha = val;
            obj.draw();
        end
        
        % get/set spread line colors
        
        function val = get.spreadBorderLineColor(obj)
            val = obj.parseColor(obj.spreadBorderLineColor);
        end
        
        function set.spreadBorderLineColor(obj,val)
            obj.spreadBorderLineColor = val;
            obj.draw();
        end
        
        % set spread border line width
        
        function set.spreadBorderLineWidth(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:functionalPlot:invalidspreadBorderLineWidth', ...
                '''SPREADBORDERLINEWIDTH'' must be a numeric scalar')
            obj.spreadBorderLineWidth = val;
            obj.draw();
        end
        
        % set/get spread color
        
        function val = get.spreadColor(obj)
            val = obj.parseColor(obj.spreadColor);
        end
        
        function set.spreadColor(obj,val)
            obj.spreadColor = val;
            obj.draw();
        end
        
        % whisker
        
        function set.whiskers(obj,val)
            assert(islogical(val) && isscalar(val), 'iosr:functionalPlot:invalidWhiskers', ...
                '''whisker'' must be a scalar and logical');
            obj.whiskers = val;
            obj.draw();
        end
        
    end
    
    methods (Access = protected)
        
        % Draw the plot
        function draw(obj)
            
            if isfield(obj.handles,'axes')
            
                axes(obj.handles.axes);
                hold on;

                if length(obj.outDims) > 2
                    nlines = obj.outDims(3);
                else
                    nlines = 1;
                end
                for n = nlines:-1:1
                    % draw inner patch
                    obj.handles.spreads(n) = patch(...
                        [obj.x obj.x(end:-1:1)], ...
                        [obj.statistics.inner_u(:,:,n) obj.statistics.inner_l(:,end:-1:1,n)], ...
                        obj.spreadColor{n} , ...
                        'FaceAlpha', obj.spreadAlpha, ...
                        'EdgeColor', obj.spreadBorderLineColor{n}, ...
                        'LineWidth', obj.spreadBorderLineWidth ...
                    );
                end

                for n = nlines:-1:1

                    % draw outer limits
                    obj.handles.outer_u(n) = line(obj.x, obj.statistics.outer_u(:,:,n), ...
                        'linestyle', obj.outerLineStyle{n}, ...
                        'linewidth', obj.outerLineWidth{n}, ...
                        'color', obj.outerLineColor{n} ...
                    );
                    obj.handles.outer_l(n) = line(obj.x, obj.statistics.outer_l(:,:,n), ...
                        'linestyle', obj.outerLineStyle{n}, ...
                        'linewidth', obj.outerLineWidth{n}, ...
                        'color', obj.outerLineColor{n} ...
                    );

                    % draw centre line
                    obj.handles.main(n) = line(obj.x, obj.statistics.main(:,:,n), ...
                        'linestyle', obj.mainLineStyle{n}, ...
                        'linewidth', obj.mainLineWidth{n}, ...
                        'color', obj.mainLineColor{n} ...
                    );

                end

                % draw whiskers
                if obj.whiskers
                    x_range = max(obj.x) - min(obj.x);
                    whisker_xlim = [(min(obj.x) + (1/3)*x_range) (max(obj.x) - (1/3)*x_range)];
                    if nlines==1
                        whisker_x = mean(whisker_xlim);
                    else
                        whisker_x = linspace(whisker_xlim(1),whisker_xlim(2),nlines);
                    end
                    for n = nlines:-1:1
                        % interpolate to find y
                        [~,x_ix] = min(abs(obj.x-whisker_x(n)));
                        x_interp_ix = [max(1,x_ix-1) x_ix min(length(obj.x),x_ix+1)];
                        x_interp_ix = unique(x_interp_ix);
                        y_i_l = interp1(obj.x(x_interp_ix), obj.statistics.inner_l(:,x_interp_ix,n), whisker_x(n));
                        y_i_u = interp1(obj.x(x_interp_ix), obj.statistics.inner_u(:,x_interp_ix,n), whisker_x(n));
                        y_o_l = interp1(obj.x(x_interp_ix), obj.statistics.outer_l(:,x_interp_ix,n), whisker_x(n));
                        y_o_u = interp1(obj.x(x_interp_ix), obj.statistics.outer_u(:,x_interp_ix,n), whisker_x(n));

                        % plot
                        obj.handles.whiskers_l(n) = line([whisker_x(n) whisker_x(n)],[y_i_l y_o_l], ...
                            'linestyle', obj.outerLineStyle{n}, ...
                            'linewidth', obj.outerLineWidth{n}, ...
                            'color', obj.outerLineColor{n});
                        obj.handles.whiskers_u(n) = line([whisker_x(n) whisker_x(n)],[y_i_u y_o_u], ...
                            'linestyle', obj.outerLineStyle{n}, ...
                            'linewidth', obj.outerLineWidth{n}, ...
                            'color', obj.outerLineColor{n});
                    end
                else
                    try
                        delete(obj.handles.whiskers_l(n))
                        delete(obj.handles.whiskers_u(n))
                    catch
                    end
                end
            end
                
        end
        
    end
    
    methods (Access = protected)
        
        % Make properties in to a cell array
        function val = parseProps(obj, prop, addP)
            
            nLines = obj.outDims(3);
            if ~addP
                nAdd = 1;
            else
                nAdd = length(obj.addPrctiles);
            end
            val = cell(nLines, nAdd);
            if isnumeric(prop)
                cellprop = num2cell(prop);
            elseif ischar(prop)
                cellprop = cellstr(prop);
            else
                cellprop = prop;
            end
            try
            for n = 1:nLines
                if ~addP
                    val(n) = cellprop(mod(n-1, numel(cellprop)) + 1);
                else
                    for p = 1:nAdd
                        ixn = mod(n-1, size(cellprop, 1)) + 1;
                        ixp = mod(p-1, size(cellprop, 2)) + 1;
                        val(n,p) = cellprop(ixn, ixp);
                    end
                end
            end
            catch
                keyboard
            end
            
        end
        
        % Put colors in to a cell array
        function val = parseColor(obj, inColor)
            
            nLines = obj.outDims(3);
            val = cell(nLines, 1);
            if isnumeric(inColor)
                if size(inColor, 2) ~= 3
                    error('iosr:functionalPlot:colorInvalid','Color must be an N-by-3 array, where N is any positive integer.')
                end
                for n = 1:nLines
                    val(n,:) = {inColor(mod(n, size(inColor, 1)) + 1, :)};
                end
            elseif ischar(inColor)
                if strcmp(inColor, 'auto')
                    colors = lines(nLines);
                    for n = 1:nLines
                        val(n,:) = {colors(n, :)};
                    end
                else
                    colors = cellstr(inColor);
                    for n = 1:nLines
                        val(n,:) = colors(mod(n, size(inColor, 1)) + 1);
                    end
                end
            elseif isa(inColor,'function_handle')
                colors = inColor(nLines);
                for n = 1:nLines
                    val(n,:) = {colors(mod(n, numel(inColor)) + 1)};
                end
            elseif iscellstr(inColor)
                for n = 1:nLines
                    val(n,:) = inColor(mod(n, numel(inColor)) + 1);
                end
            else
                error('iosr:functionalPlot:invalidColorFormat','Invalid color format specified.')
            end
            
        end
        
    end
    
end
