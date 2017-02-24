function makeFigure(fh)
%
% fh must return structure array (one element per figure  with the following fields
%
% n(needed): the data to plot, intended as the y-axis for 'plot' and 'errorbar'
% plots, the height of the bars for 'hist' plots, and the data to be
% colorcoded for 'image' plots
% x(optional): data for the x-axis
% y(optional): data for the y-axis for image plots
% style(optional): style for the line object
% lineProperties(optional): properties that will be applied to the line
% objects
%
% figureName (optional): the name to use for all output
% figureType (needed): the type of figure accepted
% xLabel (optional): label for the x axis
% ylabel (optional): label for the y axis
% nFigure(optional): figure number
% xLim(optional): limits for the x-axis
% yLim(optional): limits for the y-axis
% axesProperties(optional): properties that will be set on the
% xTick(optional): ticks for the x axis
% xTickLabel(optional): tick label for the x axis
% noXTick(optional): if set to non-zero, gets rid of xticks
% yTick(optional): ticks for the y axis
% yTickLabel(optional): tick label for the y axis
% noYTick(optional): if set to non-zero, gets rid of yticks
% title(optional): title for the figure
% legend(optional): legend for the figure
% position(optional): position for the figure


if nargin < 2
    renew = 0;
end





global FIGURE_DIR;


if isempty(FIGURE_DIR)
    FIGURE_DIR = pwd;
end

figure_dir = FIGURE_DIR;



global N_FIGURE




if iscell(fh)
    fig_st_array = fh;
else
    fig_st_array = feval(fh);
end
for nf = 1:length(fig_st_array)

    if isempty(N_FIGURE)
        N_FIGURE = 1;
    end

    n_figure = N_FIGURE;
    fig_st = fig_st_array{nf};

    figure_name = checkField('figureName', 0, ['mkf_figure' num2str(n_figure)]);


    mat_name = [figure_dir filesep figure_name '.mat'];
    fig_name = [figure_name  '.fig'];
    png_name = [figure_name  '.png'];




    font_name = 'Helvetica';
    font_size = 14;
    font_weight = 'bold';
    line_width = 2;


    fig_type = checkField('figureType');


    xlab = checkField('xLabel', 0);
    ylab = checkField('yLabel', 0);

    figTitle = checkField('title', 0);

    figLegend = checkField('legend', 0);
    position = checkField('position', 0);
    
    xtick = checkField('xTick', 0);
    xticklabel = checkField('xTickLabel', 0);
    no_xtick = checkField('noXTick', 0);
    ytick = checkField('yTick', 0);
    yticklabel = checkField('yTickLabel', 0);
    no_ytick = checkField('noYTick', 0);

    xlim = checkField('xLim', 0);
    ylim = checkField('yLim', 0);
    axes_properties = checkField('axesProperties', 0);


    n = checkField('n');


    nfigg = checkField('nFigure', 0);
    if ~isempty(nfigg)
        n_figure = nfigg;
    end



    switch fig_type

        case 'image'
            x = checkField('x', 0 , 1:(size(n,2)));

        otherwise
            x = checkField('x', 0 , 1:(size(n,1)));
    end



    y = checkField('y', 0, 1:(size(n,1)));

    if strcmp(fig_type, 'errorbar') | strcmp(fig_type, 'histerror')
        e = checkField('e');
        ehi = checkField('eHi', 0);
    end

    if strcmp(fig_type, 'histerror')
        eTop = checkField('eBarsOnTop', 0);
    end
    
    lp = (size(n,1));

    for i = 1:lp
        default_style{i} = 'k-';
    end


    style = checkField('style', 0, default_style);

    for i = 1:length(style)
        if isempty(regexp(style{i}, '[bgrcmyk]'))
            style{i} = ['k' style{i}];
        end
    end


    line_properties = checkField('lineProperties', 0);







    






    % put drawing code here

    lh = {};
    figure(n_figure);
    clf



    switch(fig_type)
        case 'hist'

            nbars = size(n,2);
            gr = linspace(0.2, 1, nbars);


            h = bar(x, n);

            for i = 1:nbars

                set(h(i), 'FaceColor', [gr(i) gr(i) gr(i)]);
                set(h(i), 'LineWidth', line_width);
            end
        case 'image'

            [n1, n2] = size(n);

            if ~exist('x')
                x = 1:n1;
            end

            if ~exist('y')
                y = 1:n2;
            end

            imagesc(x, y, n);
           % caxis(cax);

            axis equal
            axis([(min(x)-0.5) (max(x)+0.5) (min(y)-0.5) (max(y)+0.5)]);




        case 'plot'
            nplot = length(x);

            for i = 1:nplot
                lh {i} = plot(x{i}, n{i}, style{i});
                hold on
            end

            hh = findobj(gca, 'type', 'line');
            set(hh, 'LineWidth', line_width);

        case 'errorbar'
            nplot = length(x);
            for i = 1:nplot
                if ~isempty(e{i})
                    if ~isempty(ehi)
                        lh{i} = errorbar(x{i}, n{i}, e{i}, ehi{i}, style{i});
                    else
                        lh{i} = errorbar(x{i}, n{i}, e{i}, style{i});
                    end
                else
                    plot(x{i}, n{i}, style{i});
                end

                hold on
            end

            hh = findobj(gca, 'type', 'line');
            set(hh, 'LineWidth', line_width);


        case 'histerror'

            nbars = size(n,2);
            if nbars > 1
                gr = linspace(0.2, 1, nbars);
            else
                gr = 0.5;
            end


            if eTop
                ehi = e;
                e = zeros(size(e));
            end
            
            h = bar(x, n);

            for i = 1:nbars

                set(h(i), 'FaceColor', [gr(i) gr(i) gr(i)]);
                set(h(i), 'LineWidth', line_width);
                hh = get(h(i), 'children');
                xd = get(hh, 'XData');
                xd = (xd(1,:)+xd(end,:))/2;
                
                
                xe{i} = xd';
                
                ne{i} = n(:,i);
                err{i} = e(:,i);
                
                if ~isempty(ehi)
                    ehirr{i} = ehi(:,i); 
                end
                
            end

            e = err;
            if ~isempty(ehi)
                ehi = ehirr;
            end
            
            
            hold on
            
            for i = 1:nbars
                if ~isempty(ehi)
                    lh{i} = errorbar(xe{i}, ne{i}, e{i}, ehi{i}, style{i});
                else
                    lh{i} = errorbar(xe{i}, ne{i}, e{i}, style{i});
                end
                set(lh{i}, 'LineStyle',  'none');
                set(lh{i}, 'LineWidth' ,2);
            end

    end



    if ~isempty(xlim)
        set(gca, 'XLim', xlim);
    end

    if ~isempty(ylim)
        set(gca, 'YLim', ylim);
    end




    set(gca, 'LineWidth', line_width);

    if ~isempty(line_properties) & ~isempty(lh)
        for i = 1:length(lh)
            if ~isempty(line_properties{i})
                lp = line_properties{i};
                set(lh{i}, lp{:});
            end
        end
    end



    if ~isempty(xlab)
        xlabel(xlab, 'FontName', font_name, 'FontWeight', font_weight, ...
            'FontSize', font_size);
    end

    if ~isempty(ylab)
        ylabel(ylab, 'FontName', font_name, 'FontWeight', font_weight, ...
            'FontSize', font_size);
    end

    if ~isempty(figTitle)
        title(figTitle, 'FontName', font_name, 'FontWeight', font_weight, ...
            'FontSize', font_size);
    end

    if ~isempty(figLegend)
        legend(figLegend);
        % allow an extra 20% of vertical space for the legend, since
        % MATLAB doesn't do it itself... 
        yl = get(gca, 'YLim');
        yl2 = yl(2) + (yl(2)-yl(1))*0.2;
        yl = [yl(1) yl2];
        set(gca, 'YLim', yl);
    end


    if ~isempty(xticklabel)
        set(gca, 'XTickLabel', xticklabel);
    end

    if ~isempty(yticklabel)
        set(gca, 'YTickLabel', yticklabel);
    end
    
    if ~isempty(xtick)
        set(gca, 'XTick', xtick);
    end

    if ~isempty(ytick)
        set(gca, 'YTick', ytick);
    end

    if no_xtick
        set(gca, 'XTickLabel', {});
    end

    if no_ytick
        set(gca, 'YTickLabel', {});
    end

    if ~isempty(position)
        set(gcf, 'position', position);
    end
    



    if ~isempty(axes_properties)
        set(gca, axes_properties{:});
    end


    set(gca, 'FontName', font_name);
    set(gca, 'FontWeight', font_weight);
    set(gca, 'FontSize', font_size);


    % saving everything

dbstop error
    saveas(gcf, [figure_dir filesep fig_name], 'fig');
    print('-dpng', [figure_dir filesep png_name]);

    N_FIGURE = N_FIGURE + 1;
          
end








function value = checkField(field_name, needed, default)

if nargin < 2
    needed = 1;
end

if nargin < 3
    default = [];
end



fig_st = evalin('caller', 'fig_st;');

if isfield(fig_st, field_name)
    value = eval(['fig_st.' field_name ';']);
elseif needed
    error(['missing needed field ' field_name]);
else
    value = default;
end