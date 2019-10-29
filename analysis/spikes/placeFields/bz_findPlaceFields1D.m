function [mapStats] = bz_findPlaceFields1D(firingMaps,varargin)
%   [stats] = bz_findPlaceFields1D(firingMaps)
%   Find place fields from 1D firing maps. Reads the output of bz_firingMapAvg 

%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this size are considered spurious
%                   and ignored (default = 10)
%     'maxSize'     fields larger than this size are considered spurious
%                   and ignored (default = 50)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1)
%     'verbose'     display processing information (default = 'off')
%    =========================================================================

%% Parse inputs 
p=inputParser;
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',10,@isnumeric);
addParameter(p,'minPeak',1,@isnumeric);
addParameter(p,'verbose','off',@isstr);

parse(p,varargin{:});
threshold = p.Results.threshold;
minSize = p.Results.minSize;
minPeak = p.Results.minPeak;
verbose = p.Results.verbose;


%% Find place fields
for unit = 1:length(firingMaps.rateMaps)
    for c = 1:length(firingMaps.rateMaps{1})
        
        % Default values
            mapStats{unit,1}{c}.x = NaN;
            mapStats{unit,1}{c}.field = [];
            mapStats{unit,1}{c}.size = 0;
            mapStats{unit,1}{c}.peak = 0;
            mapStats{unit,1}{c}.mean = 0;
            mapStats{unit,1}{c}.fieldX = [NaN NaN];
            mapStats{unit,1}{c}.specificity = 0;
            mapStats{unit,1}{c}.m = nan;
            mapStats{unit,1}{c}.r = nan;
            mapStats{unit,1}{c}.mode = nan;
            mapStats{unit,1}{c}.k = nan;

    % Determine the field as the connex area around the peak where the value or rate is > threshold*peak
    % There can be two or more fields
    z = firingMaps.rateMaps{unit}{c};
    x = 1:1:length( firingMaps.rateMaps{1}{1});

    if max(max(z)) == 0,
      mapStats{unit,1}{c}.field = logical(zeros(size(z)));
      return;
    end

    nBinsX = max([1 length(x)]);	% minimum number of bins is 1
    circX = 0; circY = 0;
% Each time we find a field, we will remove it from the map; make a copy first
% Try to find more fields until no remaining bin exceeds min value
    i=1;
    while true,
        % Are there any candidate (unvisited) peaks left?
        [peak,idx] = max(z(:));
        if peak < minPeak, break; end
        % Determine coordinates of largest candidate peak
        [y,x] = ind2sub(size(z),idx);
        % Find field (using min threshold for inclusion)
        field1 = FindFieldHelper(z,x,y,peak*threshold,circX,circY);
        size1 = sum(field1(:));
        % Does this field include two coalescent subfields?
        % To answer this question, we simply re-run the same field-searching procedure on the field
        % we then either keep the original field or choose the subfield if the latter is less than
        % 1/2 the size of the former
        m = peak*threshold;
        field2 = FindFieldHelper(z-m,x,y,(peak-m)*threshold,circX,circY);
        size2 = sum(field2(:));
        if size2< 1/2*size1,
            field = field2;
            tc = ' ';sc = '*'; % for debugging messages
        else
            field = field1;
            tc = '*';sc = ' '; % for debugging messages
        end
        % Display debugging info
%         if verbose,
%             disp([int2zstr(i,2) ') peak  ' num2str(peak) ' @ (' int2str(x) ',' int2str(y) ')']);
%             disp([' ' tc ' field size       ' int2str(size1)]);
%             disp([' ' sc ' subfield size    ' int2str(size2)]);
%             disp(' ');
%             figure;
%             if nDims == 1,
%                 plot(z);hold on;
%                 PlotIntervals(ToIntervals(field1),'rectangles');
%                 PlotIntervals(ToIntervals(field2),'bars');
%                 ylabel(tc);
%             else
%                 subplot(3,1,1);imagesc(z);xlabel('Data');
%                 subplot(3,1,2);imagesc(field1);clim([0 1]);xlabel('Field');
%                 subplot(3,1,3);imagesc(field2);clim([0 1]);ylabel(tc);xlabel('Subfield');
%             end
%         end
        fieldSize = sum(field(:));
        % Keep this field if its size is sufficient
        if fieldSize > minSize
            mapStats{unit,1}{c}.field(:,i) = field;
            mapStats{unit,1}{c}.size(i) = fieldSize;
            mapStats{unit,1}{c}.peak(i) = peak;
            mapStats{unit,1}{c}.mean(i) = mean(z(field));
            idx = find(field & z == peak);
            [mapStats{unit,1}{c}.y(i),mapStats{unit,1}{c}.x(i)] = ind2sub(size(z),idx(1));
            [x,y] = FieldBoundaries(field,circX,circY);
            [mapStats{unit,1}{c}.fieldX(i,:),mapStats{unit,1}{c}.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
            i = i + 1;
        end
        % Mark field bins as visited
        z(field) = NaN;
        if all(isnan(z)), break; end
    end
    
    end
end

%% Refine place fields 
    %%% this is not working yet
secondaryPF = 0.6; %secondary PF should have peak FR hihger than secondayPF*peakFR of primary PF

for unit = 1:length(firingMaps.rateMaps)
    for c = 1:length(firingMaps.rateMaps{1})
        if mapStats{unit}{c}.peak ~= 0
           for k=1:length(stats{unit}{1}.peak) % for each PF of this cell
               if mapStats{unit}{c}.peak(k)>mapStats{unit}{c}.peak(1)*secondaryPF % destroy small secondary PFs
                   peakLoc(k)=firingMaps.rateMaps{unit}{c}.x(mapStats{unit}{c}.x(k));
                   scatter(peakLoc(k),mapStats{unit}{c}.peak(k),'MarkerFaceColor','r','LineWidth',1);
                   fieldIni(k)=firingMaps.rateMaps{unit}{c}.x(mapStats{unit}{c}.fieldX(k,1));
                   fieldFin(k)=firingMaps.rateMaps{unit}{c}.x(mapStats{unit}{c}.fieldX(k,2));
                   x1=fieldIni(k);y1=get(gca,'ylim');hold on; plot([x1 x1],y1,'-.k','LineWidth',1);
                   x2=fieldFin(k);y1=get(gca,'ylim');hold on; plot([x2 x2],y1,'-.k','LineWidth',1);
               else
                   mapStats{unit}{c}.x(k) = NaN;      mapStats{unit}{c}.y(k) = NaN; 
                   mapStats{unit}{c}.field(k) = 0;    mapStats{unit}{c}.size(k) = 0; 
                   mapStats{unit}{c}.peak(k) = 0;     mapStats{unit}{c}.mean(k) = 0;  
                   mapStats{unit}{c}.fieldX(k,1:2) = [NaN,NaN]; mapStats{unit}{c}.fieldY(k,1:2) = [NaN,NaN];                 
               end
           end
           for k=1:length(mapStats{unit}{c}.peak)% destroy too large PF (CHANGE for recalculate!!)
               if mapStats{unit}{c}.peak(k)~=0 && ... 
                       firingMaps.rateMaps{unit}{c}.x(stats{i}{1}.fieldX(k,2))-firingMaps.rateMaps{unit}{c}.x(mapStats{unit}{c}.fieldX(k,1)) > 0.5
                   mapStats{unit}{c}.x(k) = NaN;      mapStats{unit}{c}.y(k) = NaN; 
                   mapStats{unit}{c}.field(k) = 0;    mapStats{unit}{c}.size(k) = 0; 
                   mapStats{unit}{c}.peak(k) = 0;     mapStats{unit}{c}.mean(k) = 0;  
                   mapStats{unit}{c}.fieldX(k,1:2) = [NaN,NaN]; mapStats{unit}{c}.fieldY(k,1:2) = [NaN,NaN];   
               end     
           end
           clear temp peakLoc fieldIni fieldFin x1 x2 y1 y2 k
        end
    end
end


%% Plot
for unit = 1:length(firingMaps.rateMaps)
    figure;
    for c = 1:length(firingMaps.rateMaps{1})
        subplot(2,2,c);plot(firingMaps.rateMaps{unit}{c});
        ylabel('FR(Hz)');xlabel('track (cm)');hold on;
        
             suptitle([basename ' cell' num2str(i)]);
             saveas(gcf,[pwd '\newPCs\cell_' num2str(i) '.png'],'png'); close all;            
        
        
    end
end

 %% Write output   
    
 
end


%%
% ------------------------------- Helper functions -------------------------------

% Field boundaries (circumscribed rectangle)

function [x,y] = FieldBoundaries(field,circX,circY)

% Find boundaries
x = find(any(field,1));
if isempty(x),
	x = [NaN NaN];
else
	x = [x(1) x(end)];
end
y = find(any(field,2));
if isempty(y),
	y = [NaN NaN];
else
	y = [y(1) y(end)];
end

% The above works in almost all cases; it fails however for circular coordinates if the field extends
% around an edge, e.g. for angles between 350° and 30°

if circX && x(1) == 1 && x(2) == size(field,2),
	xx = find(~all(field,1));
	if ~isempty(xx),
		x = [xx(end) xx(1)];
	end
end
if circY && y(1) == 1 && y(2) == size(field,1),
	yy = find(~all(field,2));
	if ~isempty(yy),
		y = [yy(end) yy(1)];
	end
end
end
