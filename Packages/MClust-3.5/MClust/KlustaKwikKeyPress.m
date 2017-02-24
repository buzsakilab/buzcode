function KlustaKwikKeyPress(c)

% KlustaKwikKeyPress
%
% Callbacks for KK decision window
%
% ADR 2006-2008
%
% Version control M3.5

%---------------------------
% get startup info

global KKC_IDSet

figHandle = get(0, 'PointerWindow');
KKFigHandle = findobj('Type', 'figure', 'Tag', 'KKDecisionWindow');  

if nargin == 0
    c = get(KKFigHandle, 'CurrentKey');
end

% find current visible page
for iC = 1:length(KKC_IDSet)
	if strcmpi(get(KKC_IDSet{iC}.ui(1), 'visible'),'on')
		break;
	end
end
minVisIC = iC;
currentPageNumber = KKC_IDSet{minVisIC}.page;
maxVisIC = iC;
for iC = minVisIC:length(KKC_IDSet)
	if (KKC_IDSet{iC}.page == currentPageNumber) && (iC > maxVisIC)
        maxVisIC = iC;
	end
end

switch c
	
	case {'k','K'}
		SelectedCluster = get(findobj(figHandle,'Tag','SelectCluster','BackgroundColor','c'),'UserData');
		keepHandle = findobj(figHandle,'Tag','KeepCluster','UserData', SelectedCluster);
		set(keepHandle, 'Value', ~get(keepHandle, 'Value'));
		KlustaKwikCallbacks('KeepCluster', keepHandle);
  
	case 'downarrow'
		SelectedCluster = get(findobj(figHandle,'Tag','SelectCluster','BackgroundColor','c'),'UserData');		
		if maxVisIC > SelectedCluster
			newClusterHandle = findobj(figHandle,'Tag','SelectCluster','UserData',SelectedCluster+1);
			KlustaKwikCallbacks('SelectCluster', newClusterHandle)
        else
            KlustaKwikKeyPress('pagedown');
		end

	case 'uparrow'
		SelectedCluster = get(findobj(figHandle,'Tag','SelectCluster','BackgroundColor','c'),'UserData');		
		if minVisIC < SelectedCluster
			newClusterHandle = findobj(figHandle,'Tag','SelectCluster','UserData',SelectedCluster-1);
			KlustaKwikCallbacks('SelectCluster', newClusterHandle)
        else
            KlustaKwikKeyPress('pageup');
        end

		
	case 'pageup'
        if currentPageNumber > 1
            currentPageNumber = currentPageNumber - 1;
            for iC = 1:length(KKC_IDSet)
                if KKC_IDSet{iC}.page == currentPageNumber
                    newClusterHandle = findobj(figHandle,'Tag','SelectCluster','UserData',iC);
                    for iK = 1:length(KKC_IDSet{iC}.ui)
                        set(KKC_IDSet{iC}.ui(iK), 'visible', 'on', 'enable', 'on');
                    end
                else
                    for iK = 1:length(KKC_IDSet{iC}.ui)
                        set(KKC_IDSet{iC}.ui(iK), 'visible', 'off', 'enable', 'off');
                    end
                end
            end
            KlustaKwikCallbacks('SelectCluster', newClusterHandle)
        end
        
    case 'pagedown'
        newClusterHandle = [];
        if currentPageNumber < KKC_IDSet{end}.page
            currentPageNumber = currentPageNumber +1;
            for iC = 1:length(KKC_IDSet)
                if KKC_IDSet{iC}.page == currentPageNumber
                    if isempty(newClusterHandle)
                        newClusterHandle = findobj(figHandle,'Tag','SelectCluster','UserData',iC);
                    end    
                    for iK = 1:length(KKC_IDSet{iC}.ui)
                        set(KKC_IDSet{iC}.ui(iK), 'visible', 'on', 'enable', 'on');
                    end

                else
                    for iK = 1:length(KKC_IDSet{iC}.ui)
                        set(KKC_IDSet{iC}.ui(iK), 'visible', 'off', 'enable', 'off');
                    end
                end
            end
            KlustaKwikCallbacks('SelectCluster', newClusterHandle)
        end
end




