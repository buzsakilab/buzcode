function [  ] = EE_Initiate( FO )
%This function initializes the EventExplorer figure


bz_GetLFP

FO.fig = figure;
set(FO.fig, 'numbertitle', 'off', 'name', ['Recording: ', FO.baseName,'. Events: ',FO.EventName]);
set(FO.fig, 'Tag', 'EventExplorerMaster');

%From StateEditor
% FO.fig = figure('KeyReleaseFcn', {@DefKey}, 'Position', posvar);
% set(FO.fig, 'numbertitle', 'off', 'name', ['Event: ', FO.baseName]);
% set(FO.fig,'WindowButtonDownFcn', {@MouseClick}, 'WindowButtonUpFcn', {@unMouseClick}, 'Units', 'normalized');
% set(FO.fig, 'WindowButtonMotionFcn', {@Nothing}, 'WindowScrollWheelFcn', {@MouseScroll});
% set(FO.fig, 'CloseRequestFcn', {@CloseDialog});
% set(FO.fig, 'Tag', 'EventExplorerMaster');


%Store the data in the figure
guidata(FO.fig, FO);
end

