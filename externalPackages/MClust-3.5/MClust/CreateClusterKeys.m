function CreateClusterKeys(figHandle, iC, XLoc, YLoc, CallbackFunc, varargin)

% CreateClusterKeys(figHandle, iC, XLoc, YLoc, CallbackFunc, varargin)
% 
%
% INPUTS
%     figHandle - figure to draw on
%     iC - cluster number
%     XLoc - starting x location
%     YLoc - y location
%     Callback function 
%     varargin - Strings/Tags to add pushbuttons at end
%
% ADR 1999
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


width = 0.10;
height = 0.05; 

% color
global MClust_Colors MClust_Hide MClust_UnaccountedForOnly
global MClust_Clusters

uicontrol('Parent', figHandle, ...
   'Units', 'Normalized', 'Position', [XLoc YLoc width/2 height], ...
   'Style', 'frame', 'BackgroundColor', MClust_Colors(iC+1,:), 'Enable', 'inactive', ...
   'Tag', 'ChooseColor', 'UserData', iC, 'ButtonDownFcn', CallbackFunc, ...
   'TooltipString', 'Change cluster color');

uicontrol('Parent', figHandle, ...
   'Units', 'Normalized', 'Position', [XLoc+width/2 YLoc width/2 height], ...
   'Style', 'text', 'String', num2str(iC), 'UserData', iC);

% show/hide
uicontrol('Parent', figHandle, ...
   'Units', 'Normalized', 'Position', [XLoc+width YLoc width height], ...
   'Style', 'checkbox', 'Tag', 'HideCluster', 'String', 'Hide', 'Value', MClust_Hide(iC+1), 'UserData', iC, ...
   'Callback', CallbackFunc, ...
   'TooltipString', 'Hide/Show cluster');

if iC ==0
	uicontrol('Parent', figHandle, ...
		'Units', 'Normalized', 'Position', [XLoc+2*width YLoc 3*width height],...
		'Style', 'checkbox', 'Tag', 'UnaccountedForOnly', 'String', 'Unaccounted for points only', ...
		'UserData', iC, 'Callback', CallbackFunc, 'Value', MClust_UnaccountedForOnly, ...
		'TooltipString', 'If checked, only show points not in any other cluster');
else % iC > 0
	if which(fullfile(class(MClust_Clusters{iC}), 'GetName')) % ADR 2008 getName correction
		Name = GetName(MClust_Clusters{iC});
	else
		Name = 'functions';
	end
	uicontrol('Parent', figHandle, ...
		'Units', 'Normalized', 'Position', [XLoc+2*width YLoc 3*width height],...
		'Style', 'popupmenu', 'Tag', 'ClusterFunctions', 'String', cat(2,Name,varargin), ...
		'UserData', iC, 'Callback', CallbackFunc);
end
