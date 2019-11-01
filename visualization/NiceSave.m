function [  ] = NiceSave(figname,figfolder,baseName,varargin)
%NiceSave(figname,figfolder,baseName,figtype) formats the figure for best viewing
%and saves it as a .pdf in figfolder with name recname_figname.pdf
%
%INPUTS
%   figname     string with title of figure
%   figfolder   folder to save figure to
%   baseName    name of the recording ([] for none)
%   'figtype'     'pdf','jpeg',etc
%   'includeDate' default: false. include date in filename
%
%DLevenstein Fall 2016
%TODO: add option to not overwrite old figure, but instead save as
%'figname02' etc.
%%
p = inputParser;
addParameter(p,'figtype','pdf')
addParameter(p,'includeDate',false)
parse(p,varargin{:})
figtype = p.Results.figtype;
includeDate = p.Results.includeDate;


%%

if ~exist(figfolder,'dir')
    mkdir(figfolder)
end

if includeDate
   figname = [figname,'_',date];
end
%set(gcf,'TickDir','out')

%if strcmp(figtype,'jpg')
set(gcf,'PaperOrientation','landscape');
%end
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition',[0 0 1 1]);
%orient portrait
saveas(gcf,[figfolder,'/',baseName,'_',figname,'.',figtype],figtype) ;

end

