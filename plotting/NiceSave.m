function [  ] = NiceSave(figname,figfolder,baseName,figtype)
%NiceSave(figname,baseName,recname) formats the figure for best viewing
%and saves it as a .pdf in figfolder with name recname_figname.pdf
%
%DLevenstein Fall 2016
%%

if ~exist('figtype','var')
    figtype = 'pdf';
end

if ~exist(figfolder,'dir')
    mkdir(figfolder)
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition',[0 0 1 1]);
saveas(gcf,[figfolder,'/',baseName,'_',figname,'.',figtype],figtype) ;

end

