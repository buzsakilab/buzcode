
function [classEvs] = checkEvents(ts)
% check events
ch = 'all';
basepath = pwd;
noPrompts = true;
eventsType = 'UDStates';
win = 5;
shanks = 1;
sw = [];

sessionInfo = bz_getSessionInfo(pwd, 'noPrompts', noPrompts);
lfp = bz_GetLFP(ch,'basepath',basepath,'noPrompts',noPrompts);

classEvs = ones(size(ts));

fprintf('%3.i Events to check. \n',length(ts));
close all
fig=figure;
% set(gcf,'Position',[100 100 2000 1000])
ii = 1;
uicontrol('Style','pushbutton','Position',[515 370 40 25],'Units','normalized'...
            ,'String','>','Callback',@(src,evnt)advance);
uicontrol('Style','pushbutton','Position',[515 340 40 25],'Units','normalized'...
    ,'String','<','Callback',@(src,evnt)back);
uicontrol('Style','pushbutton','Position',[515 270 40 25],'Units','normalized','BackgroundColor',[.6 .8 .6]...
    ,'String','Good','Callback',@(src,evnt)buttonCallback('1'));
uicontrol('Style','pushbutton','Position',[515 240 40 25],'Units','normalized','BackgroundColor',[.8 .6 .6]...
    ,'String','Bad','Callback',@(src,evnt)buttonCallback('2'));
uicontrol('Style','pushbutton','Position',[515 210 40 25],'Units','normalized','BackgroundColor',[.8 .8 .8]...
    ,'String','NA','Callback',@(src,evnt)buttonCallback('0')); 
t1 = uicontrol('Style','text','Position',[515 180 40 25],'Units','normalized',...
    'String',strcat('1/',num2str(length(ts))),'Callback',[]);
t2 = uicontrol('Style','text','Position',[515 160 40 25],'Units','normalized','ForegroundColor',[.4 .8 .4]...
    ,'String',strcat(num2str(length(find(classEvs==1)))),'Callback',[]); 


while ii <= length(ts)
    winPlot = find(lfp.timestamps >= ts(ii) - win/2 & lfp.timestamps <= ts(ii) + win/2);
    for jj = shanks
    	subplot(1,length(shanks),jj);
        cla
        mm = 1;
        for kk = sessionInfo.spikeGroups.groups{jj} + 1 
            hold on
            plot(lfp.timestamps(winPlot),double(lfp.data(winPlot,kk))-(mm)*2000);
            mm = mm + 1;
            xlim([lfp.timestamps(winPlot(1)) lfp.timestamps(winPlot(end))]);
            set(gca, 'YTick',[]); ylabel('ch'); xlabel('s');
            ax = axis;
            plot([ts(ii) ts(ii)],[ax(3) ax(4)],'--r');
        end
        if classEvs(ii)==1
            title('Good','FontWeight','normal','Color',[.4 .8 .4]);
        elseif classEvs(ii)==2
            title('Bad','FontWeight','normal','Color',[.8 .4 .4]);
        elseif classEvs(ii)==0
            title('NA','FontWeight','normal','Color',[.6 .6 .6]);
        end
        t1.String = strcat(num2str(ii),'/',num2str(length(ts)));
        t2.String = strcat(num2str(length(find(classEvs==1))));
        uiwait(fig);
    end
end
close(fig);

save(strcat(date,'ChekEvents.mat'),'classEvs');

function buttonCallback(newString)
      classEvs(ii) = str2num(newString);
      ii = ii + 1;
      uiresume(fig);
end

function advance()
      ii = ii + 1;
      uiresume(fig);
end

function back()
      ii = ii - 1;
      if ii < 1
          ii = 1;
      end
      uiresume(fig);
end

end
