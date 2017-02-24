function [redraw, rekey, undoable] = ShowAllWaveforms

% ShowAllWaveforms
%  Shows the waveforms of clusters specified by user in a popup window.
%
% INPUTS: None (user interface)
%
% OUTPUTS: None (modifies MClust global variables)
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/GeneralizedCutterOptions folder

% JCJ Sept 2007

redraw = false; rekey = false; undoable = false;

global MClust_Clusters MClust_Colors MClust_Hide MClust_TTData
ShrinkSize=.4;

% User input
prompt={'Cluster(s) to show (comma or space separated list):','Show Hides:','Show Shows:','Shrink Ratio (compared to standard figure size)'};
def={'','no','no',num2str(ShrinkSize)};
dlgTitle='Show waveforms of cluster(s)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);

if ~isempty(answer)
    ShrinkSize=str2num(answer{4});

    iClust1=str2num(answer{1});

    iClust2=[];
    if strncmpi(answer{2},'yes',1) | str2num(answer{2})==1
        nClust=length(MClust_Clusters);
        iClust2=find(MClust_Hide(2:(nClust+1)));
    end

    iClust3=[];
    if strncmpi(answer{3},'yes',1) | str2num(answer{3})==1
        nClust=length(MClust_Clusters);
        iClust3=find(~MClust_Hide(2:(nClust+1)));
    end

    iClust2Show=unique([shiftdim(iClust1);iClust2;iClust3]);

    if ~isempty(iClust2Show)
        nClust2Show=length(iClust2Show);
        ShowWVFig=[];
        for iD=1:nClust2Show
            iClust=iClust2Show(iD);
            [f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});
            if isempty(f) 
                msgbox('No points in cluster.')
            else
                [clustTT was_scaled] = ExtractCluster(MClust_TTData, f);

                WV = clustTT;
                Color = MClust_Colors(iClust + 1,:);

                WVD = Data(WV);
                
                % Calculate waveform means
                mWV = squeeze(mean(WVD,1));

                if length(WVD(:,1,1)) > 1000
                    maxPts = 1000;
                    Pts = randperm(size(WVD,1));
                    Pts = Pts(1:1000);
                else
                    maxPts = size(WVD,1);
                    Pts = 1:size(WVD,1);
                end

                ShowWVFig(end+1) = figure;
                for it = 1:4
                    xrange = ((size(WVD,3) + 2) * (it-1)) + (1:size(WVD,3));
                    figure(ShowWVFig(end));
                    hold on;
                    plot(xrange,squeeze(WVD(Pts,it,:))','color',Color);
                    plot(xrange, mWV(it,:),'w');
                end

                axis off
                axis([0 4*(size(WVD,3) + 2) -2100 2100])
                title({['Cluster ' num2str(iClust) ' Waveforms']; [num2str(maxPts) ' of ' num2str(length(WVD(:,1,1))) ' waveforms shown']},'FontSize',14*ShrinkSize);
                set(ShowWVFig(end),'Name',['Cluster ' num2str(iClust) ])
                hold off
            end
        end
        ShrinkFig(ShowWVFig,ShrinkSize)
    end

end
