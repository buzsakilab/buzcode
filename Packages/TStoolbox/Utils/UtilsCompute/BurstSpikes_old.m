function [burstCell,burstLim,Burst] = BurstSpikes(S)

    dispFig=0; 
    Burst = cell(length(S),1);
    burstCell = zeros(length(S),1);
    burstLim = zeros(length(S),1);
    
    disp('Cell #:')

    for c=1:length(S)              
        
        fprintf([num2str(c) ' ']);
        if ~mod(c,10)
            fprintf('\n');
        end
        
        rg = Range(S{c});
        dr = diff(rg);
        
        drl = log10(dr/10000);
        drl = drl(drl<3);
        
        [p htab] = MultimodalTest(drl,2);
        keyboard
        if p<0.05

            burstCell(c)=1;
            
            if nargout>2
            b = min(drl):htab(1,2)/2:max(drl);
            h = hist(drl,b);
            
            mxModes = LocalMinima(-h,1,max(-h));
            
            if length(mxModes)~=2
                h1 = h(mxModes);
                [dummy,ix] = sort(h(mxModes),'descend');
                mxModes = [mxModes(ix(2)) mxModes(ix(1))];
            end
             
            h1 = h(mxModes(1):mxModes(2));
            b1 = b(mxModes(1):mxModes(2));
            
            mnModes = LocalMinima(h1,1,max(h1));
            mnModes = mnModes(1);
            
            minVal = b1(mnModes);
            burstLim(c) = 10^(minVal+4);

            b = min(drl):(max(drl)-min(drl))/41:max(drl);
            if dispFig==1
                figure(1),clf
                    h = hist(drl(drl>minVal),b);
                    bar(b,h,1,'EdgeColor','k','FaceColor','k')
                    hold on
                    col = [0.4 0.4 0.4];
                    h = hist(drl(drl<=minVal),b);
                    bar(b,h,1,'EdgeColor',col,'FaceColor',col)
                    yl = ylim;
                    line([minVal minVal],[0 yl(2)],'Color','r','LineWidth',2)
                    xlabel('log ISI')
                    ylabel('# spikes')
            end
            end
        else
           if dispFig==1 
           b = min(drl):(max(drl)-min(drl))/41:max(drl);

            minVal = log10(0.04);
            
             figure(1),clf
             h = hist(drl,b);
             bar(b,h,1,'EdgeColor','k','FaceColor','k')
                xlabel('log ISI')
                ylabel('# spikes')
           end
        end
        
        if 0
        set(1,'PaperPositionMode','auto')
        name = [figDir 'ISI_Cell_' num2str(c)];
        print('-f1','-painters','-depsc2',[name '.eps']);
        %print('-f1','-painters','-dpng',[name '.png']);
        end
        
        if nargout==3
            shortDt = dr<10^(minVal+4);	

            spkBurstIx = zeros(length(rg),1);
            firstSpkBurst = [];
            burstLg = [];

            u=1;
            while u<length(shortDt)
              if shortDt(u)
                firstSpkBurst = [firstSpkBurst;rg(u)];
                spkBurstIx(u) = 1;
                spkBurstIx(u+1) = 2;
                v=1;
                while shortDt(u+v) && (u+v)<length(shortDt)
                  spkBurstIx(u+v+1) = v+2;
                  v=v+1;
                end
                burstLg = [burstLg;v];
                u = u+v;
              else
                spkBurstIx(u) = -1;
              end
              u = u+1;
            end

            Burst{c} = spkBurstIx;
        end
        
    end
    
    
    fprintf('\n');

end
    