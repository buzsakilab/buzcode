
function [armChoice] = getArmChoice(varargin)
% Compute T maze performance and gets timestamps
%
% USAGE
%
%   [armChoice, behavior] = getArmChoice(varargin)
%
% INPUTS
% 
% digitalIn                     digitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left Alternation,  4.Righ Alternation
%                                 5. Home Delay,        6. Is alternation forzed?
%                                 7. 0 is cue left, 1 is cue right (for
%                                       cueSide task)
%                               If not called, look for it.
% task                          'alternation' and 'cudeSide'
% force                         Force detection (boolean, default false)
% verbose                       default false
%
% OUTPUT
%       - armChoice.behaviour output structure, with the fields:
% armChoice.timestamps          Choice timestamps, in seconds
% armChoice.visitedArm                 Choosed arm, 0 is left, 1 is right
% armChoice.delay.timestamps    Delay intervals, in seconds
% armChoice.delay.dur           Delay duration, in seconds
% armChoice.choice              Performance vector, 1 is right choice, 0 is
%                                   wrong. First choice is Nan.
% armChoice.performance         Alternation probability (#alternation/#trials)
% armChoice.forzed              1 if forzed alternation, 0 if spontaneous
%                                   alternation
% armChoice.task                'alternation' and 'cudeSide'  
% armChoice.expectedArm          In 'cueSide', it specifias the right choice. 
%
%  Manuel Valero 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'digitalIn',[],@ischar);
addParameter(p,'task',[]);
addParameter(p,'force',false,@islogical)
addParameter(p,'verbose',false,@islogical)
parse(p,varargin{:});
digitalIn = p.Results.digitalIn;
task = p.Results.task;
force = p.Results.force;
verbose = p.Results.verbose;

if ~isempty(dir('*.ArmChoice.Events.mat')) && ~force 
    disp('Arm choice already computed! Loading file.');
    file = dir('*.ArmChoice.Events.mat');
    load(file.name);
    return
end
%% Get basler TTL
disp('Loading digital In...');
if isempty(digitalIn)
    digitalIn = bz_getDigitalIn;
    if isempty(digitalIn)
        armChoice = [];
        return
    end
end
if isempty(task)
    switch size(digitalIn.timestampsOn,2)
        case 5
            task = 'alternation';
        case 6
            task = 'alternation';
        case 7
            task = 'cueSide';
    end
end
if strcmpi(task,'cueSide')
    if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
        armChoice.timestamps = [digitalIn.timestampsOn{3}; digitalIn.timestampsOn{4}]; 
        % 0 is left, 1 is right
        armChoice.visitedArm = [zeros(size(digitalIn.timestampsOn{3})); ones(size(digitalIn.timestampsOn{4}))];
        [armChoice.timestamps, idx] = sort(armChoice.timestamps);
        armChoice.visitedArm = armChoice.visitedArm(idx);
        % armChoice.delay.ints = digitalIn.ints{5};
        armChoice.delay.timestamps = digitalIn.ints{5};
        armChoice.delay.dur = nanmean(armChoice.delay.timestamps(2,:) - armChoice.delay.timestamps(1,:));
        
        % reconstruct side vector
        sideVector = zeros(1,int32(max(armChoice.timestamps)));
        onPulses = digitalIn.timestampsOn{7};
        onPulses(:,2) = nan;
        offPulses = digitalIn.timestampsOff{7};
        if  offPulses(1) < onPulses(1,1)
            offPulses(1) = [];
        end
        onPulses(1:length(offPulses),2) = offPulses;
        onPulses(isnan(onPulses)) = armChoice.timestamps(end) + 10;      
      
        for ii = 1:size(onPulses,1)
            sideVector(round(onPulses(ii,1)):round(onPulses(ii,2))) = 1;
        end
        armChoice.expectedArm = sideVector(round(armChoice.timestamps))'; % real choice;
        armChoice.choice = sideVector(round(armChoice.timestamps))' == armChoice.visitedArm; % 1 is right, 0 is wrong
        armChoice.performance = nansum(armChoice.choice)/(length(armChoice.choice));

        if size(digitalIn.timestampsOn,2) >=6
            if digitalIn.timestampsOn{6}>digitalIn.timestampsOff{6}
                armChoice.forzed = 1;
                desc = 'forzed side';
            else
                armChoice.forzed = 0;
                desc = 'spontaneous side';
            end
        else
            if armChoice.performance == 1 
                armChoice.forzed = 1;
                desc = 'forzed side';
            else
                armChoice.forzed = 0;
                desc = 'spontaneous side';
            end
        end

        h = figure;
        subplot(2,2,[1 2])
        hold on
        plot(armChoice.timestamps, armChoice.expectedArm,'color',[.7 .7 .7]);
        scatter(armChoice.timestamps(find(armChoice.choice == 0)),...
            armChoice.expectedArm(find(armChoice.choice == 0)),100,[.9 .6 .7],'filled','MarkerFaceAlpha',.5);
        scatter(armChoice.timestamps(find(armChoice.choice == 1)),...
            armChoice.expectedArm(find(armChoice.choice == 1)),100,[.6 .9 .7],'filled','MarkerFaceAlpha',.5);
        for ii = 1:size(armChoice.delay.timestamps,1)
            fill([armChoice.delay.timestamps(:,ii); flip(armChoice.delay.timestamps(:,ii))],[1 1 1.2 1.2]',...
            [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
        end
        xlabel('seconds'); ylim([-.2 1.2]);
        text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
            desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
            num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
            's'));
        set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});
        subplot(2,2,3)
        plot(cumsum(~armChoice.choice));
        set(gca,'TickDir','out'); xlabel('Trials'); ylabel('Errors');
 
        subplot(2,2,4)
        hold on
        for ii = 1:length(armChoice.choice)- 9
            lcurve(ii)=sum(armChoice.choice(ii:ii+9))/10;
        end
        lcurve_m = [];
        g = 1:10:length(armChoice.choice);
        for ii = 1:length(g)
            if g(ii) + 10 <= length(armChoice.choice)
                lcurve_m(ii) = mean(armChoice.choice(g(ii):g(ii)+10));
            else
                lcurve_m(ii) = mean(armChoice.choice(g(ii):length(armChoice.choice)));
            end
        end
        ecurve = cumsum(armChoice.choice)./(1:length(armChoice.choice))';
        area(1:length(armChoice.choice),ecurve,'FaceColor',[.8 .5 1],'EdgeColor','none');
        scatter(1:10:length(armChoice.choice),lcurve_m,100,[.9 .6 .7],'filled','MarkerFaceAlpha',.5);
        ylim([0 1]); xlabel('Trials'); ylabel('Performance'); 
        
        armChoice.lcurve_m = lcurve_m;
        armChoice.task = task;

        mkdir('Behavior');
        saveas(h,'Behavior\armChoice.png');
        if ~verbose
            close(h);
        end

        C = strsplit(pwd,'\');
        save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
    else
        warning('DigitalIn format does not match. Was Cue Side Maze performed? ');
    end
elseif strcmpi(task,'alternation')
    % score for alternation task
    if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
        armChoice.timestamps = [digitalIn.timestampsOn{3}; digitalIn.timestampsOn{4}]; 
        % 0 is left, 1 is right
        armChoice.visitedArm = [zeros(size(digitalIn.timestampsOn{3})); ones(size(digitalIn.timestampsOn{4}))];
        [armChoice.timestamps, idx] = sort(armChoice.timestamps);
        armChoice.visitedArm = armChoice.visitedArm(idx);
        % armChoice.delay.ints = digitalIn.ints{5};
        armChoice.delay.timestamps = digitalIn.ints{5};
        armChoice.delay.dur = nanmean(armChoice.delay.timestamps(2,:) - armChoice.delay.timestamps(1,:));
        armChoice.choice = [NaN; abs(diff(armChoice.visitedArm))]; % 1 is right, 0 is wrong
        armChoice.performance = nansum(armChoice.choice)/(length(armChoice.choice) - 1);
        armChoice.task = task;
        armChoice.expectedArm = [NaN; ~xor(armChoice.visitedArm(2:end), armChoice.choice(2:end))];

        if size(digitalIn.timestampsOn,2) >=6
            if digitalIn.timestampsOn{6}>digitalIn.timestampsOff{6}
                armChoice.forzed = 1;
                desc = 'forzed alternation';
            else
                armChoice.forzed = 0;
                desc = 'spontaneous alternation';
            end
        else
            if armChoice.performance == 1 
                armChoice.forzed = 1;
                desc = 'forzed alternation';
            else
                armChoice.forzed = 0;
                desc = 'spontaneous alternation';
            end
        end

        h = figure;
        hold on
        plot(armChoice.timestamps, armChoice.visitedArm,'color',[.7 .7 .7]);
        scatter(armChoice.timestamps(isnan(armChoice.choice)),...
            armChoice.visitedArm(isnan(armChoice.choice)),100,[.8 .8 .8],'filled');
        scatter(armChoice.timestamps(find(armChoice.choice == 1)),...
            armChoice.visitedArm(find(armChoice.choice == 1)),100,[.6 .9 .7],'filled');
        scatter(armChoice.timestamps(find(armChoice.choice == 0)),...
            armChoice.visitedArm(find(armChoice.choice == 0)),100,[.9 .6 .7],'filled');
        for ii = 1:size(armChoice.delay.timestamps,1)
            fill([armChoice.delay.timestamps(:,ii); flip(armChoice.delay.timestamps(:,ii))],[1 1 1.2 1.2]',...
            [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
        end
        xlabel('seconds'); ylim([-.2 1.2]);
        text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
            desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
            num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
            's'));
        set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});

        mkdir('Behavior');
        saveas(h,'Behavior\armChoice.png');

        C = strsplit(pwd,'\');
        save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
    else
        warning('DigitalIn format does not match. Was T maze performed? ');
    end
    
end

end

%     h1 = figure;
%     subplot(3,1,[1 2])
%     hold on
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     scatter(mazeVirtual(:,1), mazeVirtual(:,2),5,linearizeMazeVirtual);
%     ylabel('cm'); 
%     subplot(3,1,[3])
%     hold on
%     scatter(linspace(0,max(vlinMaze)+50,length(vlinMaze)),...
%         vlinMaze,5,vlinMaze);
%     scatter(linspace(0,max(vlinMaze)+50,length(vlinMazeCont)),...
%         vlinMazeCont,5,vlinMazeCont);
%     xlabel('cm');  ylabel('cm');
%     colormap jet
%     saveas(h1,'Behavior\virtualMaze.png');

%     centroidR = [30, 22]; 
%     centroidL = [30, 55];
%     x = behavior.positions.x';
%     y = behavior.positions.y';
%     t = behavior.timestamps;
%     [rangle,rradius] = cart2pol(x-centroidR(1), y - centroidR(2));
%     [langle,lradius] = cart2pol(x-centroidL(1), y - centroidL(2));
%     rangle = wrapTo2Pi(rangle - 0.9458 + .15);
%     langle = wrapTo2Pi(-langle + 5.3658 + .15);
%     
%     rangle = (rad2deg(rangle))/(360/lenghTrial); % +180 to go from 0 to 360
%     langle = (rad2deg(langle))/(360/lenghTrial) + lenghTrial; % Distance traveled in cm, left trial after righ trial
%     
%     h1 = figure;
%     hold on
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     xlabel('cm'); ylabel('cm'); 
%     plot(centroidR(1),centroidR(2),'o','MarkerFaceColor',[.8 .4 .3],'MarkerEdgeColor','none');
%     plot(centroidL(1),centroidL(2),'o','MarkerFaceColor',[.8 .4 .3],'MarkerEdgeColor','none');
%     saveas(h1,'Behavior\centroids.png');
%     
%     % 
%     h2 = figure;
%     subplot(3,1,[1 2])
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     prev = 0;
%     linearized = []; arm = [];
%     for ii = 1:size(armChoice.delay.timestamps,1)
%         winTrial = [prev armChoice.delay.timestamps(ii)];
%         prev = armChoice.delay.timestamps(ii);
%         xspam = behavior.timestamps >= winTrial(1) & behavior.timestamps <= winTrial(2);
%         subplot(3,1,[1 2])
%         hold on
%         p = plot(x(xspam),y(xspam),'lineWidth',2);
%         ylabel('cm'); 
%         
%         if armChoice.arm(ii) == 1
%             linearized = [linearized rangle(xspam)];
%             arm = [arm ones(size(rangle(xspam)))];
%         elseif armChoice.arm(ii) == 0
%             linearized = [linearized langle(xspam)];
%             arm = [arm zeros(size(rangle(xspam)))]; 
%         end
%         subplot(3,1,3)
%         hold on
%         plot(t(xspam),linearized(xspam),'lineWidth',2);
%         xlabel('samples');
%         drawnow;
%         pause;    
%         delete(p);
%     end
%     % close(h2);
%       
% %     for ii = 1:size(behavior.timestamps)
% %     end    
% while 1
%     [xg,yg]=ginput(1);
%     fprintf('x: %3.2f, y: %3.2f \n',xg,yg);
%     
%     [rangle,rradius] = cart2pol(xg-centroidR(1), yg - centroidR(2));
%     [langle,lradius] = cart2pol(xg-centroidL(1), yg - centroidL(2));
%     rangle = wrapTo2Pi(rangle - 0.9458 + 0.15);
%     langle = wrapTo2Pi(-langle + 5.3658 + 0.15);
%     
%     rangle = (rad2deg(rangle))/(360/lenghTrial); % +180 to go from 0 to 360
%     langle = (rad2deg(langle))/(360/lenghTrial) + lenghTrial ;
%     fprintf('linearized R: %3.2f \n',rangle);
%     fprintf('linearized L: %3.2f \n\n',langle);
% end
