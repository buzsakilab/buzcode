function funcsynapses = FindSynapseWrapper_NOTSToolbox(basename)

funcsynapses = Make_FindSynapse_bw_NOTSToolbox(basename);%a wrapper which uses the usual lab CCG.m rather than yours and stores away a few parameters into a struct
funcsynapses = Make_FindZeroLagCorr_bw_NOTSToolbox(basename,funcsynapses);%to find zero-lag interactions

funcsynapses = FindSynapseToStruct(funcsynapses);%interprets FindSynapse.m output to classify connections and cells as E or I

FindSynapse_ReviewOutput(funcsynapses, 'funcsynapses');% a gui-ish means to allow users to review all connections and nominate bad ones with clicks on axes
% 
% h = figure('Visible','Off','Name','BWWaitFig');
% waitfor(h,'Name','DELETEMENOW')
% delete(h)
% 
% FindSynapse_ReviewZeroAndWide(funcsynapses, 'funcsynapses');% a gui-ish means to allow users to review all connections and nominate bad ones with clicks on axes








% 
% [T,G] = oneSeries(tsdspikes);
% T = TimePoints(T);
% 
% % 
% % for a = 1:size(tsdspikes,1)
% %     T = cat(1,TimePoints(tsdspikes{a});
% %     T = cat(T,res 
% % end
% % 
% % res2 = TimePoints(tsdspikes{clus2});
% % T = cat(1,res1,res2);
% % 
% % G = cat(1,1+zeros(size(res1)),2+zeros(size(res2)));
% % 
% SampleRate = 10000;%%Based on usual TSD object sample rate<< NEED TO GENERALIZE THIS, HAVEN'T FIGURED OUT HOW YET
% numBinsBinSize = BinSize*SampleRate/1000;FindSynapse_ReviewOutput(funcsynapses, 'funcsynapses')
% HalfBins = round(300/numBinsBinSize);
% 
% [ccgR, tR, Pairs] = CCG(T, G, numBinsBinSize, HalfBins, SampleRate, unique(G), 'count'); %calc cross correlograms, output as counts... 3D output array
% 
% numpairs = size(tsdspikes,1) * (size(tsdspikes,1)-1)/2;%% reshape ccgs to a series of columnar ccgs to pass into FindSynapse.m
% ccgs = zeros(length(tR),numpairs);
% ccgindices = zeros(numpairs,2);
% counter = 0;
% for a = 1:size(tsdspikes,1)
%     for b = a+1:size(tsdspikes,1)
%         counter = counter+1;
%         ccgs(:,counter) = ccgR(:,a,b);
%         ccgindices(counter,:) = [a,b];
%     end
% end
% % [excPairs,inhPairs,pkTimes,trTimes,hiBoundAll,loBoundAll] = FindSynapse...
% %     (ccgs,'bins',BinSize,'alpha',0.01,'synwin',[1 4],'convWin',12);
% 
% [synLat,synStrZ,synStrR,sigBounds] = FindSynapse(ccgs,'bins',BinSize,'alpha',0.01,'synwin',[1 4],'convWin',12);
% [synLatZero,synStrZZero,synStrRZero,sigBoundsZero] = FindSynapse(ccgs,'bins',BinSize,'alpha',0.01,'synwin',[-.5 .5],'convWin',12);
% 
% 
% 
% 
% %review cell by cell
% FindSynapse_ReviewOutput(tR, ccgs, ccgindices, numpairs, cellIx, shank,excPairs,inhPairs,pkTimes,trTimes,hiBoundAll,loBoundAll, 'SynapseStruct')

