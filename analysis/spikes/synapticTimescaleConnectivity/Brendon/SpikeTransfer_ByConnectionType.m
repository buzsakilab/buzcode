function SpikeTransfer_ByConnectionType(funcsynapses,S)

% ccg = funcsynapses.fullCCGMtx;

for a=1:size(funcsynapses.ConnectionsE,1)
    pre = funcsynapses.ConnectionsE(a,1);
    post = funcsynapses.ConnectionsE(a,2);
    
%     thresh = funcsynapses.PairUpperThreshs(max([pre post]),min([pre post]));
%     if funcsynapses.CellShanks(pre)==funcsynapses.CellShanks(post)
%         
%     else
    cstart = funcsynapses.CnxnStartTimesVsRefSpk(pre,post); %start time
    cend = funcsynapses.CnxnEndTimesVsRefSpk(pre,post);%end time
    if length(cstart) == length(cend) & ~isempty(cstart)
        EStrengths(a) = SpikeTransfer(TimePoints(S{pre}),TimePoints(S{post}),funcsynapses.BinMs,[cstart cend]);
%         figure;bar(funcsynapses.CCGbins,ccg(:,pre,post))
%         hold on;plot(cnxntimes,[mean(ccg(:,pre,post)) mean(ccg(:,pre,post))],'m','LineWidth',5)
    end
end

for a=1:size(funcsynapses.ConnectionsI,1)
    pre = funcsynapses.ConnectionsI(a,1);
    post = funcsynapses.ConnectionsI(a,2);
    
%     thresh = funcsynapses.PairUpperThreshs(max([pre post]),min([pre post]));
%     if funcsynapses.CellShanks(pre)==funcsynapses.CellShanks(post)
%         
%     else
    cstart = funcsynapses.CnxnStartTimesVsRefSpk(pre,post); %start time
    cend = funcsynapses.CnxnEndTimesVsRefSpk(pre,post);%end time
    if length(cstart) == length(cend) & ~isempty(cstart)
        IStrengths(a) = SpikeTransfer(TimePoints(S{pre}),TimePoints(S{post}),funcsynapses.BinMs,[cstart cend]);
%         figure;bar(funcsynapses.CCGbins,ccg(:,pre,post))
%         hold on;plot(cnxntimes,[mean(ccg(:,pre,post)) mean(ccg(:,pre,post))],'m','LineWidth',5)
    end
end

%Go through each zerolag
1;
q