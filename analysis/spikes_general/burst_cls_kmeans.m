function [bn,burstiness,mburstf,mibsn,mbl,mibfr,mibif,burstf,ibsn,burstlen,ibfr,ibif,bp,normbnbdiff,burst,burstspikes,nonburstspikes] = burst_cls_kmeans(segspikes,sro)

% K-means clustering of ISIs
tivs = cell(1,size(segspikes,2));
ivnum = zeros(1,size(segspikes,2));
segfirstiv = zeros(1,size(segspikes,2));
seglastiv = zeros(1,size(segspikes,2));
for b1 = 1 : size(segspikes,2)
    tivs{b1} = diff(segspikes{b1});
    ivnum(b1) = length(tivs{b1});
    if b1 == 1
        seglastiv(b1) = ivnum(b1);
    else
        seglastiv(b1) = sum(ivnum(1 : b1 - 1)) + ivnum(b1);
    end
    if b1 == 1
        segfirstiv(b1) = 1;
    else
        segfirstiv(b1) = seglastiv(b1 - 1) + 1;
    end
end
ivs = cat(2,tivs{:});

%k-means starts
if length(ivs) > 49
    replino = 5;
    maxclust = 5;
    qual = zeros(1,maxclust);
    clustid = zeros(length(ivs),maxclust);
    for kmc1 = 2 : maxclust
        clustid(:,kmc1) = kmeans(ivs,kmc1,'emptyaction','drop','replicates',replino);
        qual(kmc1) = mean(silhouette(ivs',clustid(:,kmc1)));
    end
    ici = find(qual == max(qual));
    idealclust = clustid(:,ici);
    isiclusters = cell(1, ici);
    isiclustmean = zeros(1,ici);
    for kmc2 = 1 : ici
        isiclusters{kmc2} = ivs(idealclust == kmc2);
        isiclustmean(kmc2) = mean(isiclusters{kmc2});
    end
    ibclust = isiclusters{isiclustmean == min(isiclustmean)};
%k-means ends

%bf
    if isempty(ibclust) == 0
        firstspikes = cell(1,size(segspikes,2));
        lastspikes = cell(1,size(segspikes,2));
        segburstivs = cell(1,size(segspikes,2));
        for li1 = 1 : size(segspikes,2)
            segburstivs{li1} = find(tivs{li1} <= max(ibclust));
            sbiv1 = diff(segburstivs{li1});
            sbiv2 = find(sbiv1 > 1);
            firstspikes{li1} = segspikes{li1}(segburstivs{li1}(sbiv2 + 1));
            lastspikes{li1} = segspikes{li1}(segburstivs{li1}(sbiv2) + 1);
            firstspikes{li1} = [segspikes{li1}(1) firstspikes{li1}];
            lastspikes{li1} = [lastspikes{li1} segspikes{li1}(end)];
        end
        burst = [cat(2,firstspikes{:});cat(2,lastspikes{:})];
        burstlen = (burst(2,:) - burst(1,:)) / (sro / 1000); %length of bursts in ms
        burst(:,burstlen > 400) = [];
        burstlen(burstlen > 400) = [];
%bf

% maxthetacycle = 400;
% nonbursts = burst(:,bli >= maxthetacycle);
% burst(:,bli >= maxthetacycle) = [];

%Selection of burst and non-burst spikes
        spktm = cat(2,segspikes{:});
        burstspikes = cell(1,size(burst,2));
        ibsn = zeros(1,size(burst,2));
        for bj = 1 : size(burst,2)
            burstspikes{bj} = spktm(spktm >= burst(1,bj) & spktm <= burst(2,bj));
            ibsn(bj) = length(burstspikes{bj});
        end
        burstspkvect = cat(2,burstspikes{:});
        nbs = zeros(1,length(spktm));
        for nbi = 1 : length(spktm)
            for nbj = 1 : length(burstspkvect)
                if spktm(nbi) == burstspkvect(nbj)
                nbs(nbi) = 1;
                end
            end
        end
        nonburstspikes = spktm(nbs == 0);
%selection of burst & non-burst intervals
        nbss = zeros(1,length(ivs));
        for ni = 1 : length(ivs)
            for nj = 1 : length(ibclust)
                if ivs(ni) == ibclust(nj)
                    nbss(ni) = 1;
                end
            end
        end
        nbivs = ivs(nbss ~= 1);

%Comparison of burst and non-burst intervals
        if length(nbivs) >= 5
            bp = ranksum(ibclust,nbivs,0.01,'tail','both'); %ranksum test for burst - non-burst intervals
        else
            bp = -1;
        end
        normbnbdiff = (median(nbivs) - median(ibclust)) / (median(nbivs) + median(ibclust)); %normalized diff of burst - non-burst intervals

%burst parameter definitions
        mbl = median(burstlen); %median burst length
        ibif = 1./ibclust * sro; %intraburst inst.frequency
        mibif = median(1./ibclust*sro); %median intraburst inst.frequency
        ibfr = (ibsn ./ (burstlen/1000)); %intraburst firing rate per burst
        mibfr = median(ibfr); %median intraburst firing rate
        burstiness = (size(burstspkvect,2) / size(cat(2,segspikes{:},2),2)) * 100; %proportion of burst spikes
        covibsi = std(ibclust) / mean(ibclust); %COV of intraburst intervals
        first_isis = diff(cat(2,firstspikes{:})); %inter-firstspike intervals
        mfisis = median(first_isis) / (sro / 1000); % median interburst interval (ms)
        burstf = (1./first_isis)*sro; %inter-firstspike IF        
        mburstf = median(burstf); %burst occurrence frequency
        bn = size(burst,2); %number of bursts
        mibsn = median(ibsn); %median intraburst spike number
    else
        bn = NaN;
        burstiness = NaN;
        mburstf = NaN;
        mibsn = NaN;
        mbl = NaN;
        mibfr = NaN;
        mibif = NaN;
        burstf = NaN;
        ibsn = NaN;
        burstlen = NaN;
        ibfr = NaN;
        ibif = NaN;
        bp = NaN;
        normbnbdiff = NaN;
        burst = NaN;
        burstspikes = NaN;
        nonburstspikes = NaN;
    end
else
    bn = NaN;
    burstiness = NaN;
    mburstf = NaN;
    mibsn = NaN;
    mbl = NaN;
    mibfr = NaN;
    mibif = NaN;
    burstf = NaN;
    ibsn = NaN;
    burstlen = NaN;
    ibfr = NaN;
    ibif = NaN;
    bp = NaN;
    normbnbdiff = NaN;
    burst = NaN;
    burstspikes = NaN;
    nonburstspikes = NaN;
end