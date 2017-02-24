function [clu mergehistory] = mergeclu_slow(clu,res,fet,tR,tC,rogueSpk2Merge)


%Parameters
% Maximal fraction of rogue spike in a merged cluster
% rogueSpk2Merge = 0.1;
showtree = 0;

mergehistory=[];
sametree=1;
step=1;

% Worst programming ever
while step
    
    %OK, here, we want to recompute te errormatrix as rarely as possible.
    %We'll do it once first, and then we'll see if after merging, the new
    %errormatrix shows some others cluster to merge or not. 

    step=0;    
    
    % Compute the errormatrix
    em = errormatrix(fet,clu);
    
    % Make it symetric
    ems = max(em,em');
    
    % Some tricks to transform it as a vector for linkage analysis
    y = squareform((1-ems)-diag(diag(1-ems)),'tovector');
    Z = linkage(y);
    
    % Let's plot the binary tree
    if showtree
        figure(1),clf
        [H,T,perm] = dendrogram(Z,0);
    else
        [T,perm] = dendrogram_struct(Z,0);
    end
    T0 = T;
    iniClu = size(Z,1)+1;
    cluIx = unique(clu);
    
    % There are N-2 possible merged clusters. In other word, there are 2N-2
    % leaves in the tree.    
    cluIx = [cluIx;zeros(length(cluIx)-2,1)];

    % Once again, worst programming ever
    while sametree
        % let's do whatever we can do with this errormatrix. If no clusters
        % should be merged, sametree is still at 0 and the progam goes one
        % loop above to recompute a new errormatrix
        sametree=0;
        tomerge = [];
        for ii=1:size(Z,1)
            x=Z(ii,1);
            y=Z(ii,2);
            if any(T==x) && any(T==y) 
                rgx = double(res(clu==cluIx(x)))./20000;
                rgy = double(res(clu==cluIx(y)))./20000;
                %Compute the fraction of rogue spike fxy (fcorr is not
                %used)
                [fxy fcorr] = crossrefract(rgx,rgy,tR,tC);
                rgxRefract = FractionRogueSpk(rgx,tR,tC);
                rgyRefract = FractionRogueSpk(rgy,tR,tC);
                %Common inter-spike intervals
                Icsi = icsi(rgx,rgy);
                
                %compute some statistics
                [h px py] =  IcsiStat(rgx,rgy,Icsi,tR,tC,.05);
                temp = corrcoef(mean(fet(clu==cluIx(x),:)),mean(fet(clu==cluIx(y),:)));
                cor = temp(2);

                if fxy<rogueSpk2Merge ... % is the common refractory bad? (consider removing this)
                    && h ...   % are the ISI distributinos similar?
                    && fxy < 1.2*rgxRefract && fxy < 1.2*rgyRefract ...  % ensure we don't make a refractory period worse..
                    && cor > .95  % correlation between feature spaces should be quite high if its a single unit..
                    
                    if showtree
                        set(H(ii),'Color','r')
                    end
                    tomerge = [tomerge;ii];
                    %At least two leaves should be merged!
                    sametree=1;
                    step=1;
                else
                    T(T==x)=[];
                    T(T==y)=[];
                end
            end
        end
        % Go across all the pairs that should be merged.
        for ii=1:size(tomerge,1)
            x = Z(tomerge(ii),1);
            y = Z(tomerge(ii),2);
            T(T==x)=[];
            T(T==y)=[];
            T = [T;tomerge(ii)+iniClu];
            clu = updateclu(clu,[cluIx(x);cluIx(y)],max(cluIx)+1);
            cluIx(tomerge(ii)+iniClu) = max(cluIx)+1;
            mergehistory=[mergehistory;[cluIx(x) cluIx(y) max(cluIx)]];
        end
       
    end
    if ~step
        cluIx = unique(clu);
        for ii=1:length(cluIx)
            clu = updateclu(clu,cluIx(perm(ii)),max(cluIx)+ii+1);
        end
    end
end




