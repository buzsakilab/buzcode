function [biasmean1, biasstd1, bias_rp1] = CrossCorrBootstrap2(Sdata, S, E)
nbins = 10;
binsize = 10;

n_cells = length(Sdata);


bias = ones(n_cells,n_cells, 50) * NaN;
bias_rp = ones(n_cells,n_cells) * NaN;


for l = 1:50
  [SM, EM] = NonRippleTimes(Data(S), Data(E));
  for i = 1:n_cells
    Sdatam{i} = Restrict(Sdata{i}, SM, EM);
  end
  
  for i = 1:n_cells
    for j = 1:i-1
      [C{i,j}, B] = CrossCorr(Data(Sdatam{i}), Data(Sdatam{j}), binsize, nbins);

      if(sum(C{i,j}) == 0)
	bi = 0;
      else
	bi = (sum(C{i,j}(1:floor(nbins/2))) - ...
	    sum(C{i,j}(floor(nbins/2)+2:end)))...
	    / sum(C{i,j});
      end

      if isnan(bi) 
	error('Numerical Problem...');
      end
      bias(i,j,l) = bi;

    end
  end
end
biasmean = mean(bias,3);
biasstd = std(bias, 0, 3);
biasmean1 = biasmean(find(~isnan(biasmean)));
biasstd1 = biasstd(find(~isnan(biasstd)));


for i = 1:n_cells
  Sdatarp{i} = Restrict(Sdata{i}, Data(S), Data(E));
end


for i = 1:n_cells
  for j = 1:i-1

    [Crp{i,j}, B] = CrossCorr(Data(Sdatarp{i}), Data(Sdatarp{j}),...
	binsize, nbins);
    if(sum(Crp{i,j}) == 0)
      birp = 0;
    else
      
      birp = (sum(Crp{i,j}(1:floor(nbins/2))) -...
	  sum(Crp{i,j}(floor(nbins/2)+2:end))) / ...
	  sum(Crp{i,j});
    end
    if  isnan(birp)
      error('Numerical Problem...');
    end
    bias_rp(i,j) = birp;    

  end
end
bias_rp1 = bias_rp(find(~isnan(bias_rp)));
    
    
    
    
    


