function [bias1] = CrossCorrBias2(Sdata, nbins, binsize)



n_cells = length(Sdata);



bias = ones(n_cells,n_cells) * NaN;






for i = 1:n_cells
  for j = 1:i-1

    [Crp{i,j}, B] = CrossCorr(Data(Sdata{i}), Data(Sdata{j}),...
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
    bias(i,j) = birp;    

  end
end
bias1 = bias(find(~isnan(bias)));
    
    
    
    
    


