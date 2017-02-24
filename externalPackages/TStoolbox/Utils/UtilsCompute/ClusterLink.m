function KL = ClusterLink(L)



ncl = size(L,1);


KL = cell(1, ncl);

done = 0;

for i = 1:ncl
  KL{i} = [i];
end



while ~done 
  done = 1;
  
  for i = 1:ncl
    for j = (i+1):ncl
      if(sum(sum(L(KL{i},KL{j}))) > 0 )
	KL{i} = unique([KL{i}, KL{j}]);
	KL{j} = [];
	done = 0;
      end
    end
  end
  

  ncl = length(KL);
  nl = zeros(1,ncl);
  for i = 1:ncl
    nl(i) = length(KL{i});
  end
  KL = KL(find(nl));
  ncl = length(KL);
end

      
  