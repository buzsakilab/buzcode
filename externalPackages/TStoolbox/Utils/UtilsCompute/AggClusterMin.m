function KL = AggClusterMin(D)



ncl = size(L,1);


KL = cell(1, ncl);

done = 0;

for i = 1:ncl
  KL{i} = [i];
end

iter = 1;

while ncl > 1
  done = 1;
  
  DK = (NaN) * ones(ncl, ncl);
  [J, I] =meshgrid(1:ncl, 1:ncl);
  I = reshape(I, 1, ncl*ncl);
  J = reshape(J, 1, ncl*ncl);
  
  
  for i = 1:ncl
    for j = (i+1):ncl
      D_mn = D(KL{i}, KL{j});
      DK(i,j) = nanmin(nanmin(D_mn))
    end
  end
  
      
  [min_D, ix] = min(reshape(DK, 1, ncl*ncl));
  
  i = I(ix);
  j = J(ix);
 
  KL{i} = unique([KL{i}, KL{j}]);
  KL{j} = [];
  done = 0;
 

  ncl = length(KL);
  nl = zeros(1,ncl);
  for i = 1:ncl
    nl(i) = length(KL{i});
  end
  KL = KL(find(nl));
  ncl = length(KL);
  KL_save{iter} = KL;
  
end

      
  