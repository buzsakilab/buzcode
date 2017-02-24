function z = zTransform(obj)

%  Y = zLogTransform(X)
%  z-transform of the log of the column vector X. If X is a matrix, z-transform of each vector is computed

nbVect = size(obj,2);
z = zeros(size(obj,1),nbVect);

for i=1:nbVect
	
	
	vect = obj(:,i);
	stdev = std(vect);
	
	if (stdev==0)
		z(:,i) = 0;
	else
		z(:,i) = (vect-mean(vect))/stdev;
	end;
end;