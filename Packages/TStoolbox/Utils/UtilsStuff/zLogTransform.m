function z = zLogTransform(obj)

%  Y = zLogTransform(X)
%  z-transform of the log of the column vector X. If X is a matrix, z-transform of each vector is computed

nbVect = size(obj,2);
z = zeros(size(obj,1),nbVect);

for i=1:nbVect
	
	vect = obj(:,i);
	stdev = std(vect);
	
	if (stdev==0)
		zTransform = 0;
	else
		zTransform(:,i) = (vect-mean(vect))/stdev;
	end;
end;