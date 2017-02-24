function M = convnNormBorder(A,G)

%  Normalised convn(.,.,'same) function. This function avoids border effect too. 
%  Adrien Peyrache 2007
%  
%  M = convn(A,G,'same');
%  n = convn(ones(size(A)),G,'same');
%  M = M./n;

a = size(A)
M = zeros(3*a);
M(1:a(1),1:a(2)) = A(end:-1:1,end:-1:1);
M(1:a(1),a(2)+1:2*a(2)) = A(end:-1:1,:);
M(1:a(1),2*a(2)+1:3*a(2)) = A';
M(a(1)+1:2*a(1),2*a(2)+1:3*a(2)) = A(:,end:-1:1);
M(2*a(1)+1:3*a(1),2*a(2)+1:3*a(2)) = A(end:-1:1,end:-1:1);
M(2*a(1)+1:3*a(1),a(2)+1:2*a(2)) = A(end:-1:1,:);
M(2*a(1)+1:3*a(1),1:a(2)) = A(end:-1:1,end:-1:1);
M(a(1)+1:2*a(1),1:a(2)) = A(:,end:-1:1);

M(a(1)+1:2*a(1),a(2)+1:2*a(2)) = A;

M = convn(M,G,'same');
M = M(a(1)+1:2*a(1),a(2)+1:2*a(2));

