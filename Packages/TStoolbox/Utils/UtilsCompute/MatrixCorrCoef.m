function r = MatrixCorrCoef(A, B)
% r = MatrixCorrCoef(A, B) 
%
% Computes the correlation coefficient between the vectors of the
% off-diagonal elements of the square matrices A and B
% INPUTS
% A, B: two square matrices of the same size
% OUTPUTS
% r: the correlation coefficient 

% batta 2000

[N, N1] = size(A);
if N ~= N1
  error('A must be square');
end

[NB, N1] = size(B);
if NB ~= N1
  error('B must be square');
end

if N ~= NB
  error('A and B must be of the same size');
end

cA = [];
cB = [];



for k = 1:N-1
  cA = [cA;diag(A,k)];
  cB = [cB;diag(B,k)];
end

R = corrcoef(cA, cB);

r = R(1,2);
