function p = fPolyFit(x, y, n)
x = x(:);
V = ones(length(x), n + 1);   % Vandermonde matrix
for j = n:-1:1
   V(:, j) = V(:, j + 1) .* x;
end
[Q, R] = qr(V, 0);            % Solve least squares problem
p      = transpose(R \ (transpose(Q) * y(:)));