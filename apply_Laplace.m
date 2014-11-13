function [u] = apply_Laplace(f, lambda1, lambda2, m, n)
% u = (lambda1*Laplace1 + lambda2*Laplace2)f
% in square domain of size (m, n)
% where Laplacek is the Laplace operator in k-th dimension
c1 = [-1 zeros(1, n-2) 1];
c1T = [-1 1 zeros(1, n-2)];
c2 = [-1 zeros(1, m-2) 1]';
c2T = [-1 1 zeros(1, m-2)]';
C1 = zeros(m, n);
C2 = zeros(m, n);
for ii=1:m
	C1(ii,:) = fft(c1T) .* fft(c1);
end
for jj=1:n
	C2(:,jj) = fft(c2T) .* fft(c2);
end
F = reshape(f, m, n);
U = ifft2(fft2(F).* (lambda1 * C1 + lambda2 * C2));
u = reshape(U, m*n, 1);
end