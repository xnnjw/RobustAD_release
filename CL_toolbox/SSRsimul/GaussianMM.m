function A = GaussianMM(m,n)
% Generates Gaussian measurement matrix of size m x n 

A = (1/2)*(randn(m,n)+1i*randn(m,n));
len = diag(sqrt(A'*A)); 
A = A*diag(1./len);      %  normalized to  unit-norm columns 