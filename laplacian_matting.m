function [t_est] = laplacian_matting( t_est, I, lambda)
%LAPLACIAN_MATTING implementation of Natural Image Matting as proposed by
%Levin et al

[m,n] = size(t_est);
L = getLaplacian1(I,zeros(m,n));
[Lm Ln] = size(L);
Atmp = (L+lambda*speye(Lm,Ln));
clear L
btmp = (lambda.*reshape(t_est,m*n,1));
t_est = double(Atmp)\double(btmp);
clear Atmp

t_est = reshape(t_est,m,n);


end

