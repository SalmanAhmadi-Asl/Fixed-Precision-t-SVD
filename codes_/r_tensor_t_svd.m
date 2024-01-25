
function [U,S,V] = r_tensor_t_svd(A,R)
  [m,n,k]=size(A);
  C=t_prod(A,randn(n,R,k));
  [Q,~]=tqr(C);
  [U_1,S,V]=tsvds(t_prod(t_trans(Q),A),R);
  U=t_prod(Q,U_1);
end