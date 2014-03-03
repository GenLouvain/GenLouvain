function [B,twom] = bipartite(A,gamma)
%BIPARTITE [B,twom]=BIPARTITE(A,gamma)
%
% Input: A: MxN adjacency matrix for an undirected bipartite network
%        gamma: resolution parameter

[M,N]=size(A);

k=sum(A,2);
d=sum(A);
m=sum(k);

B1=A-gamma*k*d/m;

B=sparse(M+N,M+N);
B(1:M,M+1:M+N)=B1;
B(M+1:M+N,1:M)=B1';

twom=2*m;

end

