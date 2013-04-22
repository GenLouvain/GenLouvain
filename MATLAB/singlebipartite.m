function [B,twom] = singlebipartite( A,gamma)
%SINGLEBIPARTITE [B,twom]=SINGLEBIPARTITE(A,gamma)
%
% Input: A: MxN adjacency matrix an undirected bipartite network
%        gamma: resolution parameter

[M,N]=size(A);

k=sum(A,2);
d=sum(A);
m=sum(k);
if sum(d)~=m
   warning('singlebipartite:edgeweights','edgeweights don''t agree, difference: %f',full(abs(m-sum(d))));
end


B1=A-gamma*k*d/m;

B=sparse(M+N,M+N);
B(1:M,M+1:M+N)=B1;
B(M+1:M+N,1:M)=B1';

twom=2*m;

end

