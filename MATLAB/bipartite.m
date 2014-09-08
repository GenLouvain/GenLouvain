function [B,twom] = bipartite(A,gamma)
%BIPARTITE [B,twom]=BIPARTITE(A,gamma)
%
% Input: A: MxN adjacency matrix for an undirected bipartite network
%        gamma: resolution parameter
%
% Output: B: (M+N)x(M+N) modularity matrix using the bipartite null-model
%            of "Barber, M. Modularity and community detection in bipartite networks. 
%            Phys. Rev. E 76, 066102 (2007)".
%         twom: normalisation constant
% 
% Usage: [B,twom]=bipartite(A,gamma);
%        [S,Q]=genlouvain(B);
%        Q=Q/twom;
%
% Lucas Jeub
% jeub@maths.ox.ac.uk

[m,n]=size(A);
N=m+n;

k=sum(A,2);
d=sum(A,1);
mm=sum(k);

B1=A-gamma*k*d/mm;

B=sparse(N,N);
B(1:m,m+1:N)=B1;
B(m+1:N,1:m)=B1';

twom=2*mm;

end

