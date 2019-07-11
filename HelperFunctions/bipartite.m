function [B,twom] = bipartite(A,gamma)
% BIPARTITE returns monolayer Barber modularity matrix for undirected bipartite networks, matrix version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
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
%   Notes:
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try BIPARTITE_F for undirected bipartite networks.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different null models).
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
% References:
%       Barber, M. Modularity and community detection in bipartite networks.
%           Phys. Rev. E 76, 066102 (2007).

if nargin<2||isempty(gamma)
    gamma=1;
end

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
