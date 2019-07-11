function [B,twom]=bipartite_f(A,gamma)
% BIPARTITE_F returns monolayer Barber modularity matrix for undirected bipartite networks, function handle version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% Input: A: MxN adjacency matrix for an undirected bipartite network
%        gamma: resolution parameter
%
% Output: B: function handle such that B(i) returns the ith column of the
%            modularity matrix using the bipartite null-model
%            of "Barber, M. Modularity and community detection in bipartite networks.
%            Phys. Rev. E 76, 066102 (2007)".
%         twom: normalisation constant
%
% Usage: [B,twom]=bipartite(A,gamma);
%        [S,Q]=genlouvain(B);
%        Q=Q/twom;
%
%   Notes:
%     For smaller systems, it is potentially more efficient (and easier) to
%     directly use the sparse quality/modularity matrix B in BIPARTITE.
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

twom=2*mm;

    function modi=modf(i)

        if i<=m
            indx=(m+1:N);
            v=A(i,:)-gamma*k(i)*d/mm;
            modi=sparse(indx,1,v,N,1);
        else
            indx=(1:m);
            v=A(:,i-m)-gamma*k*d(i-m)/mm;
            modi=sparse(indx,1,v,N,1);
        end

    end

if mm==0
    B=@(i) sparse(N,1);
else
    B=@modf;
end

end
