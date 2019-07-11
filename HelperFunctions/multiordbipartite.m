function [B,twomu] = multiordbipartite(A,gamma,omega)
% MULTIORDBIPARTITE  returns multilayer Barber modularity matrix for ordered undirected bipartite networks, matrix version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% MULTIORDBIPARTITE [B,twomu] = MULTIORDBIPARTITE(A,gamma,omega)
%
% Input: A: Cell array of MxN adjacency matrices for each layer of a
%           multilayer undirected bipartite network
%        gamma: resolution parameter
%        omega: interlayer coupling strength
%
% Output: B: [(M+N)xT]x[(M+N)xT] flattened modularity tensor for the
%            multilayer bipartite network with uniform ordinal coupling (T is
%            the number of layers of the network)
%         twomu: normalisation constant
%
% Usage: [B,twomu]=multiordbipartite(A,gamma,omega);
%        [S,Q]=genlouvain(B); % see iterated_genlouvain.m and
%          postprocess_ordinal_multilayer.m for how to improve output
%          multilayer partition
%        Q=Q/twom;
%        S=reshape(S,M+N,T);
%
%  [B,twom] = MULTIORDBIPARTITE(A,GAMMA, OMEGA) with A a cell array of
%   matrices of equal size each representing an undirected bipartite network
%  "layer" computes the Barber multilayer modularity matrix using the quality
%   function described in Mucha et al. 2010, with intralayer resolution
%   parameter GAMMA, and with interlayer coupling OMEGA connecting
%   nearest-neighbor ordered layers.  Once the mulilayer modularity matrix
%   is computed, optimization can be performed by the generalized Louvain
%   code GENLOUVAIN or ITERATED_GENLOUVAIN. The sparse output matrix B can
%   be used with other heuristics, provided the same mapping is used to go
%   from the multilayer tensor to the multilayer flattened matrix. That is,
%   the node-layer tuple (i,s) is mapped to i + (s-1)*(M+N). [Note that we
%   can define a mapping between a multilayer partition S_m stored as an
%   (M+N) by T matrix and the corresponding flattened partition S stored as
%   an MNT by 1 vector. In particular S_m = reshape(S,M+N,T) and S = S_m(:).
%   Note that nodes i=1:M correspond to the first class (i.e. the rows of A)
%   and nodes i=M+1:M+N correspond to the second class (i.e. the columns of
%   A) of the bipartite network.]
%
%
%
%   Notes:
%     The matrices in the cell array A are assumed to be of equal size.
%     This assumption is not checked here.
%
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try MULTIORDBIPARTITE_F.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different intralayer null models,
%     different numbers of nodes from layer-to-layer, or systems which are
%     both multiplex and longitudinal).  That is, this code is only a
%     starting point; it is by no means exhaustive.
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
% References:
%       Barber, M. Modularity and community detection in bipartite networks.
%           Phys. Rev. E 76, 066102 (2007).
%
%       Mucha, P. J., Richardson, T., Macon, K., Porter, M. A. & Onnela, J.-P.
%           Community structure in time-dependent, multiscale, and multiplex networks.
%           Science 328, 876-878 (2010).

if nargin<2||isempty(gamma)
    gamma=1;
end

if nargin<3
    omega=1;
end

[m,n]=size(A{1});
N=m+n;
T=length(A);
if length(gamma)==1
    gamma=repmat(gamma,T,1);
end

B=spalloc(N*T,N*T,2*m*n*T+2*N*T);
mu=0;
for s=1:T
    k=sum(A{s},2);
    d=sum(A{s},1);
    mm=sum(k);

    if mm==0
        B1=sparse(m,n);
    else
        B1=A{s}-gamma(s)*k*d/mm;
    end

    indx=(s-1)*N;
    B(indx+1:indx+m,indx+m+1:indx+N)=B1;
    B(indx+m+1:indx+N,indx+1:indx+m)=B1';

    mu=mm+mu;

end

B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
twomu=2*mu+2*N*(T-1)*omega;
end
