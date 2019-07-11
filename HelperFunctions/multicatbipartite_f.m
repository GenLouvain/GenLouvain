function [B,twom] = multicatbipartite_f(A,gamma,omega)
% MULTICATBIPARTITE_F  returns multilayer Barber modularity matrix for unordered undirected bipartite networks, function handle version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% MULTICATBIPARTITE_F [B,twomu] = MULTICATBIPARTITE_F(A,gamma,omega)
%
% Input: A: Cell array of MxN adjacency matrices for each layer of a
%           multilayer undirected bipartite network
%        gamma: resolution parameter
%        omega: interlayer coupling strength
%
% Output: B: function handle where B(i) returns the ith column of the
%            [(M+N)xT]x[(M+N)xT] flattened modularity tensor for the
%            multilayer bipartite network with uniform ordinal coupling (T
%            is the number of layers of the network)
%         twomu: normalisation constant
%
% Usage: [B,twomu]=multicatbipartite_f(A,gamma,omega);
%        [S,Q]=genlouvain(B); % see iterated_genlouvain.m and
%          postprocess_categorical_multilayer.m for how to improve output
%          multilayer partition
%        Q=Q/twom;
%        S=reshape(S,M+N,T);
%
%  [B,twom] = MULTICATBIPARTITE_F(A,GAMMA, OMEGA) with A a cell array of
%   matrices of equal size each representing an undirected bipartite network
%   "layer" computes the multilayer Barber modularity matrix using the quality
%   function described in Mucha et al. 2010, with intralayer resolution
%   parameter GAMMA, and with interlayer coupling OMEGA connecting all-to-all
%   layers. Once the mulilayer modularity matrix is computed,
%   optimization can be performed by the generalized Louvain code GENLOUVAIN
%   or ITERATED_GENLOUVAIN. The  output B can be used with other heuristics,
%   provided the same mapping is used to go from the multilayer tensor to
%   the multilayer flattened matrix. That is, the node-layer tuple (i,s)
%   is mapped to i + (s-1)*(M+N). [Note that we can define a mapping between
%   a multilayer partition S_m stored as an (M+N) by T matrix and the
%   corresponding flattened partition S stored as an MNT by 1 vector. In
%   particular S_m = reshape(S,M+N,T) and S = S_m(:). Note that nodes i=1:M
%   correspond to the first class (i.e. the rows of A) and nodes i=M+1:M+N
%   correspond to the second class (i.e. the columns of A) of the bipartite
%   network.]
%
%   Notes:
%     The matrices in the cell array A are assumed to be of equal size.
%     This assumption is not checked here.
%
%     For smaller systems, it is potentially more efficient (and easier) to
%     directly use the sparse quality/modularity matrix B in MULTICATBIPARTITE.
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

k=zeros(m,T);
d=zeros(T,n);
mm=zeros(T,1);

twom=0;
for j=1:T
    twom = twom + sum(sum(A{j}));
    k(:,j)=sum(A{j},2);
    d(j,:)=sum(A{j});
    mm(j)=sum(k(:,j));
end

%interslice connections
all2all= N*[(-T+1):-1,1:(T-1)];
C=omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);



%bipartite modularity matrix
    function modi=modf(i)

        s=ceil(i/(N+eps));
        if mm(s)~=0
            ii=i-(s-1)*N;
            if ii<=m
                indx=(m+1:N)+(s-1)*N;
                v=A{s}(ii,:)-gamma(s)*k(ii,s)*d(s,:)/mm(s);

                modi=sparse(indx,1,v,N*T,1,n+2);
            else
                indx=(1:m)+(s-1)*N;
                v=A{s}(:,ii-m)-gamma(s)*k(:,s)*d(s,ii-m)/mm(s);

                modi=sparse(indx,1,v,N*T,1,m+2);
            end

            modi=modi+C(:,i);
        else
            modi=C(:,i);

        end
    end

B=@modf;
twom=2*twom+2*N*(T-1)*T*omega;
end
