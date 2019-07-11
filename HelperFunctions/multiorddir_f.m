function [B,twom]=multiorddir_f(A,gamma,omega)
%MULTIORDDIR_F  returns multilayer Leicht-Newman modularity matrix for ordered directed layers, function handle version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
%   Input: A: Cell array of NxN adjacency matrices for each layer of an
%          ordered directed multilayer network
%          gamma: intralayer resolution parameter
%          omega: interlayer coupling strength
%
%   Output: B: function handle where B(i) returns the ith column of
%          [NxT]x[NxT] flattened modularity tensor for the
%           multilayer network with uniform ordinal coupling (T is
%           the number of layers of the network)
%           twom: normalisation constant
%
%   Example of usage: [B,twom]=multiorddir_f(A,gamma,omega);
%          [S,Q]= genlouvain(B); % see iterated_genlouvain.m and
%          postprocess_ordinal_multilayer.m for how to improve output
%          multilayer partition
%          Q=Q/twom;
%          S=reshape(S,N,T);
%
%   [B,twom] = MULTIORDDIR_F(A,GAMMA, OMEGA) with A a cell array of square
%   matrices of equal size each representing an directed network "layer"
%   computes the Leicht-Newman multilayer modularity matrix Susing the
%   quality function described in Mucha et al. 2010, with intralayer
%   resolution parameter GAMMA, and with interlayer coupling OMEGA
%   connecting nearest-neighbor ordered layers. Once the mulilayer modularity
%   matrix is computed, optimization can be performed by the generalized
%   Louvain code GENLOUVAIN or ITERATED_GENLOUVAIN. The output B can be used
%   with other heuristics, provided the same mapping is used to go from the
%   multilayer tensor to the multilayer flattened matrix. That is, the
%   node-layer tuple (i,s) is mapped to i + (s-1)*N. [Note that we can
%   define a mapping between a multilayer partition S_m stored as an N by T
%   matrix and the corresponding flattened partition S stored as an NT by 1
%   vector. In particular S_m = reshape(S,N,T) and S = S_m(:).]
%
%   See also
%       genlouvain heuristics:      GENLOUVAIN, ITERATED_GENLOUVAIN
%       multilayer wrappers:        MULTICAT, MULTICATF, MULTIORD
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%
%   Notes:
%     The matrices in the cell array A are assumed to be square,
%     and of equal size.  These assumptions are not checked here.
%
%     For smaller systems, it is potentially more efficient (and easier) to
%     directly use the sparse quality/modularity matrix B in MULTIORD. For
%     large systems with undirected layer networks, use MULTIORD_F.
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
%   References:
%     Blondel, Vincent D., Jean-Loup Guillaume, Renaud Lambiotte, and
%     Etienne Lefebvre, "Fast unfolding of communities in large networks,"
%     Journal of Statistical Mechanics: Theory and Experiment, P10008
%     (2008).
%
%     Fortunato, Santo, "Community detection in graphs," Physics Reports
%     486, 75-174 (2010).
%
%     Good, Benjamin H., Yves-Alexandre de Montjoye, and Aaron Clauset,
%     "Performance of modularity maximization in practical contexts,"
%     Physical Review E 81, 046106 (2010).
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Elizabeth A. Leicht and Mark E. J. Newman. "Community structure in
%     Directed Networks", Physical Review Letters 100, 118703 (2008).
%
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     Thank you to Dani Bassett, Jesse Blocher, Bruce Rogers, and Simi Wang
%     for their collaborative help which led to significant cleaning up
%     of earlier versions of our multilayer community detection codes.



if nargin<2||isempty(gamma)
    gamma=1;
end

if nargin<3||isempty(omega)
    omega=1;
end

N=length(A{1});
T=length(A);
if length(gamma)==1
    gamma=repmat(gamma,T,1);
end
m=zeros(T,1);
for i=1:T
    m(i)=sum(A{i}(:));
end
A=blkdiag(A{:});
kout=sum(A,1);
koutmat=sparse(1:(N*T),kron(1:T,ones(1,N)),kout);
kin=sum(A,2);
kinmat=sparse(1:(N*T),kron(1:T,ones(1,N)),kin);
A=(A+A')./2;
A=A+omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);

B=@(i) A(:,i)-gamma(ceil(i./(N+eps))).*(kout(i).*kinmat(:,ceil(i./(N+eps)))+kin(i).*koutmat(:,ceil(i./(N+eps))))./(2*m(ceil(i./(N+eps))));

twom=sum(m)+omega*2*N*(T-1);
end
