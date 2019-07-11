function [B,twom] = multicat(A,gamma,omega)
%MULTICAT  returns multilayer Newman-Girvan modularity matrix for unordered layers, matrix version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
%   Input: A: Cell array of NxN adjacency matrices for each layer of an
%          unordered multilayer undirected network
%          gamma: intralayer resolution parameter
%          omega: interlayer coupling strength
%
%   Output: B: [NxT]x[NxT] flattened modularity tensor for the
%           multilayer network with uniform categorical coupling (T is
%           the number of layers of the network)
%           twom: normalisation constant
%
%   Example of usage: [B,twom]=multicat(A,gamma,omega);
%          [S,Q]= genlouvain(B); % see iterated_genlouvain.m and
%          postprocess_categorical_multilayer.m for how to improve output
%          multilayer partition
%          Q=Q/twom;
%          S=reshape(S,N,T);
%
%   [B,twom] = MULTICAT(A,GAMMA, OMEGA) with A a cell array of square
%   symmetric matrices of equal size each representing an undirected network
%   "layer" computes the multilayer modularity matrix using the quality
%   function described in Mucha et al. 2010, with intralayer resolution
%   parameter GAMMA, and with interlayer coupling OMEGA connecting
%   all-to-all categorical layers. Once the mulilayer modularity matrix is
%   computed, optimization can be performed by the generalized Louvain code
%   GENLOUVAIN or ITERATED_GENLOUVAIN. The sparse output matrix B can be used
%   with other heuristics, provided the same mapping is used to go from the
%   multilayer tensor to the multilayer flattened matrix. That is, the
%   node-layer tuple (i,s) is mapped to i + (s-1)*N. [Note that we can define
%   a mapping between a multilayer partition S_m stored as an N by T matrix
%   and the corresponding flattened partition S stored as an NT by 1 vector.
%   In particular S_m = reshape(S,N,T) and S = S_m(:).]
%
%   See also
%       genlouvain heuristics:      GENLOUVAIN, ITERATED_GENLOUVAIN
%       multilayer wrappers:        MULTICATF, MULTIORD, MULTIORDF
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%
%   Notes:
%     The matrices in the cell array A are assumed to be square,
%     symmetric, and of equal size.  These assumptions are not checked here.
%
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try MULTICAT_F for undirected layer networks.
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
%     Newman, Mark E. J. and Michelle Girvan. "Finding and Evaluating
%     Community Structure in Networks", Physical Review E 69, 026113 (2004).
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
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

if nargin<3
omega=1;
end

N=length(A{1});
T=length(A);

if length(gamma)==1
gamma=repmat(gamma,T,1);
end

B=spalloc(N*T,N*T,N*N*T+2*N*T);
twom=0;
for s=1:T
    kout=sum(A{s},1);
    kin=sum(A{s},2);
    mm=sum(kout);
	twom=twom+mm;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=(A{s}+A{s}')/2-gamma(s)/2.*((kin*kout+kout'*kin')/mm);
end
all2all = N*[(-T+1):-1,1:(T-1)];
B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
twom=twom + (N*T*(T-1)*omega);
