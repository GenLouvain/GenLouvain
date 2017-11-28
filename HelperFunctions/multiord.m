function [B,twom] = multiord(A,gamma,omega)
%MULTIORD  returns multilayer Newman-Girvan modularity matrix for ordered layers, matrix version
% Works for directed or undirected networks
%
% Version: 2.1.2
% Date: Tue Nov 28 14:20:20 EST 2017
% 
%   Input: A: Cell array of NxN adjacency matrices for each layer of an
%          ordered multilayer (directed or undirected) network
%          gamma: intralayer resolution parameter
%          omega: interlayer coupling strength
%
%   Output: B: [NxT]x[NxT] flattened modularity tensor for the
%           multilayer network with uniform ordinal coupling (T is
%           the number of layers of the network)
%           mm: normalisation constant
%
%   Example of usage: [B,mm]=multiord(A,gamma,omega);
%          [S,Q]= genlouvain(B); % see iterated_genlouvain.m and 
%          postprocess_temporal_multilayer.m for how to improve output
%          multilayer partition
%          Q=Q/mm;
%          S=reshape(S,N,T);
%
%   [B,mm] = MULTIORD(A,GAMMA, OMEGA) with A a cell array of square
%   (symmetric or assymetric) matrices of equal size each representing a 
%   directed or undirected network "layer" computes the Newman Girvan multilayer 
%   modularity matrix using the quality function described in Mucha et al. 
%   2010, with intralayer resolution parameter GAMMA, and with interlayer 
%   coupling OMEGA connecting nearest-neighbor ordered layers.  The null 
%   model used for the quality function is the Newman-Girvan null model 
%   (see e.g. Bazzi et al. for other possible null models). Once the 
%   mulilayer modularity matrix is computed, optimization can be performed
%   by the generalized Louvain code GENLOUVAIN or ITERATED_GENLOUVAIN. The 
%   sparse output matrix B can be used with other heuristics, provided the 
%   same mapping is used to go from the multilayer tensor to the multilayer 
%   flattened matrix. That is, the node-layer tuple (i,s) is mapped to 
%   i + (s-1)*N. [Note that we can define a mapping between a multilayer 
%   partition S_m stored as an N by T matrix and the corresponding flattened 
%   partition S stored as an NT by 1 vector. In particular S_m = reshape(S,N,T) 
%   and S = S_m(:).] 
%
%   See also
%       genlouvain heuristics:      GENLOUVAIN, ITERATED_GENLOUVAIN
%       multilayer wrappers:        MULTICAT, MULTICATF, MULTIORDF
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%
%   Notes:
%     The matrices in the cell array A are assumed to be square,
%     and of equal size.  These assumptions are not checked here.
%
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try MULTIORD_F for undirected layer networks and MULTIORDDIR_F
%     for directed layer networks. 
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different intralayer null models [see eg 
%     Bazzi et al. 2016 for examples], different numbers of nodes from 
%     layer-to-layer, or systems which are both multiplex and longitudinal). 
%     That is, this code is only a starting point; it is by no means 
%     exhaustive.
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
%     Elizabeth A. Leicht and Mark E. J. Newman. "Community structure in
%     Directed Networks", Physical Review Letters 100, 118703 (2008). 
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Bazzi, Marya, Mason A. Porter, Stacy Williams, Mark McDonald, Daniel
%     J. Fenn, and Sam D. Howison. "Community Detection in Temporal 
%     Multilayer Networks, with an Application to Correlation Networks", 
%     MMS: A SIAM Interdisciplinary Journal 14, 1-41 (2016). 
%
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     Thank you to Dani Bassett, Jesse Blocher, Bruce Rogers, and Simi Wang
%     for their collaborative help which led to significant cleaning up
%     of earlier versions of our multilayer community detection codes.
%
%   Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2016).

if nargin<2
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
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
twom=twom+2*N*(T-1)*omega;

