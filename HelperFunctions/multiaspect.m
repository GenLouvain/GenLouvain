function [B, twom] = multiaspect(A, gamma, omega, type)
% MULTIASPECT  returns multilayer Newman-Girvan modularity matrix for multiple aspects.
%
% Aspects can be ordered or unordered, and the function supports different
% intralayer resolution parameters for each layer and different interlayer
% coupling parameters for each aspect.
%
% Version: 2.1.2
% Date: Tue Nov 28 14:20:20 EST 2017
%
%   Input: A: Cell array of NxN adjacency matrices for each layer of a
%             multilayer network. Each dimension of A corresponds to an aspect
%
%          gamma: intralayer resolution parameter (scalar or array with the
%                 same size as A)
%
%          omega: interlayer coupling strength (scalar or vector of length
%                 ndims(A))
%
%          type: string specifying the coupling type for each aspect.
%                Categorical (mulitplex) coupling is specified using the
%                letters 'c' or 'm' and ordinal (temporal) coupling is
%                specified using the letters 'o' or 't'.
%
%
%   Output: B: [N x numel(A)]x[N x numel(A)] flattened modularity
%              tensor for the multilayer network (note that numel(A) is the
%              number of layers of the network)
%
%           twom: normalisation constant
%
%   Example of usage:
%          [B,twom]=multicat(A,gamma,omega,type);
%          [S,Q]= genlouvain(B); % see iterated_genlouvain.m for how to
%          improve output multilayer partition
%          Q=Q/twom;
%          S=reshape(S,[N,size(A)]);
%
%          For a multiplex temporal network, one should have `ndims(A)==2`
%          and `type='co'` (or `type='mt'`)
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
%
%   Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2016).

    aspects = size(A);
    if numel(omega) == 1
        omega = repmat(omega, numel(aspects), 1);
    end
    if numel(gamma) == 1
        gamma = repmat(gamma, aspects);
    end
    na = prod(aspects);
    [B, twom] = arrayfun(@(A, g) modularity(A{1}, g), A, gamma, ...
                         'UniformOutput', false);
    B = blkdiag(B{:});
    C = zeros(prod(aspects));
    for a = 1:numel(aspects)
        C = C + kron(kron(eye(prod(aspects(a+1:end))),coupling(aspects(a), type(a))*omega(a)),...
                     eye(prod(aspects(1:a-1))));
    end
    C = kron(C, eye(length(A{1})));
    B = B + C;
    twom = sum([twom{:}]) + sum(sum(C));
end

function C = coupling(size, type)
    switch type
        case {'m', 'c'}
            C = ones(size)-eye(size);
        case {'t', 'o'}
            C = diag(ones(size-1,1),1);
            C = C + C';
        otherwise
            error('unknown aspect type %s', type)
    end
end
