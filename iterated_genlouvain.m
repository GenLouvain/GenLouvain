function [S,Q,n_it]=iterated_genlouvain(B,limit,verbose,randord,randmove,S0,postprocessor)
% Optimise modularity-like quality function by iterating GenLouvain until convergence.
% (i.e., until output partition does not change between two successive iterations)
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:43 CEST
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B) with matrix B repeatedly implements
%   GenLouvain until convergence to an output partition. After the first
%   iteration, the starting partition of each subsequent iteration is the
%   output partition of the previous iteration. The output vector S encodes
%   the community assignments of the output partition after the final
%   iteration, with S(i) identifying the community to which node i has been
%   assigned. The output Q gives the quality of the partition S of the
%   network. The output n_it gives the observed number of iterations until
%   convergence.
%   NOTE: The matrix represented by B must be both symmetric and square.
%   This condition is not checked thoroughly if B is a function handle, but
%   is essential to the proper use of this routine. When B is a matrix,
%   non-symmetric input is symmetrised (B=(B+B')/2), which preserves the
%   quality function.
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B) with function handle B such that
%   B(i) returns the ith column of the modularity/quality matrix uses this
%   function handle (to reduce the memory footprint for large networks)
%   until the number of groups is less than 10000 and then builds the B
%   matrix corresponding to the new aggregated network in subsequent passes.
%   Use [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit) to change this
%   default=10000 limit.
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,0) suppresses displayed text
%   output.
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,verbose,0) forces index-ordered
%   (cf. randperm-ordered) consideration of nodes, for deterministic results
%   with randord = 'move'.
%
%   [S,Q,n_it]=ITERATED_GENLOUVAIN(B,limit,verbose,randord,randmove) controls
%   additional randomization to obtain a broader sample of the quality
%   function landscape. The possible values for 'randmove' are
%       'move': always move node under consideration to the community that
%           results in maximal improvement in modularity (default)
%       'moverand': move the node under consideration to a community chosen
%           uniformly at random from all moves that increase the qualilty
%           function
%       'moverandw': move the node under consideration to a community chosen
%           at random from all moves that increase the qualilty where the
%           probability of choosing a particular move is proportional to
%           its increase in the quality function
%       0: equivalent to 'move' (provided for backwards compatibility)
%       1: equivalent to 'moverand' (provided for backwards compatibility)
%
%   'moverand', and 'moverandw' mitigate some undesirable behavior for
%   "multilayer" modularity with ordinal coupling ('moverandw' tends to be
%   better behaved for large values of the interlayer coupling). With
%   'move', the algorithm exhibits an abrupt change in behavior when the
%   strength of the interlayer coupling approaches the maximum value of the
%   intralayer modularity matrices (see Bazzi et al. 2016 for more detail).
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,verbose,randord,randmove,S0)
%   uses S0 as an inital partition for the first iteration (the starting
%   partition of subsequent iterations is the output partition of the previous
%   iteration). The default choice for S0 is all singletons (and given by a
%   length(B) by 1 vector). If a user specifies S0, it needs to satisfy
%   numel(S0) = length(B). Note that the size of S will match the size of S0.
%   In a multilayer setting, a user may have to reshape S appropriately
%   (e.g., reshape(S,N,T), where N is the number of nodes in each layer and
%   T is the number of layers).
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,verbose,randord,randmove,S0,
%   postprocessor) applies the function handle postprocessor to the
%   output partition of GENLOUVAIN at each iteration. The function handle
%   postprocessor takes a partition as input and returns a partition as
%   output. To ensure convergence, the function postprocessor should
%   increase the value of the (multilayer) quality function for a partition
%   (e.g, kernigan-lin type algorithms in a monolayer or multilayer setting,
%   and "postprocess-ordinal-multilayer.m",
%   "postprocess-categorical-multilayer.m" in HelperFunctions for a multilayer
%   setting)
%
%   Example on multilayer network quality function of Mucha et al. 2010
%   (using multilayer cell A with A{s} the adjacency matrix of layer s)
%
%   gamma = 1; omega = 0.1;
%   N=length(A{1});
%   T=length(A);
%
%   B = multiord(A,gamma,omega); % multiord.m computes multilayer
%   modularity matrix B with homogeneous ordinal interlayer coupling w and
%   a Newman-Girvan null model on each layer (see more detail in
%   documentation of multiord.m)
%
%   PP = @(S) postprocess_ordinal_multilayer(S,T); % define postprocessing
%   function handle that increases multilayer modularity without changing
%   intralayer partitions in an ordered multilayer networks (see
%   postprocess_ordinal_multilayer.m for more detail)
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,10000,0,1,'moverandw',[], PP);
%   S = reshape(S, N, T);
%
%   Notes:
%
%     Under default options, this routine can return different results from
%     run to run because it considers nodes in pseudorandom (randperm)
%     order.  Because of the potentially large number of nearly-optimal
%     partitions (Good et al. 2010), one is encouraged to investigate
%     results of repeated applications of this code (and, if possible, of
%     other computational heuristics).  To force deterministic behavior with
%     randord = 'move', ordering nodes by their index, pass zero as the
%     fourth input: ITERATED_GENLOUVAIN(B,limit,verbose,0).
%
%     This algorithm is only "Louvain-like" in the sense that the two
%     phases are used iteratively in the same manner as in the Louvain
%     algorithm (Blondel et al. 2008).  Because it operates on a general
%     quality/modularity matrix B, it does not include any analytical
%     formulas for quickly identifying the change in modularity from a
%     proposed move nor any improved efficiency obtained by their use.  If
%     your problem uses one of the well-used null models included in other
%     codes, those codes should be much faster for your task.
%
%     Past versions had a problem where accumulated subtraction error might
%     lead to an infinite loop with each pass oscillating between two or
%     more partitions yet incorrectly identifying increases in quality.  We
%     believe this problem has been corrected by the relative change checks
%     in lines 178 and 269.  If you encounter a similar problem, notify
%     Peter Mucha (<a href="mailto:mucha@unc.edu">mucha@unc.edu</a>).
%
%     The output Q provides the sum over the appropriate elements of B
%     without any rescaling.
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
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
%     A special thank you to Stephen Reid, whose greedy.m code was the
%     original version that has over time developed into the present code.
%
%     Thank you also to Dani Bassett, Jesse Blocher, Mason Porter and Simi
%     Wang for inspiring improvements to the code.
%
%   Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and
%       Peter J. Mucha (2011-2019).
%
%   See also genlouvain HelperFunctions


%set default for maximum size of modularity matrix
if nargin<2
    limit = [];
end

%set level of reported/displayed text output
if nargin<3||isempty(verbose)
    verbose = false;
end


%set randperm- v. index-ordered
if nargin<4
    randord = [];
end

%set move function (maximal (original Louvain) or random improvement)
if nargin<5
    randmove=[];
end

% set initial partition
if nargin<6
    S0=[];
end

% set postprocessing function
if nargin<7
    postprocessor=[];
end

% verbose output switch
if verbose
    mydisp = @(s) disp(s);
else
    mydisp = @(s) [];
end

S_old=[];
n_it=1;
mydisp('Iteration 1');
[S,Q]=genlouvain(B,limit,verbose,randord,randmove,S0);

mydisp('');

Q_old=-inf;
while ~isequal(S,S_old)&&(Q-Q_old)>10*eps
    n_it=n_it+1;
    S_old=S;
    Q_old=Q;

    mydisp(sprintf('Iteration %u',n_it));
    if ~isempty(postprocessor)
        S=postprocessor(S);
    end
    [S,Q]=genlouvain(B,limit,verbose,randord,randmove,S);
    mydisp(sprintf('Improvement in modularity: %f\n',Q-Q_old));
end

end
