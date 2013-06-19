function [B,mm] = multiord(A,gamma,omega)
%MULTIORD  Multislice community detection for ordered slices, matrix version
%   Version 0.99, August 26, 2011.
%
%   [S,Q] = MULTIORD(A,OMEGA) with A a cell array of square symmetric
%   matrices of equal size each representing an undirected network "slice"
%   performs multislice community detection using the quality function
%   described in Mucha et al. 2010, with interslice coupling OMEGA
%   connecting nearest-neighbor ordered slices.  Optimization is performed
%   by the generalized Louvain code GENLOUVAINRAND by Jutla & Mucha
%   (following from Blondel et al. 2008).  The output matrix S encodes the
%   obtained community assignments, with S(i,s) identifying the community
%   to which node i in slice s has been assigned.  The output Q gives the
%   quality of the resulting partition of the network.  The sparse output
%   matrix B is the equivalent quality/modularity matrix that can be used
%   with other heuristics, ordering multislice nodes looping over identify
%   followed by slices [corresponding to S(:) reordering].
%
%   See also
%       multislice wrappers:        MULTICAT, MULTICATF, MULTIORDF
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%
%   Notes:
%     The matrices in the cell array A are assumed to be symmetric, square,
%     and of equal size.  These assumptions are not checked here.
%
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try MULTIORDF.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different intraslice null models,
%     different numbers of nodes from slice-to-slice, or systems which are
%     both multiplex and longitudinal).  That is, this code is only a
%     starting point; it is by no means exhaustive.
%
%     This version can return different results from run to run because it
%     considers nodes in pseudorandom (randperm) order.  Because of the
%     potentially large number of nearly-optimal partitions (Good et al.
%     2010), one is encouraged to investigate results of repeated
%     applications of this and other codes.
%
%     The quality Q output here does not include the 1/(2*mu) prefactor.
%
%     Resolution parameters within each slice can be added in line 90.
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
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     Thank you to Dani Bassett, Jesse Blocher, Bruce Rogers, and Simi Wang
%     for their collaborative help which led to significant cleaning up
%     of earlier versions of our multislice community detection codes.
%
%   Citation: If you use this code, please cite as
%       Inderjit S. Jutla and Peter J. Mucha, "A generalized Louvain method
%       for community detection implemented in MATLAB,"
%       http://netwiki.amath.unc.edu/GenLouvain (2011).

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
mm=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
	mm=mm+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma(s).*(k'*k/twom);
end
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
mm=mm+2*N*T*omega;

