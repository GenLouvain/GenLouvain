function [B,twom] = multiordbipartite_f(A,omega,gamma)
%MULTIORD  Multislice bipartite community detection for ordered slices, function
%version
%   Version 1.0, August 26, 2011.
%
%   [S,Q] = multiordbipartite_f(A,OMEGA) with A a cell array of bipartite
%   matrices of equal size each representing an undirected network "slice"
%   performs multislice community detection using the quality function
%   described in Mucha et al. 2010, with interslice coupling OMEGA
%   connecting nearest-neighbour ordered slices. 
%   
%   Outputs a function handle for the B, so that B(i) is the ith column of
%   the modularity matrix, that can be used to obtain communities using the
%   genlouvainrand algorithm.
%
%   See also
%       multislice wrappers:        MULTICATF, MULTIORD, MULTIORDF
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%
%   Notes:
%     The matrices in the cell array A are assumed to be of equal size.
%     These assumptions are not checked here.
%
%  
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different intraslice null models,
%     different numbers of nodes from slice-to-slice, or systems which are
%     both multiplex and longitudinal).  That is, this code is only a
%     starting point; it is by no means exhaustive.
%
%     twom gives the 1/(2*mu) prefactor used to normalise the output of
%     genlouvainrand.
%
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
%       Inderjit S. Jutla, Lucas Jeub, and Peter J. Mucha, "A generalized Louvain method
%       for community detection implemented in MATLAB,"
%       http://netwiki.amath.unc.edu/GenLouvain (2011).

if nargin<3
	gamma=1;
end


[m,n]=size(A{1});
N=m+n;
T=length(A);
k=zeros(m,T);
d=zeros(n,T);
mm=zeros(T,1);

twom=0;
for j=1:T
    twom = twom + sum(sum(A{j}));
    k(:,j)=sum(A{j},2);
    d(:,j)=sum(A{j});
    mm(j)=sum(k(:,j));
end

k=sparse(k);
d=sparse(d);
mm=full(mm);

%bipartite modularity matrix        
B= @(i) modf_bipord(i,m,n,T,A,k,d,mm,omega,gamma);

%normalising prefactor                
twom=2*twom+2*N*T*omega;
end