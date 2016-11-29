function [B,twom] = modularity_f(A,gamma)
%MODULARITY_F returns monolayer Newman-Girvan modularity matrix for undirected network given by adjacency matrix A, function handle version
%
% Version: 2.1
% Date: Tue 29 Nov 2016 15:29:57 EST
% 
% Only works for undirected networks
%
%   Input: A:  NxN adjacency matrices of an undirected network
%          gamma: resolution parameter
%
%   Output: B: function handle where B(i) returns the ith column of 
%          [N]x[N] modularity matrix of the monolayer network 
%           with adjacency matrix A
%           twom: normalisation constant
%
%   Example of usage: [B,twom]=modularity_f(A,gamma);
%          [S,Q]= genlouvain(B); 
%          Q=Q/twom;
%   Notes:
%
%     The matrix A is assumed to be symmetrix and square. These assumptions
%     are not checked here.
%
%     For smaller systems, it is potentially more efficient (and easier) to
%     directly use the sparse quality/modularity matrix B in MODULARITY. For
%     large systems with directed networks, use MODULARITYDIR_F.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different null models).  
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%    References:
%     Newman, Mark E. J. and Michelle Girvan. "Finding and Evaluating 
%     Community Structure in Networks", Physical Review E 69, 026113 (2004). 
%
%   Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2016).

if nargin<2||isempty(gamma)
	gamma=1;
end

k=sum(A,2);
twom=sum(k);
    
B=@(i) full(A(:,i)-gamma*(k*k(i))/twom);

end
