function [B,twom] = modularitydir_f(A,gamma)
% MODULARITYDIR_F returns monolayer Leicht-Newman modularity matrix for directed network given by adjacency matrix A, function handle version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
%
%   Input: A:  NxN adjacency matrices of a directed network
%          gamma: resolution parameter
%
%   Output: B: function handle where B(i) returns the ith column of
%          [N]x[N] modularity matrix of the monolayer network
%           with adjacency matrix A
%           twom: normalisation constant
%
%   Example of usage: [B,twom]=modularitydir_f(A,gamma);
%          [S,Q]= genlouvain(B);
%          Q=Q/twom;
%   Notes:
%     The matrix A is assumed to be square. This assumption is not checked
%     here.
%
%     For smaller systems, it is potentially more efficient (and easier) to
%     directly use the sparse quality/modularity matrix B in MODULARITY. For
%     large systems with undirected networks, use MODULARITY_F.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different null models).
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%   References
%     Elizabeth A. Leicht and Mark E. J. Newman. "Community structure in
%     Directed Networks", Physical Review Letters 100, 118703 (2008).


if nargin<2||isempty(gamma)
	gamma=1;
end

k=sum(A,2);
d=sum(A,1);
twom=sum(k);
A=(A+A')/2;

B=@(i) full(A(:,i)-gamma/2*(k*d(i)+d'*k(i))/twom);

end
