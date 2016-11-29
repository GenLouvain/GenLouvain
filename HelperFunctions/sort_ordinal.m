function [S_sorted,s]=sort_ordinal(S)
% SORT_ORDINAL reorders nodes to emphasize persistent structure in an ordered multilayer partition
%
% Version: 2.1
% Date: Tue 29 Nov 2016 15:29:58 EST
% 
% Nodes are reordered using the optimal leave order for the 
% average linkage hierarchical clustering tree based on Hamming distance
% between community assignments
%
% Call as:
%
%     [S_sorted, s] = sort_ordinal(S)
%
% Input: 
% 
%     S: multilayer partition (matrix of size NxT, where N is the number
%        of nodes and T is the number of layers)
%
% Output:
%
%     S_sorted: reordered multilayer partition
%
%     s: mapping of nodes to reordered nodes
%
% Note that S_sorted=S(s, :)
%
% Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2016).
   

d1=pdist(S,'hamming');
Z1=linkage(d1,'average');
s=optimalleaforder(Z1,d1);
S_sorted=S(s,:);
