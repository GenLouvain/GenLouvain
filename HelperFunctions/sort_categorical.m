function [S_sorted,s1,s2]=sort_categorical(S)
% SORT_CATEGORICAL reorders nodes and layers to emphasize persistent structure in an unordered multilayer partition
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% Nodes and layers are reordered using the optimal leave order for the
% average linkage hierarchical clustering tree based on Hamming distance
% between community assignments
%
% Call as:
%
%     [S_sorted, s1, s2] = sort_categorical(S)
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
%     s1: mapping of nodes to reordered nodes
%
%     s2: mapping of layer to reordered layers
%
% Note that S_sorted=S(s1,s2)


d1=pdist(S,'hamming');
d2=pdist(S','hamming');

Z1=linkage(d1,'average');
Z2=linkage(d2,'average');

s1=optimalleaforder(Z1,d1);
s2=optimalleaforder(Z2,d2);
S_sorted=S(s1,s2);
