function pers=categorical_persistence(S)
% CATEGORICAL_PERSISTENCE computes the persistence of an unordered multilayer partition
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% pers = categorical_PERSISTENCE(S) with a single multilayer partition or a
% cell of multilayer partitions S (with S{i} the ith multilayer partition of
% cell S) computes the "persistence" of multilayer partitions. The output pers
% is a vector such that pers(i) is the persistence of multilayer partition
% S{i} (stored as an N*T matrix where N is the number of nodes in each layer
% and T is the number of layers). For a multilayer network with unordered
% layers, the value of persistence is the sum over all pairs of layers of
% the number of nodes that do not change community assignments. Categorical
% persistence varies between 0 and 1. Categorical persistence is related to
% the multilayer quality function developed in Mucha et al. 2010 when one
% uses categorical interlayer coupling.
%
%   References:
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Bazzi, Marya, Mason A. Porter, Stacy Williams, Mark McDonald, Daniel
%     J. Fenn, and Sam D. Howison. "Community Detection in Temporal
%     Multilayer Networks, with an Application to Correlation Networks",
%     MMS: A SIAM Interdisciplinary Journal 14, 1-41 (2016). 

if ~iscell(S)
    S={S};
end

[N,T]=size(S{1});
all2all = N*[(-T+1):-1,1:(T-1)];
A=spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
pers=zeros(length(S),1);
for i=1:length(S)

    G=sparse(1:length(S{i}(:)),S{i}(:),1);
    pers(i)=trace(G'*A*G)/(N*T*(T-1));
end

end
