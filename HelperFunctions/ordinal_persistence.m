function pers=ordinal_persistence(S)
% ORDINAL_PERSISTENCE computes the persistence of an ordered multilayer partition
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% pers = ORDINAL_PERSISTENCE(S) with a single multilayer partition or a
% cell of multilayer partitions S (with S{i} the ith multilayer partition
% of cell S) computes the "persistence" of ordered multilayer partitions. The
% output pers is a vector such that pers(i) is the persistence of multilayer
% partition S{i} (stored as an N*T matrix where N is the number of nodes in
% each layer and T is the number of layers). For a multilayer partition with
% ordered layers, the value of persistence is the number of nodes that do
% not change community assignments between consecutive layers (see Bazzi
% et al. 2016 for more detail). Ordinal persistence varies between 0 and 1.
% Ordinal persistence is related to the multilayer quality function
% developed in Mucha et al. 2010 when one uses ordinal interlayer coupling.
%
%   References:
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Bazzi, Marya, Mason A. Porter, Stacy Williams, Mark McDonald, Daniel
%     J. Fenn, and Sam D. Howison. "Community Detection in Temporal Multilayer
%     Networks, with an Application to Correlation Networks", MMS: A SIAM
%     Interdisciplinary Journal 14, 1-41 (2016).


if ~iscell(S)
    S={S};
end

[N,T]=size(S{1});
pers=zeros(length(S),1);
for i=1:length(S)
    pers(i)=sum(sum(S{i}(:,1:end-1)==S{i}(:,2:end)))/(N*(T-1));
end

end
