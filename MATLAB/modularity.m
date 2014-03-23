function [B,twom] = modularity(A,gamma)
%MODULARITY returns modularity matrix for network given by Adjacency matrix A.


if nargin<2
	gamma=1;
end

k=full(sum(A));
twom=sum(k);

B=full(A-gamma*(k'*k)/twom);

end

