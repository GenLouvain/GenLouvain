function [B,twom] = modularity(A,gamma)
%MODULARITY returns modularity matrix for network given by Adjacency matrix A.
% Works for directed and undirected networks


if nargin<2
	gamma=1;
end

k=sum(A,2);
d=sum(A,1);
twom=sum(k);

B=full((A+A')/2-gamma/2*(k*d+d'*k)/twom);

end

