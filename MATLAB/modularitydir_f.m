function [B,twom] = modularitydir_f(A,gamma)
%MODULARITY_F returns modularity matrix as a function handle for network given by Adjacency matrix A.
%works for directed networks

if nargin<2
	gamma=1;
end

k=sum(A,2);
d=sum(A,1);
twom=sum(k);
A=(A+A')/2;
    
B=@(i) full(A(:,i)-gamma/2*(k*d(i)+d'*k(i))/twom);

end
