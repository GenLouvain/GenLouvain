function [B,twom] = modularity_f(A,gamma)
%MODULARITY_F returns modularity matrix as a function handle for network given by Adjacency matrix A.


if nargin<2
	gamma=1;
end

k=sum(A,2);
twom=sum(k);
    
B=@(i) full(A(:,i)-gamma*(k*k(i))/twom);

end
