function [ B,twom ] = mod_matrix_f( A ,gamma)
%MOD_MATRIX returns modularity matrix for network given by Adjacency matrix A.
%
%   supports giving A as a function handle, returning ith column of Adjacency matrix
%	and resolution parameter GAMMA (defaults to 1).

if nargin<2
	gamma=1;
end

    
    
    
   
        k=sum(A,2);

    
    m2=sum(k);
    
    B=@(i) A(:,i)-gamma*(k*k(i))/m2;
    twom=m2;



end
