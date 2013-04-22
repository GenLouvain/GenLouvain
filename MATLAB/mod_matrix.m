function [ B,twom ] = mod_matrix( A ,gamma)
%MOD_MATRIX returns modularity matrix for network given by Adjacency matrix A.
%
%   supports giving A as a function handle, returning ith column of Adjacency matrix
%	and resolution parameter GAMMA (defaults to 1).

if nargin<2
	gamma=1;
end

if isa(A,'function_handle')
    n=length(A(1));
    
    k=zeros(n,1);
    parfor i=1:n
        k(i)=sum(feval(A,i));
    end
    
    m2=sum(k);
    
    B=@(i) A(i)-gamma*(k*k(i))/m2;
    twom=m2;

else
k=sum(A);
m2=sum(k);


B=(A-gamma*(k'*k)/m2);
B=full(B);

twom=m2;
end
end

