function [B,twom]=bipartite_f(A,gamma)
%BIPARTITE_f [B,twom]=BIPARTITE_f(A,gamma)
%
% Input: A: MxN adjacency matrix for an undirected bipartite network
% gamma: resolution parameter

[m,n]=size(A);
N=m+n;

k=sum(A,2);
d=sum(A,1);

mm=sum(k);

twom=2*mm;

    function modi=modf(i)
        
        if i<=m
            indx=(m+1:N);
            v=A(i,:)-gamma*k(i)*d/mm;
            modi=sparse(indx,1,v,N,1);
        else
            indx=(1:m);
            v=A(:,i-m)-gamma*k*d(i-m)/mm;
            modi=sparse(indx,1,v,N,1);
        end
        
    end

if mm==0
    B=@(i) sparse(N,1);
else
    B=@modf;
end

end