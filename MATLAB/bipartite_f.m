function [B,twom]=bipartite_f(A,gamma)
%BIPARTITE_f [B,twom]=BIPARTITE_f(A,gamma)
%
% Input: A: MxN adjacency matrix an undirected bipartite network
% gamma: resolution parameter


[m,n]=size(A);
N=m+n;

k=sum(A,2);
d=sum(A,1);

mm=sum(k);

twom=2*mm;

    function modi=modf(i)
        
        
        if mm~=0
            if i<=m
                indx=(m+1:N);
                v=A(i,:)-k(i)*d/mm;
                
                modi=sparse(indx,1,v,N,1);
            else
                indx=(1:m);
                v=A(:,i-m)-gamma*k*d(i-m)/mm;
                modi=sparse(indx,1,v,N,1);
            end
            
        else
            modi=sparse(N,1);
        end
    end

B=@modf;

end