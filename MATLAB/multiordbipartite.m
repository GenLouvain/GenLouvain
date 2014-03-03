function [B,twom] = multiordbipartite(A,gamma,omega)
%MULTIORDBIPARTITE [B,TWOM]=MULTIORDBIPARTITE(A,gamma,omega)
%
% Input: A: cell array of mxn bipartite adjacency matrices for each slice
% of the network
%        gamma: resolution parameter
%        omega: interslice connection strength

[m,n]=size(A{1});
N=m+n;
T=length(A);
B=spalloc(N*T,N*T,2*m*n*T+2*N*T);
twom=0;

for s=1:T
    k=sum(A{s},2);
    d=sum(A{s});
    mm=sum(k);
    
    if mm==0
        B1=sparse(m,n);
    else
        B1=A{s}-gamma*k*d/mm;
    end
        
    indx=(s-1)*N;    
    B(indx+1:indx+m,indx+m+1:indx+N)=B1;
    B(indx+m+1:indx+N,indx+1:indx+m)=B1';
    
    twom=mm+twom;
   
end

B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
twom=2*twom+2*N*T*omega;

end