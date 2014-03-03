function [B,twom] = multiordbipartite_f_nomex(A,gamma,omega)
%MULTIORDBIPARTITE [B,TWOM]=MULTIORDBIPARTITE(A,gamma,omega)
%
% Input: A: cell array of mxn bipartite adjacency matrices for each slice
% of the network
%        gamma: resolution parameter
%        omega: interslice connection strength

[m,n]=size(A{1});
N=m+n;
T=length(A);
k=zeros(m,T);
d=zeros(T,n);
mm=zeros(T,1);

twom=0;
for j=1:T
    twom = twom + sum(sum(A{j}));
    k(:,j)=sum(A{j},2);
    d(j,:)=sum(A{j});
    mm(j)=sum(k(:,j));
end

%interslice connections
C=omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);



%bipartite modularity matrix
    function modi=modf(i)      
        
        s=ceil(i/N);
        if mm(s)~=0        
            ii=i-(s-1)*N;
            if ii<=m
                indx=(m+1:N)+(s-1)*N;            
                v=A{s}(ii,:)-gamma*k(ii,s)*d(s,:)/mm(s);
 
                modi=sparse(indx,1,v,N*T,1,n+2);
            else
                indx=(1:m)+(s-1)*N;
                v=A{s}(:,ii-m)-gamma*k(:,s)*d(s,ii-m)/mm(s);

                modi=sparse(indx,1,v,N*T,1,m+2);
            end
            
            modi=modi+C(:,i);
        else
            modi=C(:,i);
          
        end
    end

B=@modf;
twom=2*twom+2*N*T*omega;
end