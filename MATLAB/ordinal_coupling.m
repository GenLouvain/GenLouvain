function B=ordinal_coupling(B,omega)
% convert modularity matrices for each layer (given by cellarray `B` into
% multilayer modularity matrix with ordinal coupling of strength
% `omega`
T=length(B);
N=length(B{1});

B=spblkdiag(B{:})+omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);

end
