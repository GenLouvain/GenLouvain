function [B,twom]=stability(A,t)

k=sum(A,2);
twom=sum(k);
L=A./repmat(k,1,length(A))-diag(ones(length(A),1));
L=expm(t*L);
if isequal(A,A')
    B=L.*repmat(k./twom,1,length(A))-(k*k')/twom^2;
else
    'error'
end
end