function [S_res] = reshape_layers(S,T)
%[S_RES]=RESHAPE_LAYERS(S,T) reshapes each columns of S to have T layers. 

k=size(S,2);
n=size(S,1);
S_res=cell(k,1);

for i=1:k
    S_res{i}=reshape(S(:,i),n/T,T);
end

% if k==1
%     S_res=S_res{1};
% end

end

