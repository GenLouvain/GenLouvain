function [ Spar, Qpar ] = parameter_iteration( wrapper,A,param )
%PARAMETER_ITERATION [ Spar, Qpar ] = parameter_iteration( wrapper,A,param )
%
%Finds modularity maximum with GenLouvainRand for wrapper for each value of
%in vector param

iter=length(param);
% if iscell(A)
%     T=length(A);
% else
%     T=1;
% end
% 
T=1;

parfor i=1:iter
    
    [B,twom]=feval(wrapper,A,param(i));
    [Spar{i},Qpar(i)]=genlouvain(B,[],[],1,1);
    Qpar(i)=Qpar(i)/twom;
    n=size(Spar{i},1);
    Spar{i}=reshape(Spar{i},n/T,T);
    
end

if ~iscell(A)
    for i=1:iter
        Spara(:,i)=Spar{i};
    end
    Spar=Spara;
end