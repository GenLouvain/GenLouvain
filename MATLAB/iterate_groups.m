function [ Sdist,Qdist,B ] = iterate_groups(iter, wrapper,A,varargin )
%	ITERATE_GROUPS [Sdist,Qdist,B]=ITERATE_GROUPS(iter, wrapper, A, varargin)
%
%	applies wrapper #iter times to A, where wrapper is a function handle to one of the wrapper functions
%	returning the modularity matrix   



if isempty(varargin) 
[B,twom]=wrapper(A);
else
[B,twom]=wrapper(A,varargin{:});
end

parfor i=1:iter
        [Sprem,Qprem]=genlouvain(B,7000,1,1,1);
   
   Sdist(:,i)=Sprem;
   Qdist(i)=Qprem/twom;
end


%	reshaping for time-dependent optimisation (uncomment if desired)
% if iscell(A)
%     T=length(A);
%     Sdistcell=cell(T,1);
%     tot=length(Sdist(:,1));
%     for i=1:iter
%         Sdistcell{i}=reshape(Sdist(:,i),tot/T,T);
%     end
%     Sdist=Sdistcell;
% end

end

