function S=postprocess_temporal_multilayer(S,T,max_coms)

% post-process an NxT multilayer partition S with N nodes and T layers to
% maximise persistence without changing the community structure within
% layers. Uses the Hungarian algorithm to solve the optimal assignment
% problem for each consecutive pair of layers. Note that this procedure
% always increases temporal multilayer modularity for any non-zereo value
% of omega.
if nargin<2||isempty(T)
    T=size(S,2);
end
N=numel(S)/T;


if nargin<3||isempty(max_coms)
    max_coms=inf;
end

[N0,T0]=size(S);

if max(S(:))<max_coms % don't do anything if too many communities for performance
    S=reshape(S,N,T);
    max_com=max(S(:,1));
    for i=2:T
        [u1,~,e1]=unique(S(:,i-1)); % unique communities in previous layer
        [u2,~,e2]=unique(S(:,i)); % unique communities in this layer
        G1=sparse(e1,1:length(e1),1); % community assignment matrix for previous layer
        G2=sparse(1:length(e2),e2,true); % community assignment matrix for this layer
        overlap=G1*G2; % node overlap matrix between communities in the two layers
        dist=sum(overlap(:))-overlap;
        S2=assignmentoptimal(dist'); % find best assignment for communities in current layer
        for j=1:length(u2)
            if S2(j)~=0
                S(G2(:,j),i)=u1(S2(j)); % update assignment
            else
                S(G2(:,j),i)=max_com+1; % get new label if not assigned
                max_com=max_com+1;
            end
        end
    end
    
end

% return in original format
S=tidy_config(S);
S=reshape(S,N0,T0);
end


