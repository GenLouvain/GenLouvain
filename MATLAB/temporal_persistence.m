function pers=temporal_persistence(S)

if ~iscell(S)
    S={S};
end

[N,T]=size(S{1});
pers=zeros(length(S),1);
for i=1:length(S)
    pers(i)=sum(sum(S{i}(:,1:end-1)==S{i}(:,2:end)))/(N*(T-1));
end

end
