function idx_all = get_dist_idx(P, dist, rate)

idx_all=zeros(size(P,1), 2);

dP=zeros(size(P));
dP(2:end,:)=P(2:end,:)-P(1:end-1,:);
d=sqrt(dP(:,1).^2+dP(:,2).^2+dP(:,3).^2);
d=cumsum(d);

i=1;

while 1
    dd=d-d(i);
    idx=find(dd>=dist*(1-rate) & dd<=dist*(1+rate));
    if isempty(idx)
        break;
    end
    idx_all(i,:)=[i, idx(1)];
    i = i+1;
end

idx_all=idx_all(1:i-1,:);

end

