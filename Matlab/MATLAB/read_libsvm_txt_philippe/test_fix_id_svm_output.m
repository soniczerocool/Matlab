classes = unique (MASK(:));
classes([1,end],:)= classes([end,1],:);
classes(end) = [];

for c=1:length(classes)
    idx_mask_c = find(MASK == classes(c));
    ix = randperm(length(idx_mask_c));
    random_point_c = idx_mask_c(ix(1));
    [dummy,c_id] = max(posterior_matrix(:,random_point_c));
    posterior_matrix(:,[c,c_id]) = posterior_matrix(:,[c_id,c]);
end


