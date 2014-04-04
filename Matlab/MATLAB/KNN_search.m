function [IDX,D] = KNN_search(selected_space, space, k)

% space = space(:,:);   
% space_repmat = repmat(space,[1 1 size(selected_space,2)]);
% selected_space_3d = reshape(selected_space,size(selected_space,1),1,size(selected_space,2));
% selected_space_repmat = repmat(selected_space_3d,[1 size(space,2) 1]);
IDX = zeros (size(space,2),k);
D = zeros (size(space,2),k);
for i = 1:size(space,2)
    
point = space(:,i);
tmp = sqrt((selected_space-repmat(point,1,size(selected_space,2))).^2);
tmp = sum(tmp,1);
[tmp IX] = sort(tmp);
D(i,:) = tmp(1:k);
IDX (i,:)= IX(1:k);



end
