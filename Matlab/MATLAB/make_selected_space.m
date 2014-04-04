function [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected )
%% create selected space
selected_points = find(mask_idx> 0 & mask_idx<10);
healthy_points = find(mask_idx ==10);
iindex = downsample(selected_points, ceil(length(selected_points)/Nb_Selected));%sort(uint32(rand(Nb_Selected,1)*length(selected_points)));
h_index = downsample(healthy_points, ceil(length(healthy_points)/Nb_healthy));
selected_space = [space(:,iindex),space(:,h_index)];
mask_idx = mask_idx([iindex,h_index]);


