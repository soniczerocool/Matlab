function make_txt_file(space,selected_space, truth_idx,mask_idx,type, name)
folder = 'N:\User\Research\Brain_data\Brats-2_challenge\';
%folder = '/Users/uoft/Desktop/text_files/';
brain_folder = [folder,type,'_',name,'\']
mkdir(brain_folder)
save_interaction = [brain_folder, 'interaction_',type,'_',name,'.txt'];
save_allpoints = [brain_folder, 'allpoints_',type,'_',name,'.txt'];
save_background = [brain_folder,'background_',type,'_',name,'.txt'];

f_interaction = fopen(save_interaction,'w');
f_allpoints   = fopen(save_allpoints,'w');
f_background  = fopen(save_background,'w');

%  selected_points = find(mask_idx>0);
% selected_space = space(:,selected_points);
% mask_idx = mask_idx(selected_points);

%% downsample
% selected_points = find(mask_idx> 0 & mask_idx<10);
% healthy_points = find(mask_idx ==10);
% iindex = downsample(selected_points, ceil(length(selected_points)/Nb_Selected));%sort(uint32(rand(Nb_Selected,1)*length(selected_points)));
% h_index = downsample(healthy_points, ceil(length(healthy_points)/Nb_healthy));
% selected_space = [space(:,iindex),space(:,h_index)];
% mask_idx = mask_idx([iindex,h_index]);
%%
h_idx = find(mask_idx==10);
mask_idx(h_idx) = 0;


% %remove background from selected space
% background = find(sum(selected_space(1:3,:))< 0.00001);
% selected_space(:,background) = [];
% mask_idx(background) = [];


for p=1:size(space,2)
    a = space(:,p);
    label = truth_idx(p);
    flag_background =0;
    if sum(space(1:3,p))< 0.001
        flag_background = 1;
    end
    fprintf(f_allpoints,int2str(label));
    for f = 1: 6
        if a(f)~=0
            fprintf(f_allpoints,[' ',int2str(f),':',num2str(a(f))]);
        end
    end
    fprintf(f_allpoints,'\n');

    fprintf(f_background,num2str(flag_background));
    fprintf(f_background,'\n');
end

for p = 1:size(selected_space,2)
    a = selected_space(:,p);
    label = mask_idx(p);
   
    
    fprintf(f_interaction,int2str(label));
    for f = 1: 6
        if a(f)~=0
            fprintf(f_interaction,[' ',int2str(f),':',num2str(a(f))]);
        end
    end
    fprintf(f_interaction,'\n');

    
end

fclose(f_interaction);
fclose(f_allpoints);
fclose(f_background);


zip_interaction = [brain_folder,'interaction_',type,'_',name];
zip_allpoints = [brain_folder,'allpoints_',type,'_',name];
zip_background = [brain_folder,'background_',type,'_',name];


zip(zip_interaction,save_interaction )
zip(zip_allpoints, save_allpoints)
zip(zip_background,save_background )
delete(save_interaction, save_allpoints,save_background)