mainfolder='J:\User\Research\Brain_data\BRATS-1\Images\';
results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\Plots\';
if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end

list=dir(mainfolder);
list = list(3:end);
% path_brains = cell(1,length(list));
% j=1;
% for i=1:length(list)
%     dirpath =['J:\User\Research\Brain_data\BRATS-1\Images\',list(i).name];
%     if isdir(dirpath)
%      path_brains{j}= dirpath;
%      j=j+1;
%     end
% end

for i=1:length(list)
    dir_name =['J:\User\Research\Brain_data\BRATS-1\Images\',list(i).name];
    result_path_brain = [results_path_root, list(i).name];
    if exist(result_path_brain,'dir')==0
        mkdir(result_path_brain)
    end
    list_modality = dir (dir_name);
    list_modality = list_modality(2:end);
    for j=1:length(list_modality)
        name_modality = list_modality(j);
        if ~isempty(strfind(name_modality.name,'T1C.nii'))
            filepath = [dir_name,'\',name_modality.name];
            T1C = load_nii(filepath, [], 1);
            T1C=T1C.img; 
        end
        if ~isempty(strfind(name_modality.name,'T2.nii'))
            filepath = [dir_name,'\',name_modality.name];
            T2 = load_nii(filepath, [], 1);
            T2=T2.img; 
        end
        if ~isempty(strfind(name_modality.name,'truth.nii'))
            filepath = [dir_name,'\',name_modality.name];
            truth = load_nii(filepath, [], 1);
            truth=truth.img; 
        end
    end
    T1C = T1C(:);
    T2 = T2(:);
    truth = truth(:);
    
    T1C = double(T1C);
    T2 = double(T2);
    
    T1C = ((T1C - min(T1C)) / max(T1C)) * 2000; 
    T2 = ((T2 - min(T2)) / max(T2)) * 2000; 

    healthyidx=find(truth<.5);
    tumoridx=find(.5<truth & truth<=1.5);
    edemaidx=find(truth>1.5);
   %%
    [T1C_hist n] = hist(T1C,1000);
    T1C_hist= T1C_hist(2:end);
    fig1=figure(1);
    bar(n(2:end),T1C_hist)
    title('T1C hist')
    saveas(fig1,[result_path_brain,'\','T1C_all_hist.jpg'])
    
    fig2 = figure(2);
    [T1C_hist_healthy n] = hist(T1C(healthyidx),1000);
    T1C_hist_healthy = T1C_hist_healthy(2:end);
    plot(n(2:end),T1C_hist_healthy,'.blue')
    hold on
    [T1C_hist_edema n] = hist(T1C(edemaidx),1000);
    T1C_hist_edema = T1C_hist_edema(2:end);
    plot(n(2:end),T1C_hist_edema,'.green')
    hold on
    [T1C_hist_tumor n] = hist(T1C(tumoridx),1000);
    T1C_hist_tumor = T1C_hist_tumor(2:end);
    plot(n(2:end),T1C_hist_tumor,'.red')
    hold off
    legend('healthy','edema','tumor')
    title('T1C hist')
    saveas(fig2,[result_path_brain,'\','T1C_mixed_hist.jpg'])

 %%  
    
    
    [T2_hist n] = hist(T2,1000);
    T2_hist= T2_hist(2:end);
    fig1=figure(1);
    bar(n(2:end),T2_hist)
    title('T2 hist')
    saveas(fig1,[result_path_brain,'\','T2_all_hist.jpg'])
    
    fig2 = figure(2);
    [T2_hist_healthy n] = hist(T2(healthyidx),1000);
    T2_hist_healthy = T2_hist_healthy(2:end);
    plot(n(2:end),T2_hist_healthy,'.blue')
    hold on
    [T2_hist_edema n] = hist(T2(edemaidx),1000);
    T2_hist_edema = T2_hist_edema(2:end);
    plot(n(2:end),T2_hist_edema,'.green')
    hold on
    [T2_hist_tumor n] = hist(T2(tumoridx),1000);
    T2_hist_tumor = T2_hist_tumor(2:end);
    plot(n(2:end),T2_hist_tumor,'.red')
    hold off
    legend('healthy','edema','tumor')
    title('T2 hist')
    saveas(fig2,[result_path_brain,'\','T2_mixed_hist.jpg'])
  %%  
    T1Cc = T1C - mean(T1C);
    [T1Cc_hist n] = hist(T1Cc,1000);
    T1Cc_hist= T1Cc_hist(2:end);
    fig1=figure(1);
    bar(n(2:end),T1Cc_hist)
    title('T1Cc hist')
    saveas(fig1,[result_path_brain,'\','T1Cc_all_hist.jpg'])
    
    fig2 = figure(2);
    [T1Cc_hist_healthy n] = hist(T1Cc(healthyidx),1000);
    T1Cc_hist_healthy = T1Cc_hist_healthy(2:end);
    plot(n(2:end),T1Cc_hist_healthy,'.blue')
    hold on
    [T1Cc_hist_edema n] = hist(T1Cc(edemaidx),1000);
    T1Cc_hist_edema = T1Cc_hist_edema(2:end);
    plot(n(2:end),T1Cc_hist_edema,'.green')
    hold on
    [T1Cc_hist_tumor n] = hist(T1Cc(tumoridx),1000);
    T1Cc_hist_tumor = T1Cc_hist_tumor(2:end);
    plot(n(2:end),T1Cc_hist_tumor,'.red')
    hold off
    legend('healthy','edema','tumor')
    title('T1Cc hist')
    saveas(fig2,[result_path_brain,'\','T1Cc_mixed_hist.jpg'])
    %%
    
    T2c = T2 - mean(T2);
     [T2c_hist n] = hist(T2c,1000);
    T2c_hist= T2c_hist(2:end);
    fig1=figure(1);
    bar(n(2:end),T2c_hist)
    title('T2c hist')
    saveas(fig1,[result_path_brain,'\','T2c_all_hist.jpg'])
    
    fig2 = figure(2);
    [T2c_hist_healthy n] = hist(T2c(healthyidx),1000);
    T2c_hist_healthy = T2c_hist_healthy(2:end);
    plot(n(2:end),T2c_hist_healthy,'.blue')
    hold on
    [T2c_hist_edema n] = hist(T2c(edemaidx),1000);
    T2c_hist_edema = T2c_hist_edema(2:end);
    plot(n(2:end),T2c_hist_edema,'.green')
    hold on
    [T2c_hist_tumor n] = hist(T2c(tumoridx),1000);
    T2c_hist_tumor = T2c_hist_tumor(2:end);
    plot(n(2:end),T2c_hist_tumor,'.red')
    hold off
    legend('healthy','edema','tumor')
    title('T2c hist')
    saveas(fig2,[result_path_brain,'\','T2c_mixed_hist.jpg'])
    

    
 %% 
    
    h=figure(3);
    plot(downsample(T1Cc(healthyidx),20) ,downsample(T2c(healthyidx),20),'.b');
    hold on;

    plot(downsample(T1Cc(edemaidx),20),downsample(T2c(edemaidx),20),'.g');
    plot(downsample(T1Cc(tumoridx),20),downsample(T2c(tumoridx),20),'.r');
    axis([-1000 2000 -1000 2000])
    saveas(h,[result_path_brain,'\','T1T2_meansubtracted.jpg']);
    hold off
    h=figure(4);
    plot(downsample(T1C(healthyidx),20),downsample(T2(healthyidx),20),'.b');
    hold on
    plot(downsample(T1C(tumoridx),20),downsample(T2(tumoridx),20),'.r');
    plot(downsample(T1C(edemaidx),20),downsample(T2(edemaidx),20),'.g');
    axis([-200 2000 -200 2000])
    saveas(h,[result_path_brain,'\','T1T2.jpg']);
    hold off;

end
