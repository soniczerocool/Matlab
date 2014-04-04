for beta = 0.0001:0.0005:0.001
  disp(['Processing', ' with beta = ', num2str(beta)])   
    segmented =  markov_network_multilable5(T1C,energy,beta);

[h,w,d] = size(T1C);
segmented_volume = label_to_matrix(segmented,space,h,w,d);

time = toc(t);
%% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;

%[dice,jaccard,precision,recall] = compute_eval_metrics(segmented,truth);


%% save results

mha_result_path = ['findbeta\markov_segment_',num2str(beta),'.mha'];
segmented_brain=uint16(segmented_brain);
 writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);

% 

end