function [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(auto_lbl,lbl, NumberOfLabels)

% recall is a vector which elemets are recall for edema, non_enhanced ,
% enhanced, and the complete(abnormality vs healthy) recall. 
% the same indexing is adopted for precision, jaccard and dice


if nargin < 3
    NumberOfLabels = 4;
end



recall = zeros(1,NumberOfLabels);
precision = zeros(1,NumberOfLabels);
jaccard = zeros(1,NumberOfLabels);
dice = zeros(1,NumberOfLabels);
%% statistical measures for edema (label 2)
ind1 = find(auto_lbl==2);
ind2 = find(lbl==2);

dice(1)      = (100*((2*length(intersect(ind1,ind2)))/(length(ind1)+length(ind2))));
jaccard(1)   = (100*(length(intersect(ind1,ind2))/length(union(ind1,ind2))));
precision(1) = (100*(length(intersect(ind1,ind2))/length(ind1)));
recall(1)    = (100*(length(intersect(ind1,ind2))/length(ind2)));

%% statistical measures for non-enhanced (labels 1 and 3)
% get the indexes of voxels with label 1 (nicroses)
ind1_a = find(auto_lbl==1);
ind2_a = find(lbl==1);
% get the indexes of voxels with label 1 (non_enhanced tumore)
ind1_b = find(auto_lbl==3);
ind2_b = find(lbl==3);

% unify them as one class of non_enhanced
ind1 = union(ind1_a,ind1_b);
ind2 = union(ind2_a,ind2_b);

dice(2)      = (100*((2*length(intersect(ind1,ind2)))/(length(ind1)+length(ind2))));
jaccard(2)   = (100*(length(intersect(ind1,ind2))/length(union(ind1,ind2))));
precision(2) = (100*(length(intersect(ind1,ind2))/length(ind1)));
recall(2)    = (100*(length(intersect(ind1,ind2))/length(ind2)));

%% statistical measures for enhanced tumor (label 4)
ind1 = find(auto_lbl==4);
ind2 = find(lbl==4);

dice(3)      = (100*((2*length(intersect(ind1,ind2)))/(length(ind1)+length(ind2))));
jaccard(3)   = (100*(length(intersect(ind1,ind2))/length(union(ind1,ind2))));
precision(3) = (100*(length(intersect(ind1,ind2))/length(ind1)));
recall(3)    = (100*(length(intersect(ind1,ind2))/length(ind2)));

%% statistical measures for complete (labels 1,2,3,4)
ind1 = find(auto_lbl>0);
ind2 = find(lbl>0);

dice(4)      = (100*((2*length(intersect(ind1,ind2)))/(length(ind1)+length(ind2))));
jaccard(4)   = (100*(length(intersect(ind1,ind2))/length(union(ind1,ind2))));
precision(4) = (100*(length(intersect(ind1,ind2))/length(ind1)));
recall(4)    = (100*(length(intersect(ind1,ind2))/length(ind2)));


