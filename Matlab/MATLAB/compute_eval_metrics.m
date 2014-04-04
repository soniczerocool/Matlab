function [dice,jaccard,precision,recall] = compute_eval_metrics(auto_lbl,lbl)

%auto_lbl = segmentation,lbl = ground truth

ind1 = find(auto_lbl>0);
ind2 = find(lbl>0);

if(length(ind1) ~= 0 & length(ind2) ~= 0)
  dice      = (100*((2*length(intersect(ind1,ind2)))/(length(ind1)+length(ind2))));
  jaccard   = (100*(length(intersect(ind1,ind2))/length(union(ind1,ind2))));
  precision = (100*(length(intersect(ind1,ind2))/length(ind1)));
  recall    = (100*(length(intersect(ind1,ind2))/length(ind2)));
else
  dice      = 0;
  jaccard   = 0;
  precision = 0;
  recall    = 0;
end
