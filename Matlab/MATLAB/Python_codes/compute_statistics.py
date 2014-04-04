def compute_eval_multilabel_metrics(auto_lbl, lbl, NumberOfLabels = 4):
	'''recall is a vector which elemets are recall for edema, non_enhanced ,
	 enhanced, and the complete(abnormality vs healthy) recall. 
	 the same indexing is adopted for precision, jaccard and dice'''

	import numpy as np

	recall = np.zeros(NumberOfLabels)
	precision = np.zeros(NumberOfLabels)
	dice = np.zeros(NumberOfLabels)
	jaccard = np.zeros( NumberOfLabels)


	ind1 = set(np.argwhere(auto_lbl==2).flatten())
	ind2 = set(np.argwhere(lbl==2).flatten())
	
	dice[0]     = 100*(2*len(ind1 & ind2))/(len(ind1)+len(ind2));
	jaccard[0]   = 100*(len(ind1 | ind2))/len(ind1 | ind2);
	precision[0] = 100*(len(ind1 | ind2))/len(ind1);
	recall[0]    = 100*(len(ind1 | ind2))/len(ind2);

	ind1_a = set(np.argwhere(auto_lbl==1).flatten())
	ind2_a = set(np.argwhere(lbl==1).flatten())

	#get the indexes of voxels with label 1 (non_enhanced tumore)
	ind1_b = set(np.argwhere(auto_lbl==3).flatten())
	ind2_b = set(np.argwhere(lbl==3).flatten())

	ind1 = ind1_a | ind1_b
	ind2 = ind2_a | ind2_b

	dice[1]     = 100*(2*len(ind1 & ind2))/(len(ind1)+len(ind2));
	jaccard[1]   = 100*(len(ind1 | ind2))/len(ind1 | ind2);
	precision[1] = 100*(len(ind1 | ind2))/len(ind1);
	recall[1]    = 100*(len(ind1 | ind2))/len(ind2);

	ind1 = set(np.argwhere(auto_lbl==4).flatten())
	ind2 = set(np.argwhere(lbl==4).flatten())
	
	dice[2]     = 100*(2*len(ind1 & ind2))/(len(ind1)+len(ind2));
	jaccard[2]   = 100*(len(ind1 | ind2))/len(ind1 | ind2);
	precision[2] = 100*(len(ind1 | ind2))/len(ind1);
	recall[2]    = 100*(len(ind1 | ind2))/len(ind2);

	ind1 = set(np.argwhere(auto_lbl > 0).flatten())
	ind2 = set(np.argwhere(lbl > 0).flatten())
	
	dice[3]     = 100*(2*len(ind1 & ind2))/(len(ind1)+len(ind2));
	jaccard[3]   = 100*(len(ind1 | ind2))/len(ind1 | ind2);
	precision[3] = 100*(len(ind1 | ind2))/len(ind1);
	recall[3]    = 100*(len(ind1 | ind2))/len(ind2);
	return(dice, jaccard, precision, recall)