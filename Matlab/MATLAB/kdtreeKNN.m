function [ind,knndist] = kdtreeKNN(feature_training,label_training,feature_total,K)

    kdtree = vl_kdtreebuild(feature_training);
    [ind,knndist] = vl_kdtreequery(kdtree, feature_training, feature_total,'NUMNEIGHBORS',K);
    
    for i=1:size(ind,2)
        for k=1:K
            ind(k,i) = label_training(ind(k,i));
        end
    end

end