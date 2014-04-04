function matrix = make_dataterm_matrix_classifier(score,mask_idx,space,T1C)

%k = size(IDX,2);

nlabels = size(score,2);

matrix = zeros(nlabels,length(T1C(:)));
matrix(1,:) = 1;

[h,w,d] = size(T1C);
for count = 1: length(score(:,1))

        ROW = uint32(space(4,count) * h);
        COL = uint32(space(5,count) * w);
        DEP = uint32(space(6,count) * d);
       % importance(ROW,COL,DEP) = min(D(count,:));
        index = (h*w)*(DEP-1)+(COL - 1)*h+ROW;

%         M = mask_idx(IDX(count,:));
%         M(M==10)=0;
%         M(M==3)=1;
%         M(M==4)=3;
   
        matrix(:,index) = score(count,:)/sum(score(count,:));
    %    end
end 



