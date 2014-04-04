function lfFiltered = approxMedfilt3(lf,WinSizeRow,WinSizeCol,WinSizeDepth)

lfFiltered = lf;

for d=1:nbDepths
    lfFiltered(:,:,d) = medfilt2(lf(:,:,d),[WinSizeRow WinSizeCol]);
end

for r=1:nbRows
    tmp = squeeze(lfFiltered(r,:,:));
    tmp = medfilt2(tmp,[1 WinSizeDepth]);
    lfFiltered(r,:,:) = tmp;
end
