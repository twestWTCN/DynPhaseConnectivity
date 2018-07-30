function [BinMu,BinSEM] = binstats_v2(PhiCol,SegCol,Bin)
for j = 1:size(SegCol,2)
    for i = 1:length(Bin)-1
        binDat = SegCol(PhiCol>=Bin(i) & PhiCol<=Bin(i+1),j);
        binNumel(i,j) = numel(binDat);
        BinMu(:,i,j) = nanmean(binDat);%*(numel(binDat)/numel(SegCol));
        BinSEM(:,i,j) = nanstd(binDat)/size(binDat,1);
    end
    
    BinMu(:,:,j) = BinMu(:,:,j).*(binNumel(:,j)./max(binNumel(:,j)))';
end
