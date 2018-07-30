function [BinMu,BinSEM] = binstats(PhiCol,SegCol,Bin)
for j = size(SegCol,2)
for i = 1:length(Bin)-1
    binDat = SegCol(:,PhiCol>=Bin(i) & PhiCol<=Bin(i+1))';
    binNumel(i) = numel(binDat);
    BinMu(:,i) = nanmean(binDat);%*(numel(binDat)/numel(SegCol));
    BinSEM(:,i) = nanstd(binDat)/size(binDat,1);
end
BinMu = BinMu.*(binNumel./max(binNumel));
end