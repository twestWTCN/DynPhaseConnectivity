function plotExampleSub_PA(phi,phiBin,ampSegCol,segLCol,cmap)
phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
phiBinMidSm = linspace(phiBinMid(1),phiBinMid(end),24);
for i = 1:3
    subplot(3,1,i)
    [BinMu,BinSEM,binDat,BinPr2] = binstats(phi,ampSegCol(i,:),phiBin);
    Xs = [spline(phiBinMid,BinPr2(1,:),phiBinMidSm); spline(phiBinMid,BinPr2(3,:),phiBinMidSm); spline(phiBinMid,BinPr2(2,:),phiBinMidSm)];
    [hl, hp] = boundedline(phiBinMidSm', Xs(1,:)',[(Xs(1,:)-Xs(3,:));(Xs(2,:)-Xs(1,:))]','cmap',cmap(i,:));  hold on
    hl.LineWidth = 2;
    scatter(phi,ampSegCol(i,:),segLCol.*50,cmap(i,:),'filled')
    xlim([-pi pi]); grid on; xlabel('Relative Phase (rads)'); ylabel('Obs')
end
set(gcf,'Position',[680.0000   51.5000  449.0000  926.5000])
