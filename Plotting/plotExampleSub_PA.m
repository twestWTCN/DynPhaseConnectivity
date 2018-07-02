    function plotExampleSub_PA(phi,phiBin,ampSegCol,segLCol,ampSur,cmap,obs,ylimlist,titlist,cond)
phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
phiBinMidSm = linspace(phiBinMid(1),phiBinMid(end),24);
% phiBinMidSm = phiBinMid;
for i = 1:4
    subplot(1,4,i)
    [BinMu,BinSEM,binDat,BinPr2] = binstats(phi,ampSegCol(i,:),phiBin);
%     Xs = [circ_interp(phiBinMid,BinMu,phiBinMidSm); circ_interp(phiBinMid,BinPr2(3,:),phiBinMidSm); circ_interp(phiBinMid,BinPr2(2,:),phiBinMidSm)];
%     [hl, hp] = boundedline(phiBinMidSm', Xs(1,:)',[(Xs(1,:)-Xs(3,:));(Xs(2,:)-Xs(1,:))]','cmap',cmap(i,:));  hold on
    hl.LineWidth = 2;
    scatter(phi,ampSegCol(i,:),segLCol.*100,cmap(i,:),'filled')
%     plot(phiBinMidSm,circ_interp(phiBinMid,squeeze(ampSur(:,i,1))',phiBinMidSm),'k--')
%     plot(phiBinMidSm,circ_interp(phiBinMid,squeeze(ampSur(:,i,2))',phiBinMidSm),'k--')
    
    
    xlim([-pi pi]); grid on; xlabel('Relative Phase (rads)'); ylabel(obs{i});
    title(titlist{i});
    ylim(ylimlist{i})
end
if cond ==1
    set(gcf,'Position',[393         548        1438         289])
else
    set(gcf,'Position',[391         170        1438         289])
end


% y2 = y1;
