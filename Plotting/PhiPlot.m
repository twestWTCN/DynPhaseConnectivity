function a = PhiPlot(phiGroup,phibinmid,cmap,R)
for condr = 1:2
    x2 = vertcat(phiGroup{:,:,condr});
    for ii = 1:size(x2,1)
        a(2*condr - 1) = polarscatter(phibinmid,(x2(ii,:)),25,cmap(condr,:),'Marker','o'); hold on
    end
    mat = mean(x2,1); mat(end+1) = mat(1);
    phibinmid2 = phibinmid;
    phibinmid2(end+1) = phibinmid2(1);
    a(2*condr) = polarplot(phibinmid2,(mat),'color',cmap(condr,:),'LineWidth',2)
end
pa = gca;
pa.ThetaAxisUnits = 'radians';
pa = legend(a,{[R.condname{1} ' Data'],[R.condname{1} ' Mean'],[R.condname{2} ' Data'],[R.condname{2} ' Mean']});
set(pa,'Position',[0.7830    0.7660    0.1950    0.1650])
title('Distribution of Group Level Relative Phase')