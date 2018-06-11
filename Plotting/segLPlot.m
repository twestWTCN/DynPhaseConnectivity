function a = segLPlot(segLGroup,segbinmid,cmap,R)
x1 = vertcat(segLGroup{:,:,1});
for i = 1:size(x1,1)
    a(1) = scatter(segbinmid,x1(i,:),25,cmap(1,:),'Marker','o'); hold on
end
f = fit(repmat(segbinmid',size(x1,1),1),reshape(x1',1,[])','gauss1');
a(2) = plot(f) ; set(a(2),'color',cmap(1,:),'LineWidth',2);
x1 = vertcat(segLGroup{:,:,2});
for i = 1:size(x1,1)
    a(3) = scatter(segbinmid,x1(i,:),25,cmap(2,:),'Marker','x'); hold on
end
f = fit(repmat(segbinmid',size(x1,1),1),reshape(x1',1,[])','gauss1');
a(4) = plot(f) ; set(a(4),'color',cmap(2,:),'LineWidth',2);
legend(a,{[R.condname{1} ' Data'],[R.condname{1} ' Fit'],[R.condname{2} ' Data'],[R.condname{2} ' Fit']})
grid on
ylim([0 0.25]);
xlabel('log10 Frame Length (s)')
ylabel('Occurence Rate (s^-1)')
title('Distribution of Group Level Frame Lengths')