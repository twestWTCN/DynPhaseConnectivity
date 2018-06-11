function [ptt,pvart] = cagnan_stats_test(x,ysave,plotop,panlist,cmap)

for i = 1:size(ysave,2)
    x1 = ysave{1,i};
    x2 = ysave{2,i};
    
    [h ptt(:,i)] =ttest2(x1,x2);
    [h pvart(:,i)] =vartest2(x1,x2);
end

if plotop == 1
    for cond=1:size(ysave,1)
        for i = 1:size(ysave,2)
            alpha = 0.05/size(ptt,1);
            a =  subplot(4,2,panlist(cond,i));
            
            mask = pvart(:,i)<=alpha;
            plotSigBar(mask,x,a,cmap(i,:),'-',0.70,pvart(:,i))
            
            mask = ptt(:,i)<=alpha;
            plotSigBar(mask,x,a,cmap(i,:),':',0.85,ptt(:,i))
            
        end
    end
end

function    plotSigBar(mask,x,a,cmap,ls,yoff,pvart)
consecSegs = SplitVec((mask'.*(1:numel(mask))),'consecutive');
for p=1:numel(consecSegs)
    if any(consecSegs{p})>0
        sigx = x(consecSegs{p});
        dx = (x(2)-x(1))/2;
        sigx = min(sigx)-dx:dx:max(sigx)+dx;
        sigx = sigx(sigx>=min(x) & sigx<=max(x));
        plot(sigx,repmat(a.YLim(2)*yoff,size(sigx)),'color',cmap,'LineWidth',5,'LineStyle',ls)
        [figx figy] = dsxy2figxy(gca, median(sigx)-pi/6,a.YLim(2)*yoff);
        figy = figy-0.01;
        annotation(gcf,'textbox',...
            [figx figy 0.3 0.03],...
            'String',sprintf('P = %.3f   ',median(pvart(consecSegs{p}))),...
            'LineStyle','none',...
            'FontSize',6,...
            'FontWeight','bold',...
            'FontAngle','italic',...
            'FitBoxToText','off');
    end
end

function annotP(P,i,cond)
posbank = {[0.1500    0.5000    0.1560    0.0250; 0.3200    0.5000    0.1560    0.0250];
    [0.5700    0.5000    0.1560    0.0250; 0.7400    0.5000    0.1560    0.0250];
    [0.1700    0.0210    0.1560    0.0250; 0.3400    0.0210    0.1560    0.0250];
    [0.5900    0.0210    0.1560    0.0250; 0.7600    0.0210    0.1560    0.0250]};

annotation(gcf,'textbox',...
    posbank{i}(cond,:),...
    'String',sprintf('Ray. P = %.3f',P),...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontAngle','italic',...
    'FitBoxToText','off');
