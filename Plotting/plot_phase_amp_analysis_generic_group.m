function plot_phase_amp_analysis_generic_group(R,panlist,obs,phiBin,ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup,ampdata,segdata,rpdata,ampsurdata)
phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
cmap = linspecer(5);

figure(30)
set(gcf,'Position',[1190         330         627         750])
x1 = vertcat(RcoeffGroup{1,:,:});
x2 = vertcat(RcoeffGroup{2,:,:});
xl = -1:.1:0.5; xl= [ones(length(xl),1) xl'];
for i=1:4
    for cond = 1:2
        subplot(4,2,panlist(cond,i))
        if cond == 1; coefs = x1(:,:,i); else;   coefs = x2(:,:,i)'; end
        if size(coefs,2)>size(coefs,1); coefs = coefs'; end
        W = sum(coefs(:,1)~=0)./size(coefs,1);
        coefs(coefs(:,1)==0,:) = NaN;
        coefstore{cond,i}  = coefs;
        if ~isempty(coefs); coefs = (mean(coefs,1).*W)'; else; coefs = [NaN NaN]'; end
        %         plot(xl(:,2),xl*coefs,'color',cmap(i,:),'LineWidth',2)
        a(cond) = gca; uistack(findobj(a(cond).Children,'type','Line'),'top');
    end
    [xydata stats] = getscatterdata(a);
    regstats(:,i,cond) = stats;
    
    statvec(:,i) = coefstats(coefstore{1,i},coefstore{2,i});
end
% statab = table(statvec);
% statab.Properties.RowNames = {'N_{OFF}','N_{ON}', 'Chi_ON_OFF','A','X1bar_SEM_OFF','SEM1_OFF','X1bar_ON','SEM1_ON','X1bar_ON/OFF','t-test P1',...
%     'X2bar_SEM_OFF','SEM2_OFF','X2bar_ON','SEM2_ON','X2bar_ON/OFF','t-test P2'};
% statab.Properties.VariableNames  = {'A','B','C'};
% 
% statvec = [sum(~isnan(x1)) sum(~isnan(y1)) chi2stat(1) p(1) nanmean(x1) nansem2(x1) nanmean(y1) nansem2(y1) nanmean(x1)-nanmean(y1) p(3)...
%     nanmean(x2) nansem2(x2) nanmean(y2) nansem2(y2) nanmean(x2)-nanmean(y2) p(4)];

a = gca;


f = figure(1);
pause(1)
set(f,'Position',[539    81   626   999])

for cond=1:2
    y =  horzcat(ampBinGroup{cond,:,:});

    for i = 1:4

        a = subplot(5,2,panlist(cond,i)); hold on
        y1 = reshape(y(i,:),size(phiBinMid,2),[])';
        ysave{cond,i} = y1;
        [hl hp] = boundedline(phiBinMid', nanmedian(y1)',abs([prctile(y1,16)' prctile(y1,84)']-nanmedian(y1)'),'cmap',cmap(i,:)); %
        hl.LineWidth = 2;
        L = get(gca,'children'); L(2) = L(end); L(end) = hp;set(gca,'children',L)
        
        ampSur = meanSurData(ampsurdata,cond,i)
        plot(phiBinMid,squeeze(ampSur(:,1)),'k--')
        plot(phiBinMid,squeeze(ampSur(:,2)),'k--')
        
        xlim([-pi pi])
        if cond ==2
            figure(f)
            a = subplot(5,2,panlist(cond,i)); hold on
            yfit = vertcat(ysave{:,i})
            x1 = sqres(ysave{1,i},nanmedian(yfit))/1000; 
            x2 = sqres(ysave{2,i},nanmedian(yfit))/1000;
            p(i) = ranksum(x1,x2);
            [figx figy] = dsxy2figxy(a, -3,a.YLim(1)*0.95)
            g = annotation(gcf,'textbox',...
                [figx figy 0.3 0.03],...
                'String',sprintf('MW-U test = %.3f',p(i)),...
                'LineStyle','none',...
                'FontSize',8,...
                'FontWeight','bold',...
                'FontAngle','italic',...
                'FitBoxToText','off');
            for p = 1:2
                [figx figy] = dsxy2figxy(a, -3,a.YLim(1)*0.95);
                g.Position = [figx figy 0.3 0.03];
            end
        end
    end
end
cmapA = linspecer(5);
cmap = linspecer(5);
% cmapA(4,:) = cmap(5,:);
panlist(:,5) = [9; 10];

for cond=1:2
    y =  vertcat(segBinGroup{cond,:,:});
    ysave{cond,5} = y;
    subplot(5,2,panlist(cond,5)); hold on
    [hl hp] = boundedline(phiBinMid', nanmedian(y)',abs([prctile(y,16)' prctile(y,84)']-nanmedian(y)'),'cmap',cmap(5,:));
    hl.LineWidth = 2;
    L = get(gca,'children'); L(2) = L(end); L(end) = hp; set(gca,'children',L)
    xlim([-pi pi])
    if cond ==2
        x1 = sqres(ysave{1,5},nanmedian(ysave{1,5}))/1000; x2 = sqres(ysave{2,5},nanmedian(ysave{2,5}))/1000;
        p(4) = ranksum(x1,x2);
        a = gca;[figx figy] = dsxy2figxy(gca, -3,a.YLim(1)*0.95);
        annotation(gcf,'textbox',...
            [figx figy 0.3 0.03],...
            'String',sprintf('RS test = %.3f',p(4)),...
            'LineStyle','none',...
            'FontSize',8,...
            'FontWeight','bold',...
            'FontAngle','italic',...
            'FitBoxToText','off');
    end
end
clear p
[ptt,pvart] = cagnan_stats_test(phiBinMid,ysave,1,panlist,cmapA);
annotation(gcf,'textbox',...
    [ 0.1820    0.9570    0.2550    0.0350],...
    'String',R.condname{1},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.583 0.957 0.304 0.0349],...
    'String',R.condname{2},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off');
%     savefigure_v2([R.resultspathr '\Group\seganalysis\'],...
%         [R.bregname{breg} '_Group_CagnanAnalysis'],[],[],'-r100');
%     close all
clear ampBinGroup segBinGroup

% % % figure
% % % set(gcf,'Position',[823   421   785   677])
% % % x1 = vertcat(phipeakGroup{1,:,:});
% % % x2 = vertcat(phipeakGroup{2,:,:});
% % % 
% % % legloc = legloc();
% % % for i = 1:5
% % %     subplot(3,2,i)
% % %     x11 = x1(:,i);
% % %     x22 = x2(:,i);
% % %     polarhistogram(x11,6,'Normalization','Probability','FaceColor',cmapA(i,:).*0.75,'LineStyle','none','FaceAlpha',.5); hold on
% % %     a(i,1) =  polarhistogram(x11,6,'Normalization','Probability','DisplayStyle','stairs','EdgeColor',cmapA(i,:).*0.75,'LineWidth',1.2);
% % %     polarhistogram(x22,6,'Normalization','Probability','FaceColor',cmapA(i,:),'LineStyle','none','FaceAlpha',.3); hold on
% % %     a(i,2) = polarhistogram(x22,6,'Normalization','Probability','DisplayStyle','stairs','EdgeColor',cmapA(i,:),'LineWidth',1.5);
% % %     title(['RP Dist. for Peak ' obs{i}]); rlim([0 0.5])
% % %     legend(a(i,:),R.condname,'Position',legloc(i,:),'Orientation','Horizontal','Color','none','EdgeColor','none')
% % %     %         [p,m,pv] = circ_cmtest([x11;x22]',[repmat(1,size(x11));repmat(2,size(x22))]');
% % %     %circ_wwtest([x11;x22],[repmat(1,size(x11));repmat(2,size(x22))]); % ttest
% % %     p(i,1) = circ_rtest(x11);rayleg(p(i,1),i,1);
% % %     p(i,2) = circ_rtest(x22);rayleg(p(i,2),i,2);
% % %     %         [sig(i)] = disp_anova(x11,x22)
% % % end
%     savefigure_v2([R.resultspathr '\Group\seganalysis\'],...
%         [R.bregname{breg} '_Group_CagnanAnalysis_RP_rose'],[],[],'-r100'); close all


function legloc = legloc();
legloc = [0.1710    0.5210    0.2500    0.0300;
    0.6110    0.5210    0.2500    0.0300;
    0.1710    0.0530    0.2500    0.0300;
    0.6110    0.0530    0.2500    0.0300];
function rayleg(P,i,cond)
posbank = {[0.1500    0.5000    0.1560    0.0250; 0.3400    0.5000    0.1560    0.0250];
    [0.5800    0.5000    0.1560    0.0250; 0.7700    0.5000    0.1560    0.0250];
    [0.1700    0.0210    0.1560    0.0250; 0.3400    0.0210    0.1560    0.0250];
    [0.5800    0.0210    0.1560    0.0250; 0.7700    0.0210    0.1560    0.0250]};

annotation(gcf,'textbox',...
    posbank{i}(cond,:),...
    'String',sprintf('Ray. P = %.3f',P),...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontAngle','italic',...
    'FitBoxToText','off');

%%[sprintf('OFF Res. = %.1f ',nanmean(x1)) sprintf('ON Res. = %.1f ',nanmean(x2))];
function statvec = coefstats(x,y)
x1 = x(:,1); y1 = y(:,1);
x2 = x(:,2); y2 = y(:,2);

[h,p(1), chi2stat(1),df] = prop_test([sum(~isnan(x1)) sum(~isnan(y1))], [size(x1,1),size(y1,1)],'false')
% [h,p(2), chi2stat(2),df] = prop_test([sum(~isnan(y1))], [size(y1,1),size(y1,1)],'false')
try
    [h p(3)]=ttest2(x1(:,1),y1(:,1))
    [h p(4)]=ttest2(x2(:,1),y2(:,1))
catch
    p(3:4) = 1;
end

statvec = [sum(~isnan(x1)) sum(~isnan(y1)) chi2stat(1) p(1) nanmean(x1) nansem2(x1) nanmean(y1) nansem2(y1) nanmean(x1)-nanmean(y1) p(3)...
    nanmean(x2) nansem2(x2) nanmean(y2) nansem2(y2) nanmean(x2)-nanmean(y2) p(4)];

function nansem2 = nansem2(x)
nansem2 = nanstd(x)./sqrt(sum(~isnan(x(:))))

function [xydata stats] = getscatterdata(a)
for cond = 1:2
    scatob = findobj(a(cond).Children,'type','scatter');
    for i = 1:numel(scatob)
        xydata{i} = [scatob(i).XData; scatob(i).YData; repmat(i,1,size(scatob(i).YData,2)); repmat(cond,1,size(scatob(i).YData,2))]';
    end
    xydatac{cond} = vertcat(xydata{:});
end
xydata = vertcat(xydatac{:})
X = xydata(:,1); Y = xydata(:,2); I = xydata(:,3); COND = xydata(:,4);
tbl = table(xydata(:,1),xydata(:,2),xydata(:,3),xydata(:,4),'VariableNames',{'X','Y','I','COND'});
lme1 = fitlme(tbl,'Y~X+(X|I)');
lme2 = fitlme(tbl,'Y~ X + COND +(X|I)+(X|COND)');
% results = compare(lme1,lme2,'nsim',500)
stats = anova(lme2);
stats = [stats.FStat(2) stats.pValue(2) stats.FStat(3) stats.pValue(3)];
yhat = predict(lme2,tbl);
xyhat = [X yhat I COND];
for cond = 1:2
    axes(a(cond))
    for i = 1:numel(scatob)
        xy = xyhat(xyhat(:,3) == i & xyhat(:,4) == cond,1:2);
        xy = sort(xy,1);
        plot(xy(:,1),xy(:,2),'k:'); hold on
    end
end


function X = meanSurData(ampsurdata,cond,i)
X1 = []; X2 = [];
for side=1:2
    for sub = 1:size(ampsurdata,3)
        X1 = [X1 squeeze(ampsurdata{cond,side,sub}(:,i,1))]
        X2 = [X2 squeeze(ampsurdata{cond,side,sub}(:,i,2))]
    end
end
X1 = mean(X1,2); X2 = mean(X2,2);
X =[X1 X2];
%
% w = linspace(min(X),max(X));
%
% figure()
% gscatter(X,Y,I,'bgr','x.o')
% line(w,feval(lme,w,1),'Color','b','LineWidth',2)
% line(w,feval(lme,w,2),'Color','g','LineWidth',2)
% line(w,feval(lme,w,3),'Color','r','LineWidth',2)
% title('Fitted Regression Lines by Model Year')