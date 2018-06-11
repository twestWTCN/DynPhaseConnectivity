function plot_phase_amp_analysis_generic_group(R,panlist,obs,phiBin,ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup)
    phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);


cmap = linspecer(3);

    figure(30)
    set(gcf,'Position',[1190         330         627         750])
    x1 = vertcat(RcoeffGroup{1,:,:});
    x2 = vertcat(RcoeffGroup{2,:,:});
    xl = -1:.1:0.5; xl= [ones(length(xl),1) xl'];
    for cond = 1:2
        for i=1:3
            subplot(3,2,panlist(cond,i))
            if cond == 1; coefs = x1(:,:,i); else;   coefs = x2(:,:,i)'; end
            if size(coefs,2)>size(coefs,1); coefs = coefs'; end
            W = sum(coefs(:,1)~=0)./size(coefs,1);
            coefs(coefs(:,1)==0,:) = [];
            if ~isempty(coefs); coefs = (mean(coefs,1).*W)'; else; coefs = [0 0]'; end
            plot(xl(:,2),xl*coefs,'color',cmap(i,:),'LineWidth',2)
            a = gca; uistack(findobj(a.Children,'type','Line'),'top');
        end
    end
    a = gca;
        mean(x1)
        mean(x2)
    [h,p, chi2stat,df] = prop_test([sum(x2(:,1:2)>0)], [size(x2,1),size(x2,1)],'false');
    [h,p, chi2stat,df] = prop_test([sum(x1(:,1:2)>0)], [size(x1,1),size(x1,1)],'false');
    [h p]=ttest2(x1(:,2),x2(:,2));
    [h p]=ttest2(x1(:,1),x2(:,1));
    
    f = figure(1);
    pause(1)
    set(f,'Position',[539    81   626   999])
    
    for cond=1:2
        y =  horzcat(ampBinGroup{cond,:,:});
        for i = 1:3
            a = subplot(4,2,panlist(cond,i)); hold on
            y1 = reshape(y(i,:),size(phiBinMid,2),[])';
            ysave{cond,i} = y1;
            [hl hp] = boundedline(phiBinMid', nanmedian(y1)',abs([prctile(y1,16)' prctile(y1,84)']-nanmedian(y1)'),'cmap',cmap(i,:)); %
            hl.LineWidth = 2;
            L = get(gca,'children'); L(2) = L(end); L(end) = hp;set(gca,'children',L)
            xlim([-pi pi])
            if cond ==2
                figure(f)
                a = subplot(4,2,panlist(cond,i)); hold on
                x1 = sqres(ysave{1,i},nanmedian(ysave{1,i}))/1000; x2 = sqres(ysave{2,i},nanmedian(ysave{2,i}))/1000;
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
    cmapA = linspecer(3);
    cmap = linspecer(5);
    cmapA(4,:) = cmap(5,:);
    panlist(:,4) = [7; 8];

    for cond=1:2
        y =  vertcat(segBinGroup{cond,:,:});
        ysave{cond,4} = y;
        subplot(4,2,panlist(cond,4)); hold on
        [hl hp] = boundedline(phiBinMid', nanmedian(y)',abs([prctile(y,16)' prctile(y,84)']-nanmedian(y)'),'cmap',cmap(5,:));
        hl.LineWidth = 2;
        L = get(gca,'children'); L(2) = L(end); L(end) = hp; set(gca,'children',L)
        xlim([-pi pi])
        if cond ==2
            x1 = sqres(ysave{1,4},nanmedian(ysave{1,4}))/1000; x2 = sqres(ysave{2,4},nanmedian(ysave{2,4}))/1000;
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
    
    figure
    set(gcf,'Position',[823   421   785   677])
    x1 = vertcat(phipeakGroup{1,:,:});
    x2 = vertcat(phipeakGroup{2,:,:});
    
    legloc = legloc();
    for i = 1:4
        subplot(2,2,i)
        x11 = x1(:,i);
        x22 = x2(:,i);
        polarhistogram(x11,6,'Normalization','Probability','FaceColor',cmapA(i,:).*0.75,'LineStyle','none','FaceAlpha',.5); hold on
        a(i,1) =  polarhistogram(x11,6,'Normalization','Probability','DisplayStyle','stairs','EdgeColor',cmapA(i,:).*0.75,'LineWidth',1.2);
        polarhistogram(x22,6,'Normalization','Probability','FaceColor',cmapA(i,:),'LineStyle','none','FaceAlpha',.3); hold on
        a(i,2) = polarhistogram(x22,6,'Normalization','Probability','DisplayStyle','stairs','EdgeColor',cmapA(i,:),'LineWidth',1.5);
        title(['RP Dist. for Peak ' obs{i}]); rlim([0 0.5])
        legend(a(i,:),R.condname,'Position',legloc(i,:),'Orientation','Horizontal','Color','none','EdgeColor','none')
        %         [p,m,pv] = circ_cmtest([x11;x22]',[repmat(1,size(x11));repmat(2,size(x22))]');
        %circ_wwtest([x11;x22],[repmat(1,size(x11));repmat(2,size(x22))]); % ttest
        p(i,1) = circ_rtest(x11);rayleg(p(i,1),i,1);
        p(i,2) = circ_rtest(x22);rayleg(p(i,2),i,2);
        %         [sig(i)] = disp_anova(x11,x22)
    end
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
