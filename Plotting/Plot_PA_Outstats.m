function Plot_PA_Outstats(R,RP,SegStat,OVL,RPSurr,SegLsurr,OVLsurr,obs,titlist,cond)
panlist = [1 3 5 7 9 ; 2 4 6 8 10];
cmap = linspecer(5);

            
            %% Duration Histograms
            figure(1)
            subplot(2,1,cond)
            histogram(SegStat(:,5),logspace(log10(1e-2),log10(1.5),24),'Normalization','probability'); hold on
            histogram(SegLsurr,logspace(log10(1e-2),log10(1.5),24),'Normalization','probability')
            legend({'Data','Surr.'}); title(R.condnames{cond}); xlabel(['log ' obs{5}]); ylabel('P(x)')
            set(gca,'XScale', 'log'); grid on
            set(gcf,'Position',[654   324   358   750])
            
            %% RP vs Amp
            figure(2)
            sizeheads(R)
            for i = 1:size(SegStat,2)
                subplot(size(SegStat,2),2,panlist(cond,i))
                hold on
                scatter(RP,SegStat(:,i),SegStat(:,5).*10,cmap(i,:),'filled'); grid on
                xlim([-pi pi]); xlabel('Relative Phase')
                title(titlist{i}); ylabel(obs{i});
                if i<5
                %    ylim([0 200])
                else
                    %                 ylim([0 200])
                end
            end
            
            %% Seg Dur vs Amp
            figure(3)
            sizeheads(R)
            for i = 1:size(SegStat,2)-1
                subplot(size(SegStat,2),2,panlist(cond,i))
                scatter(SegStat(:,5),SegStat(:,i),50,cmap(i,:),'filled'); grid on
                [A p] = linplot(SegStat(:,5),SegStat(:,i),0.05)
                if p<0.05; text(0.15,150,sprintf('P = %.2f',p)); end
                title([obs{i} ' vs ' obs{5}]); ylabel(obs{i}); xlabel(obs{5});
               % ylim([0 200]); xlim([0 0.25])
            end

%% RP Rose Histograms
            figure(4)
            subplot(2,1,cond)
            b1 = polarhistogram(RP,linspace(-pi,pi,24),'Normalization','probability'); hold on
            [p z] = circ_rtest(RP); rlim([0 0.3])
                text(60,0.55,{'Ray Test:'; sprintf('P = %.3f',p)},'FontSize',10)
            %             b1.Values = b1.Values/TotTime;
            b2 = polarhistogram(RPSurr,linspace(-pi,pi,18),'Normalization','probability');
            %              b2.Values = b1.Values/TotTimesurr;
            title(R.condnames{cond})
            if cond ==2
                leg1 = legend({'Data','Surr.'})
                set(leg1,...
                    'Position',[754   369   354   729]);
            end
            set(gcf,'Position',[1190         330         627         750])
            
            
function [] = sizeheads(R)
set(gcf,'Position',[1190         330         627         750])
annotation(gcf,'textbox',...
    [ 0.1820    0.9570    0.2550    0.0350],...
    'String',R.condnames{1},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.583 0.957 0.304 0.0349],...
    'String',R.condnames{2},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off');