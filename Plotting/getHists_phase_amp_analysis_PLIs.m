function dwell = getHists_phase_amp_analysis_PLIs(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
close all
condcr = {'r','b'};
for breg = 1:length(R.bregname)
    for sub = 1:length(R.subname)
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences
                    
                    segL = vc_clean.PA.segL_pli_dist_save;
                    x = vc_clean.PA.timevec{1}(end);
                    dwellrat(side,cond,sub) = sum(segL)/x;
                    
                    figure(10)
                    subplot(1,2,cond)
                    phase_ang = vc_clean.PA.pA_pli_dist_save;
                    histogram(phase_ang,linspace(-pi,pi,6),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
                    xlabel('Phase Angle');ylabel('P(X)'); ylim([0 0.5]); grid on; title(R.condname{cond})
                    [N B] = histcounts(phase_ang,linspace(-pi,pi,6));
                    phaseang_dist(cond,sub).N = N;
                    phaseang_dist(cond,sub).B = B;
                    
                    figure(20)
                    subplot(1,2,cond)
                    histogram(segL,logspace(log10(0.01),log10(6),16),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
                    set(gca,'xscale','log')
                    xlabel('Segment Length');ylabel('P(X)');  grid on; title(R.condname{cond});ylim([0 0.5]); xlim([0.01 6]);
                    [N B] = histcounts(segL,linspace(0,3.5,12));
                    segL_dist(cond,sub).N = N;
                    segL_dist(cond,sub).B = B;
                    
                    ampname = {'HB CTX','HB STN','LB STN'};
                    for i = 1:3
                        amp = vc_clean.PA.amp_pli_dist_save(i,:);
                        figure(30+i)
                        subplot(1,2,cond)
                        histogram(amp,linspace(-100,600,12),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
                        xlabel(ampname{i});ylabel('P(X)'); ylim([0 0.7]); grid on; title(R.condname{cond})
                        [N B] = histcounts(amp,linspace(-100,100,14));
                        amp_dist(i,cond,sub).N  = N;
                        amp_dist(i,cond,sub).B  = B;
                    end
                    
                    H = vc_clean.PA.H_dist_save{cond}(1,:);
                    figure(40)
                    subplot(1,2,cond)
                    histogram(H,linspace(-1,1,12),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
                    xlabel('Amp Env Correlation');ylabel('P(X)'); ylim([0 0.3]); grid on; title(R.condname{cond})
                    [N B] = histcounts(segL,linspace(0,1,14));
                    H_dist(cond,sub).N = N;
                    H_dist(cond,sub).B = B;
                    Hcorr(side,cond,sub) = nanmean(abs(H(1,:)));
                    
                    figure(252)
                    subplot(1,2,cond)
                    y =  vc_clean.PA.H_dist_save{cond}(1,:); x = vc_clean.PA.segL_pli_dist_save;
                    scatter(x,y(1,:),'filled','MarkerFaceColor',condcr{cond},'MarkerFaceAlpha',0.1);hold on;[r1 p1] = corrcoef(x,y(1,:)); length_amp_corr(:,cond,sub) = [r1(2) p1(2)];
                    xlabel('Segment Length'); ylabel('Amplitude Correlation'); grid on; title(R.condname{cond})
                    set(gca,'xscale','log'); xlim([0.2 1.5]); ylim([-1 1])
                    for i=1:3
                        figure(260+i)
                        subplot(1,2,cond)
                        y =  vc_clean.PA.H_dist_save{cond}(1,:); x = vc_clean.PA.amp_pli_dist_save(i,:); %segL_pli_dist_save{cond};
                        scatter(x,y(1,:),'filled','MarkerFaceColor',condcr{cond},'MarkerFaceAlpha',0.1);hold on;[r1 p1] = corrcoef(x,y(1,:)); amp_amp_corr(:,i,cond,sub) = [r1(2) p1(2)];
                        xlabel(ampname{i}); ylabel('Amplitude Correlation'); grid on; title(R.condname{cond});xlim([-100 250]); ylim([-1 1])
                    end
                end
            end
        end
    end
    % save([R.datapathr '\results\seganalysis\groupseganaly'],'ampsave','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
    
    
    x = (squeeze(dwellrat(:,1,:))); x(x==0) = [];
    y = (squeeze(dwellrat(:,2,:))); y(y==0) = [];
    
    figure(1)
    subplot(1,2,1)
    histogram(x(:),linspace(0.01,0.3,12),'FaceColor',condcr{1},'Normalization','Probability'); %xlim([0 1])
    xlabel('Dwell/Escape'); ylabel('P(X)'); grid on; ylim([0 0.7])
    subplot(1,2,2);
    histogram(y(:),linspace(0.01,0.3,12),'FaceColor',condcr{2},'Normalization','Probability'); %xlim([0 1]);
    xlabel('Dwell/Escape'); ylabel('P(X)'); grid on; ylim([0 0.7])
    clear p
    [h p] = ttest2(x,y);
    title(num2str(p))
    dwell = {x,y,p};
    
    x = (squeeze(Hcorr(:,1,:))); x(x==0) = [];
    y = (squeeze(Hcorr(:,2,:))); y(y==0) = [];
    
    figure(2)
    subplot(1,2,1)
    histogram(x,linspace(0,1,12),'FaceColor',condcr{1},'Normalization','Probability'); xlim([0 1])
    xlabel('STN/M1 Amp Corr'); ylabel('P(X)'); grid on; ylim([0 0.7])
    subplot(1,2,2);
    histogram(y,linspace(0,1,12),'FaceColor',condcr{2},'Normalization','Probability'); xlim([0 1]);
    xlabel('STN/M1 Amp Corr'); ylabel('P(X)'); grid on; ylim([0 0.7])
    [h p] = ttest2(x,y);
    title(num2str(p))
    dwell = {x,y,p};
    close all
%     savefigure_v2([R.datapathr '\results\seganalysis\PLI\partests\'],[ num2str(idd) '_PLI_surrogate_bwid_' num2str(R.PA.bwid) '_PLVeps_' num2str(R.PA.PLVeps)],[],[],[]); close all
    
    % figure(3)
    % x = sum(squeeze(length_amp_corr(1,1,:,:)),3);x(x==0) = [];
    % y = sum(squeeze(length_amp_corr(1,2,:,:)),3);y(y==0) = [];
    % subplot(1,2,1)
    % histogram(x,linspace(-1,1,12),'FaceColor',condcr{1},'Normalization','Probability');% xlim([0 1])
    % xlabel('Frame Length'); ylabel('Amplitude Correlation'); grid on;% ylim([0 0.7])
    % subplot(1,2,2);
    % histogram(y,linspace(-1,1,12),'FaceColor',condcr{2},'Normalization','Probability'); %xlim([0 1]);
    % xlabel('Frame Length'); ylabel('Amplitude Correlation'); grid on; %ylim([0 0.7])
    % [h p] = ttest2(x,y)
    % title(num2str(p))
    
    %
    % for i = 1:3
    % figure(400+i)
    % x = sum(squeeze(amp_amp_corr(1,i,1,:,:)),3);x(x==0) = [];
    %     y = sum(squeeze(amp_amp_corr(1,i,2,:,:)),3);y(y==0) = [];
    %     subplot(1,2,1)
    %     histogram(x,linspace(-1,1,12),'FaceColor',condcr{1},'Normalization','Probability');% xlim([0 1])
    %     xlabel(ampname{i}); ylabel('Amplitude Correlation'); grid on;% ylim([0 0.7])
    %     subplot(1,2,2);
    %     histogram(y,linspace(-1,1,12),'FaceColor',condcr{2},'Normalization','Probability'); %xlim([0 1]);
    %     xlabel(ampname{i}); ylabel('Amplitude Correlation'); grid on; %ylim([0 0.7])
    %     [h p] = ttest2(x,y)
    %     title(num2str(p))
    % end
    % savefigure_v2([R.datapathr '\results\seganalysis\PLI\partests\'],[ num2str(idd) '_PLI_surrogate_.tiff'],[],[],[]); close all
    
    a = 1
end