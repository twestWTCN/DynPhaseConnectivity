function dwell = PA_Amp_DeAmp_frqvar(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
close all
QX = 8; Lmark = 0;
for breg = 2:length(R.bregname)
    for sub = 1:length(R.subname)
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences
                    Lmark = Lmark+1;
                    
                    obs = {[R.bregname{breg} ' ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{R.bregband{breg}}],['SMA ' R.bandinits{2}],['STN ' R.bandinits{2}]};
                    for i = 1:4
                        titlist{i} = ['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{i} ' Power'];
                    end
                    
                    for frqi = 1:size(vc_clean.PA.frqvar.timevec,2)
                        tottime = vc_clean.PA.frqvar.timevec{frqi}(end)-vc_clean.PA.frqvar.timevec{frqi}(1);
                        
                        pureAmp = vc_clean.PA.frqvar.pureAmp{frqi}';
                        amp_dist = vc_clean.PA.frqvar.amp_pli_dist_save{frqi};
                        seg_dist = vc_clean.PA.frqvar.segL_pli_dist_save{frqi};
                        rp_dist = vc_clean.PA.frqvar.pA_pli_dist_save{frqi};
                        surrAmp = vc_clean.PA.frqvar.surrAmp{frqi}(end);
                        % First Establish Relationship between PL and Amplitude
                        for i = 1:4
                            ampi = pureAmp(i,:); % Find the amplitude envelope
                            amp_segs = amp_dist(i,:); % Find amplitude within segs
                            thr = [prctile(ampi,25) prctile(ampi,75)]; % Treat at deamp or amp
                            %                         subplot(2,2,i)
                            %                         hist(ampi,-100:25:200)
                            amp_inds = find(amp_segs<=thr(1));
                            ampVal(1,i,cond,side,sub) = median(amp_segs(amp_inds));
                            ampDurfrq(frqi,i,1) = (sum(seg_dist(amp_inds))/tottime)*100;
                            ampDur(1,i,cond,side,sub) = ampDurfrq(frqi,i,1);
                            
                            x_draw = amp_dist(i,:);
                            draw_inds = randperm(numel(x_draw));
                            draw_inds = draw_inds(1:numel(amp_inds));
                            drawVal(1,i,cond,side,sub) = median(x_draw(draw_inds));
                            drawDur(1,i,cond,side,sub) = (sum(seg_dist(draw_inds))/tottime)*100;
                            
                            amp_inds = find(amp_segs>=thr(2));
                            ampVal(2,i,cond,side,sub) = median(amp_segs(amp_inds));
                            ampDurfrq(frqi,i,2) = (sum(seg_dist(amp_inds))/tottime)*100;
                            ampDur(2,i,cond,side,sub) = ampDurfrq(frqi,i,2);
                            x_draw = amp_dist(i,:);
                            draw_inds = randperm(numel(x_draw));
                            draw_inds = draw_inds(1:numel(amp_inds));
                            drawVal(2,i,cond,side,sub) = median(x_draw(draw_inds));
                            drawDur(2,i,cond,side,sub) = (sum(seg_dist(1,draw_inds))/tottime)*100;
                            
                            
                            % Relative Phi - Amplifying
                            phiBin = linspace(-pi,pi,QX);
                            phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
                            [Amp Ampbindata Ampprc binNumel] = binAmp(rp_dist,amp_dist(i,:),seg_dist,phiBin,thr,tottime,median(ampi)); % Find fraction of recording amp/deamp across bins
                            
                            % Recentre
                            for dm = 1:2 % over amp/deamp
                                A1 = Amp(:,dm);
                                [p a] = max(Ampprc);
                                mid = median(1:size(A1,1));
                                Amp_cs(:,dm) = circshift(A1,mid-a);
                            end
                            Amp_g(:,:,i,cond,side,sub) = Amp_cs;
%                             if cond ==1
                                A1 = Ampprc;
                                [p a] = max(Ampprc);
                                mid = median(1:size(A1,2));
                                xs = mid-a;
%                             end
                            Ampprc = circshift(Ampprc,xs);
                            
                            maxAmpFrq(frqi,i) = max(Ampprc);
                            if max(Ampprc)<=20; %vc_clean.PA.maxplv>=0.40
                                Ampprc_g(:,i,1,cond,side,sub) = Ampprc;
                            elseif max(Ampprc)>=50; %vc_clean.PA.maxplv<=0.15
                                Ampprc_g(:,i,2,cond,side,sub) = Ampprc;
                            end
                            for j=1:20
                                rpSur = wrapToPi(2*pi.*rand(1,length(rp_dist)));
                                [Ampz(:,:,j) Ampbindata ampprcz(:,j)] = binAmp(rpSur,amp_dist(i,:),seg_dist,phiBin,thr,tottime,median(ampi));
                            end
                            Ampz_pr(:,1,i,cond,side,sub) = prctile(squeeze(Ampz(:,1,:)),75,2);
                            Ampz_pr(:,2,i,cond,side,sub) = prctile(squeeze(Ampz(:,1,:)),5,2);
                            
                            Ampprcz_g(:,:,i,cond,side,sub) = ampprcz;
                            cmap = linspecer(2);
%                             if Lmark<3
%                                 figure(cond)
%                                 subplot(2,2,i)
%                                 b = bar(phiBinMid,Amp); grid on; hold on
%                                 b(1).FaceColor = cmap(cond,:); b(2).FaceColor = cmap(cond,:).*[0.75 1.25 1.25];
%                                 plot(phiBinMid,prctile(squeeze(Ampz(:,1,:)),75,2),'k--')
%                                 %                             plot(phiBinMid,prctile(squeeze(Ampz(:,1,:)),5,2),'k--')
%                                 xlabel('Relative Phase'); ylabel('% Occupancy'); legend('Deamp','Amp');title(obs{i})
%                                 set(gcf,'Position',[169         336        1358         766])
%                             end
                        end
                        
                        
                    end % End of Frqvar
                     frqlist = 4:1:34;
                      figure(100+cond)
                    subplot(3,1,1)
                    plot(frqlist,maxAmpFrq(:,3:4))
                    xlabel('Coherent Frequency'); ylabel('Maximum Amplification')
                    legend(obs{3:4})
                   maxAmpFrq_g(:,:,cond,side,sub) = maxAmpFrq(:,3:4);
                    
                    subplot(3,1,2)
                    ampDurFrq = ampDurfrq(:,3:4,2);
                    plot(frqlist,ampDurFrq)
                    xlabel('Coherent Frequency'); ylabel('Duration of Amplifying Segments')
                    legend(obs{3:4})
                    ampDurfrq_g(:,:,cond,side,sub) = ampDurFrq;
                    
                    subplot(3,1,3)
                    ampRatFrq = ampDurfrq(:,3:4,2)./ampDurfrq(:,3:4,1);
                    plot(frqlist,ampRatFrq)
                    xlabel('Coherent Frequency'); ylabel('Amp:Deamp Ratio')
                    legend(obs{3:4})
                    
                    ampRatfrq_g(:,:,cond,side,sub) = ampRatFrq;
                    
                    syncdur(cond,side,sub) = sum(seg_dist)/tottime;
%                     maxplv(cond,side,sub) = vc_clean.PA.maxplv;
                    
                    indz = [2 4];
                    for i=1:2
                        p = [3 2];
                        maxpow(i,cond,side,sub) = vc_clean.specanaly.powstats.maxpow(p(i));
                        maxSegPow(i,cond,side,sub) = ampVal(1,indz(i),cond,side,sub);
                    end
                    close all
                end
            end
            
%             for cond =1:2
%                 %                 if Lmark<3
%                 figure(cond+5)
%                 cmap = linspecer(2);
%                 for i =1:4
%                     subplot(2,2,i)
%                     b = bar(phiBinMid,squeeze(Ampprc_g(:,i,2,cond,side,sub))); grid on; hold on
%                     b(1).FaceColor = cmap(cond,:); %b(2).FaceColor = cmap(2,:);
%                     plot(phiBinMid,prctile(squeeze(Ampprcz_g(:,:,i,cond,side,sub)),75,2),'k--')
%                     plot(phiBinMid,prctile(squeeze(Ampprcz_g(:,:,i,cond,side,sub)),25,2),'k--')
%                     xlabel('Relative Phase'); ylabel('% Amplification');
%                     legend(b,R.condname{cond});title(obs{i});    ylim([-50 50])
%                 end
%                 set(gcf,'Position',[183   349   844   766])
%                 %                 end
%             end
%             
%             close all
        end
    end
    % Plot Spectral
    for cond = 1:2
        for i =1:2
            x = reshape((squeeze(maxAmpFrq_g(:,i,cond,:,:))),size(maxAmpFrq_g,1),[]);
            subplot(3,1,1)
            plot(frqlist,mean(x,2)); hold on
            xlabel('Coherent Frequency'); ylabel('Maximum Amplification')
            legend({'SMA \beta_1 OFF','STN \beta_1 OFF','SMA \beta_1 ON','STN \beta_1 ON'})
            
            x = reshape((squeeze(ampDurfrq_g(:,i,cond,:,:))),size(ampDurfrq_g,1),[]);
            subplot(3,1,2)
            plot(frqlist,mean(x,2)); hold on
            xlabel('Coherent Frequency'); ylabel('Duration of Amplifying Segments')
            legend({'SMA \beta_1 OFF','STN \beta_1 OFF','SMA \beta_1 ON','STN \beta_1 ON'})
            
            x = reshape((squeeze(ampRatfrq_g(:,i,cond,:,:))),size(ampRatfrq_g,1),[]);
            subplot(3,1,3)
            plot(frqlist,mean(x,2)); hold on
            xlabel('Coherent Frequency'); ylabel('Amp:Deamp Ratio')
            legend({'SMA \beta_1 OFF','STN \beta_1 OFF','SMA \beta_1 ON','STN \beta_1 ON'})
            
        end
    end
    
    
    for cond = 1:2
        for i = 1:4
            gdmpDur = squeeze(ampDur(1,i,cond,:,:));
            gddrawDur = squeeze(drawDur(1,i,cond,:,:));
            [gdmpDur,gddrawDur] = delzero(gdmpDur(:),gddrawDur(:));
            
            gampDur = squeeze(ampDur(2,i,cond,:,:));
            gdrawDur = squeeze(drawDur(2,i,cond,:,:));
            [gampDur,gdrawDur] = delzero(gampDur(:),gdrawDur(:));
            
            dat = [gdmpDur'; gddrawDur'; gampDur'; gdrawDur']';
            
            if  any(intersect([1 3],i))
                ylimz = [0 .25];
            else
                ylimz =[0 0.02];
            end
            figure(10+cond)
            subplot(2,2,i)
            bploter(dat,{'Power < 25%','Rnd Draw 1','Power > 75%','Rnd Draw 2'},ylimz)
            ylabel('Average Length');
            title(obs{i})
            %             ylim([ylimz(1) ylimz(2)*1.1]);
            
            
            figure(20+cond)
            subplot(2,2,i)
            x = squeeze(Amp_g(:,2,i,cond,:,:));
            x1 = reshape(x,size(x,1),size(x,2),[]);
            x1 = x1(:,sum(x1)~=0);
            boxplot(x1','labels', sprintfc('%d',1:QX-1),... %
                'BoxStyle','filled','Widths',0.8);
            hold on
            for dm =1:2
                x = squeeze(Ampz_pr(:,dm,i,cond,:,:));
                x1 = reshape(x,size(x,1),[]);
                x1 = x1(:,sum(x1)~=0);
                plot(mean(x1'),'k--')
            end
            ax = gca;
            ax.XTickLabel = round(phiBinMid,2);
            xlabel('Relative Phase (rads)');
            ylabel('% Occupancy')
            title(obs{i})
            grid on; ylim([0 3.5])
            
            for d = 1:2
                figure(30+d)
                subplot(2,2,i)
                x = squeeze(Ampprc_g(:,i,d,1,:,:));
                x = reshape(x,size(Ampprc_g,1),[]);
                x = x(:,sum(x,1)~=0)
                x1 = nanmean(x');
                %             bar(phiBinMid,nanmean(x'),'b')
                plot(phiBinMid,x,'b:');hold on
                plot(phiBinMid,x1,'b','LineWidth',3)
                
                
                x = squeeze(Ampprc_g(:,i,d,2,:,:));
                x = reshape(x,size(Ampprc_g,1),[]);
                x = x(:,sum(x,1)~=0);
                x2 = nanmean(x');
                plot(phiBinMid,x,'r:')
                plot(phiBinMid,x2,'r','LineWidth',3)
                
                xlabel('Relative Phase (rads)');
                ylabel('% Amplification')
                title(obs{i})
                grid on; ylim([-50 100])
%                 plot(phiBinMid,mean(squeeze(Ampprcz_g(:,1,i,cond,1,:)),2),'k--')
%                 plot(phiBinMid,mean(squeeze(Ampprcz_g(:,2,i,cond,1,:))),'k--')
%                 
                
                set(gcf,'Position',[702   332   842   766])
                
            end
            %             bar(phiBinMid,[x1; x2]')
        end
        figure(10+cond)
        set(gcf,'Position',[173         332        1371         766])
        figure(20+cond)
        set(gcf,'Position',[173         332        1371         766])
    end
    % Plot of sync duration vs maxplv
    [maxplv,syncduri] = delzero(maxplv(:),syncdur(:));
    figure
    linplot(syncduri,maxplv,'Sync Rate',[R.bandinits{3} 'STN/SMA PLV'],'ro');
    
    for i =1:2
        p = [3 2];
        x = maxpow(i,:,:,:);
        y = maxSegPow(i,:,:,:);
        [x,y] = delzero(x,y);
        figure
        linplot(x',y','Maximum Power in Band',[R.bandinits{p(i)} ' Frame Power'],'ro');
        if i == 1
            %             ylim([0 4]); xlim([0 0.08]);
        elseif i ==2
            ylim([0 4]); xlim([0 0.08]);
        end
    end
    
end

function bploter(dat,gtit,ylimz)
boxplot(dat,'labels',gtit,...
    'BoxStyle','filled','Widths',0.8); %,'datalim',ylimz)
%%
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',35); % Set width
idx=strcmpi(t,'Whisker');  % Find Box objects
whisker=a(idx);          % Get the children you need
set(whisker,'LineWidth',2); % Set width
idx=strcmpi(t,'Median');  % Find Box objects
whisker=a(idx);          % Get the children you need
set(whisker,'LineWidth',2); % Set width
%%
pv(1) = ranksum(dat(:,1),dat(:,2));
pv(2) = ranksum(dat(:,3),dat(:,4));
grid on; H=sigstar({gtit(1:2) gtit(3:4)},pv);
a = get(gca,'Children');
a = findobj(a,'type','Text');
set(a,'FontSize',16)

function [x,y] = delzero(x,y)
y(x==0) = []; x(x==0) = [];
x(y==0) = []; y(y==0) = [];


