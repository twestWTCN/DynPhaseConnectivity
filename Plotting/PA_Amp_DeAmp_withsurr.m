function dwell = PA_Amp_DeAmp(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
close all
QX = 6; Lmark = 0;
for breg = 2:length(R.bregname)
    for sub = 3:length(R.subname)
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences
                    Lmark = Lmark+1;
                    
                    obs = {[R.bregname{breg} ' ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{R.bregband{breg}}],['SMA ' R.bandinits{2}],['STN ' R.bandinits{2}]};
                    for i = 1:4
                        titlist{i} = ['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{i} ' Power'];
                    end
                    tottime = vc_clean.PA.timevec(end)-vc_clean.PA.timevec(1);
                    
                    pureAmp = vc_clean.PA.pureAmp';
                    amp_dist = vc_clean.PA.amp_pli_dist_save;
                    seg_dist = vc_clean.PA.segL_pli_dist_save;
                    rp_dist = vc_clean.PA.pA_pli_dist_save;
                    
                    
                    amp_dist_sur = vc_clean.PA.surr.amp_pli_dist_save;
                    seg_dist_sur = vc_clean.PA.surr.segL_pli_dist_save;
                    rp_dist_sur = vc_clean.PA.surr.pA_pli_dist_save;
                    tottime_sur =  sum(vc_clean.PA.surr.PA_time);
                    % First Establish Relationship between PL and Amplitude
                    for i = 1:4
                        ampi = pureAmp(i,:); % Find the amplitude envelope
                        amp_segs = amp_dist(i,:); % Find amplitude within segs
                        amp_segs_ss = amp_dist_sur(i,:);
                        thr = [prctile(ampi,25) prctile(ampi,75)]; % Treat at deamp or amp
                        %                         subplot(2,2,i)
                        %                         hist(ampi,-100:25:200)
                        amp_inds = find(amp_segs<=thr(1));
                        ampVal(1,i,cond,side,sub) = median(amp_segs(amp_inds));
                        ampDur(1,i,cond,side,sub) = (median(seg_dist(amp_inds))/tottime)*100;
                        
%                         x_draw = amp_dist(i,:);
%                         draw_inds = randperm(numel(x_draw));
%                         draw_inds = draw_inds(1:numel(amp_inds));
                        ss_inds =  find(amp_segs_ss<=thr(1));
                        drawVal(1,i,cond,side,sub) = median(amp_segs_ss(ss_inds));
                        drawDur(1,i,cond,side,sub) = (median(seg_dist_sur(ss_inds))/tottime)*100;
                        
                        amp_inds = find(amp_segs>=thr(2));
                        ampVal(2,i,cond,side,sub) = median(amp_segs(amp_inds));
                        ampDur(2,i,cond,side,sub) = (median(seg_dist(amp_inds))/tottime)*100;
                        
%                         x_draw = amp_dist(i,:);
%                         draw_inds = randperm(numel(x_draw));
%                         draw_inds = draw_inds (1:numel(amp_inds));

                        ss_inds =  find(amp_segs_ss>=thr(2));
                        drawVal(2,i,cond,side,sub) = median(amp_segs_ss(ss_inds));
                        drawDur(2,i,cond,side,sub) = (median(seg_dist_sur(1,ss_inds))/tottime)*100;
                        
                        
                        % Relative Phi - Amplifying
                        phiBin = linspace(-pi,pi,QX);
                        phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
                        [Amp] = binAmp(rp_dist,amp_dist(i,:),seg_dist,phiBin,thr,tottime); % Find fraction of recording amp/deamp across bins
                        
                        % 
                        for dm = 1:2 % over amp/deamp
                            A1 = Amp(:,dm);
                            [p a] = max(A1);
                            mid = median(1:size(A1,1));
                            Amp_cs(:,dm) = circshift(A1,mid-a);
                        end
                        Amp_g(:,:,i,cond,side,sub) = Amp_cs;
                        
                        for j=1:2000
                            rpSur = wrapToPi(2*pi.*rand(1,length(rp_dist)));
                            [Ampz(:,:,j)] = binAmp(rpSur,amp_dist(i,:),seg_dist,phiBin,thr,tottime);
                        end
                        Ampz_pr(:,1,i,cond,side,sub) = prctile(squeeze(Ampz(:,1,:)),95,2);
                        Ampz_pr(:,2,i,cond,side,sub) = prctile(squeeze(Ampz(:,1,:)),5,2);
                        
                        
%                        [Ampz] = binAmp(rp_dist_sur,amp_segs_ss,seg_dist_sur,phiBin,thr,tottime_sur);
%                         Ampz_pr(:,1,i,cond,side,sub) = Ampz(:,1);
%                         Ampz_pr(:,2,i,cond,side,sub) = Ampz(:,2);
                        
                        if Lmark<3
                            figure(cond)
                            subplot(2,2,i)
                            bar(phiBinMid,Amp); grid on; hold on
                            plot(phiBinMid,prctile(squeeze(Ampz(:,1,:)),95,2),'k--')
                            plot(phiBinMid,prctile(squeeze(Ampz(:,1,:)),5,2),'k--')
%                             plot(phiBinMid,Ampz(:,1),'k--')
%                             plot(phiBinMid,Ampz(:,2),'k--')
                            xlabel('Relative Phase'); ylabel('% Occupancy'); legend('Deamp','Amp');title(obs{i})
                            set(gcf,'Position',[169         336        1358         766])
                        end
                    end
                    
                    syncdur(cond,side,sub) = sum(seg_dist)/tottime;
                    maxplv(cond,side,sub) = vc_clean.PA.maxplv;
                    
                    indz = [2 4];
                    for i=1:2
                        p = [3 2];
                        maxpow(i,cond,side,sub) = vc_clean.specanaly.powstats.maxpow(p(i));
                        maxSegPow(i,cond,side,sub) = ampVal(1,indz(i),cond,side,sub);
                    end
                    
                end
            end
        end
    end
    for cond = 1:2
        for i = 1:4
             gdmpDur = squeeze(ampDur(1,i,cond,:,:));
            gddrawDur = squeeze(drawDur(1,i,cond,:,:));
            [gdmpDur,gddrawDur] = delzero(gdmpDur,gddrawDur);
           
            gampDur = squeeze(ampDur(2,i,cond,:,:));
            gdrawDur = squeeze(drawDur(2,i,cond,:,:));
            [gampDur,gdrawDur] = delzero(gampDur,gdrawDur);
            
            dat = [gdmpDur; gddrawDur; gampDur; gdrawDur]';
            
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
                boxplot(x1','labels',{'1','2','3','4','5'},... % ,'6','7'
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
        end
        figure(10+cond)
        set(gcf,'Position',[173         332        1371         766])
        figure(20+cond)
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


