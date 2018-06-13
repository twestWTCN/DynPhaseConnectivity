function [ampBinGroup segBinGroup phipeakGroup RcoeffGroup] =  plot_phase_amp_analysis_generic(R,ylimlistS,obs,breg,cond,side,sub,phiBin,...
    cmapint,segLCol,relativePhiCol,ampSegCol,...
    ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup,...
    panlist,titleR)


%                     ampsegColGroup{cond,side,sub} = ampSegCol;
%                     segLColGroup{cond,side,sub} = segLCol;
phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
[sub side cond]
for i = 1:3
    [shiftPhiCol phipeak(i) binind]= findAmpPhi(R,ampSegCol(i,:),relativePhiCol,phiBin);
    [ampBinMu(i,:) ampBinSEM(i,:)] = binstats(shiftPhiCol,ampSegCol(i,:),phiBin);
    try
        selind{i} = find(relativePhiCol>=phiBin(binind-1) & relativePhiCol<=phiBin(binind+1));
    catch
        disp('Could not compute bins'); selind{i} = NaN;
    end
end
[shiftPhiCol phipeak(4) bind]= findAmpPhi(R,segLCol,relativePhiCol,phiBin);
[segBinMu segBinSEM] = binstats(shiftPhiCol,segLCol,phiBin);
try
    selind{4} = find(relativePhiCol>=phiBin(binind-1) & relativePhiCol<=phiBin(binind+1));
catch
    disp('Could not compute bins'); selind{i} = NaN;
end
segBinSEM(isnan(segBinSEM)) = 0;


figure(30)
ylimlist = ylimlistS{3};
for i = 1:3
    subplot(3,2,panlist(cond,i))
    try
        [A stat] = linplot_PD(log10(segLCol)',ampSegCol(i,:)','Seg Length (s)','Amplitude',cmapint(i,:),0); xlim([-1.5 1])
        ylim(ylimlist{breg}{i})
        coef = stat.modcoef.*(stat.p<0.05); if size(coef,1)>size(coef,2); coef  =coef'; end
    catch
        coef = [0 0];
        disp('Correlation couldnt run!!!')
    end
    Rcoeff(:,:,i) = coef;
    title(titleR{1,i})
    ylabel(['% Change in ' obs{i}]); xlabel('Segment Length (s)')
    
end

RcoeffGroup{cond,side,sub} = Rcoeff;
ampBinGroup{cond,side,sub} = ampBinMu;
segBinGroup{cond,side,sub} = segBinMu;
phipeakGroup{cond,side,sub} = phipeak;

figure(1)
cmap = linspecer(3);
ylimlist = ylimlistS{1};
for i = 1:3
    subplot(4,2,panlist(cond,i))
    hl = plot(phiBinMid', ampBinMu(i,:)','--','color',cmap(i,:)); hold on
    % %                         [hl, hp] = boundedline(phiBinMid', ampBinMu(i,:)',ampBinSEM(i,:)','cmap',cmap(i,:)); hold on
    % %                         if cond == 1; hl.LineStyle = '--'; end
    % %                         hp.FaceAlpha = 0.4;
    % %                         [xq yq R2 exitflag] = VMfit(phiBinMid',ampBinMu(i,:)',20,[0,-pi,-100],[100,pi,1e3],0);
    % %                         if exitflag == 1
    % %                             hold on; plot(xq,yq,'color',cmap(i,:));
    % %                         end
    % %                         [xq yq R2 exitflag] = sinfit(phiBinMid', ampBinMu(i,:)',20,[25; 2*pi; 0; -500],[300;2*pi; 2*pi; 500],0);
    % %                         plot(xq,yq,'color',cmap(i,:));
    ylim(ylimlist{breg}{i})
    title(titleR{2,i})
    ylabel(['% Change in ' obs{i}]); xlabel('Relative Phase (rads)')
    grid on
end

cmap = linspecer(5);
panlist(:,4) = [7; 8];
subplot(4,2,panlist(cond,4))
hl = plot(phiBinMid', segBinMu(1,:)','--','color',cmap(5,:)); hold on
% %                     [hl, hp] = boundedline(phiBinMid', segBinMu(1,:)',segBinSEM(1,:)','cmap',cmap(5,:));
% %                     if cond == 1; hl.LineStyle = '--'; end
% %                     hp.FaceAlpha = 0.4;
% %                     [xq yq R2 exitflag] = VMfit(phiBinMid',segBinMu(1,:)',20,[0,-10,0.01],[1,10,1e3],0);
% %                     [xq yq R2 exitflag] = sinfit(phiBinMid', segBinMu(1,:)',20,[0.01; 2*pi; 0; -1],[0.5; 2*pi; 2*pi; 1],0);
% %                     if exitflag == 1
% %                     hold on; plot(xq,yq,'color',cmap(5,:));
% %                     end
ylimlist = ylimlistS{2};
ylim(ylimlist{breg}{1}); grid on
title(titleR{4,4})

ylabel(['LB Segment Length (s)']); xlabel('Relative Phase (rads)')
