function [Ampeps,PLVeps] = phase_amp_surrComp(R,Xdata,band,cohfrq,powfrq)
fsamp = Xdata.fsample;
WinSize = R.PA.slidingwindow*fsamp;
data = Xdata.trial{1};
% progressStepSize = 5;
% gcp;
% ppm = ParforProgMon('Power Estimation: ', NR, progressStepSize, 800, 300);
parfor N = 1:R.PA.AmpSurrN
    % Phase Shuffling
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = phaseran(data',1)';
    Rdata.time = Xdata.time;
    [amp phi dphi_12 dphi_12_dt] = comp_instant_angle_phase(Rdata,cohfrq,R.PA.bwid(band));
    
    %     rm = rem(numel(dphi_12_dt),100);
    %     rnblock = reshape(dphi_12_dt(1:end-rm),100,[]);
    %     rnblock = rnblock(:,randperm(size(rnblock,2)));
    %     dphi_12_dt = [rnblock(:); dphi_12_dt(end-rm:end)];
    %     dphi_12 = cumsum(dphi_12_dt);
    %     phi = [dphi_12*1 dphi_12*2];
    %     phi = wrapToPi(phi);
    PLV = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
    %     PLVbank(:,N) = PLV(:);
    PLVepsN(N) = prctile(PLV,R.PA.PLVeps_prctile);
    % Amp Shuffling
    ampscm =[];
    ampscm(1,:) =data(1,randperm(length(data(1,:))));
    ampscm(2,:) =data(2,randperm(length(data(2,:))));
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = ampscm;
    Rdata.time = Xdata.time;
    [amp,~] = comp_instant_angle_phase(Rdata,powfrq(1),R.PA.bwid(band));
    ampbank(:,:,N) = amp;
    %     ppm.increment();
end
PLVeps = median(PLVepsN(:));
Ampeps = prctile(reshape(permute(ampbank,[1 3 2]),[],2),50,1);
% SNReps = prctile(reshape(permute(ampbank,[1 3 2]),[],2),R.PA.SNReps_prctile,1);
% ppm.delete();



%     Ampshuff = Amp(1,randperm(length(Phi(1,:))));
%     ampscm(1,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(1,:))));
%     Ampshuff = Amp(2,randperm(length(Phi(2,:))));
%     ampscm(2,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(2,:))));
