function [SRPeps,Ampeps,SNReps,PLVeps] = phase_amp_surrComp(R,Xdata,band,frq,stn_lb_frq,NR)
fsamp = Xdata.fsample;
WinSize = R.PA.slidingwindow*fsamp;


surdata = Xdata.trial{1};
surfft(:,1) = fft(surdata(1,:));
surfft(:,2) = fft(surdata(2,:));
Amp = abs(surfft)';
Phi = angle(surfft)';
% progressStepSize = 5;
% gcp;
% ppm = ParforProgMon('Power Estimation: ', NR, progressStepSize, 800, 300);
parfor N = 1:NR
    % Phase Shuffling
    phasescm =[];
    Phishuff = Phi(1,randperm(length(Phi(1,:))));
    phasescm(1,:) = abs(ifft(Amp(1,:).*exp(sqrt(-1)*Phishuff)));
    Phishuff = Phi(2,randperm(length(Phi(2,:))));
    phasescm(2,:) = abs(ifft(Amp(2,:).*exp(sqrt(-1)*Phishuff)));
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = phasescm;
    Rdata.time = Xdata.time;
    [amp,phi,dphi_12,dphi_12_dt,betaS] = comp_instant_angle_phase(Rdata,frq,stn_lb_frq,R.PA.bwid,band);
    
    SRPbank(:,N) = mean(abs(dphi_12_dt));
    PLV = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
    PLVbank(:,N) = mean(PLV);
    % Amp Shuffling
    ampscm =[];
    %     Ampshuff = Amp(1,randperm(length(Phi(1,:))));
    %     ampscm(1,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(1,:))));
    %     Ampshuff = Amp(2,randperm(length(Phi(2,:))));
    %     ampscm(2,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(2,:))));
    ampscm(1,:) =surdata(1,randperm(length(surdata(1,:))));
    ampscm(2,:) =surdata(2,randperm(length(surdata(2,:))));
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = ampscm;
    Rdata.time = Xdata.time;
    [amp,phi,dphi_12,dphi_12_dt,betaS] = comp_instant_angle_phase(Rdata,frq,stn_lb_frq,R.PA.bwid,band);
    ampbank(:,:,N) = amp;
    %     ppm.increment();
end
% ppm.delete();
PLVeps = prctile(PLVbank(:),R.PA.PLVeps_prctile);
SRPeps = prctile(SRPbank(:),R.PA.SRPeps_prctile);
Ampeps = prctile(reshape(permute(ampbank,[1 3 2]),[],3),50,1);
SNReps = prctile(reshape(permute(ampbank,[1 3 2]),[],3),R.PA.SNReps_prctile,1);