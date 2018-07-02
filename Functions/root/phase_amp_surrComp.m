function [Ampeps,PLVeps] = phase_amp_surrComp(R,Xdata,band,cohfrq,powfrq,NR)
fsamp = Xdata.fsample;
WinSize = R.PA.slidingwindow*fsamp;
data = Xdata.trial{1};
% progressStepSize = 5;
% gcp;
% ppm = ParforProgMon('Power Estimation: ', NR, progressStepSize, 800, 300);
parfor N = 1:NR
    % Phase Shuffling
    
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = phaseran(data',1)';
    Rdata.time = Xdata.time;
    [~,phi] = comp_instant_angle_phase(Rdata,cohfrq(band),[2 2 2],band);
    PLV = slidingwindowPLV(WinSize,phi,0,R.PA.optimalPLFrqMeth);
    PLVbank(:,N) = PLV(:);
    % Amp Shuffling
    ampscm =[];
    %     Ampshuff = Amp(1,randperm(length(Phi(1,:))));
    %     ampscm(1,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(1,:))));
    %     Ampshuff = Amp(2,randperm(length(Phi(2,:))));
    %     ampscm(2,:) = abs(ifft(Ampshuff.*exp(sqrt(-1)*Phi(2,:))));
    ampscm(1,:) =data(1,randperm(length(data(1,:))));
    ampscm(2,:) =data(2,randperm(length(data(2,:))));
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = ampscm;
    Rdata.time = Xdata.time;
    [amp,~] = comp_instant_angle_phase(Rdata,cohfrq(band),R.PA.bwid,band);
    ampbank(:,:,N) = amp;
    %     ppm.increment();
end
PLVeps(:,band) = prctile(PLVbank(:),R.PA.PLVeps_prctile);
Ampeps(:,band) = prctile(reshape(permute(ampbank,[1 3 2]),[],2),50,1);
% SNReps = prctile(reshape(permute(ampbank,[1 3 2]),[],2),R.PA.SNReps_prctile,1);
% ppm.delete();
