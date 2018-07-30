function [PLV_Phi_Amp,Amp_vs_PLV_RP,OVL,TotTime] = phase_amp_surrComp_SWPLV(R,Xdata,band,cohfrq,powfrq)
fsamp = Xdata.fsample;
WinSize = R.PA.slidingwindow*fsamp;
data = Xdata.trial{1};
% progressStepSize = 5;
% gcp;
% ppm = ParforProgMon('Power Estimation: ', NR, progressStepSize, 800, 300);
parfor N = 1:R.PA.AmpSurrN
    % Phase Shuffling
    ampscm =[];
    ampscm(1,:) =data(1,randperm(length(data(1,:))));
    ampscm(2,:) =data(2,randperm(length(data(2,:))));
    Rdata = [];
    Rdata.label = Xdata.label;
    Rdata.fsample = Xdata.fsample;
    Rdata.trial{1} = ampscm; %phaseran(data',1)';
    Rdata.time = Xdata.time;   
    [amp phi dphi_12 dphi_12_dt] = comp_instant_angle_phase(Rdata,cohfrq,R.PA.bwid(band));
    TotTime(N) = Rdata.time{1}(end)
    % Randomize
%     rm = rem(numel(dphi_12_dt),100);
%     rnblock = reshape(dphi_12_dt(1:end-rm),100,[]);
%     rnblock = rnblock(:,randperm(size(rnblock,2)));
%     dphi_12_dt = [rnblock(:); dphi_12_dt(end-rm:end)];
%     dphi_12 = cumsum(dphi_12_dt);
%     phi = [dphi_12*1 dphi_12*2];
%     phi = wrapToPi(phi);
%     PLV = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);

    % Amp Shuffling
    ampscm =[];
    ampscm(1,:) =data(1,randperm(length(data(1,:))));
    ampscm(2,:) =data(2,randperm(length(data(2,:))));
    Adata = [];
    Adata.label = Xdata.label;
    Adata.fsample = Xdata.fsample;
    Adata.trial{1} = ampscm;
    Adata.time = Xdata.time;
    [alt_amp,~] = comp_instant_angle_phase(Adata,powfrq(1),R.PA.bwid(band));
    [PLV_Phi_Amp(N),Amp_vs_PLV_RP(N),OVL(N)] = SW_PLV_Analysis_v1(R,amp,phi,dphi_12,Rdata,alt_amp,cohfrq) ;   
    

    %     ppm.increment();
end
