function surr = phase_amp_surrCompPost(R,Xdata,band,cohfrq,powfrq,PLVeps,NR)
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
    
    [amp phi dphi_12 dphi_12_dt XdataT] = comp_instant_angle_phase(Rdata,cohfrq,R.PA.bwid,band);
    [lb_amp] = comp_instant_angle_phase(Rdata,powfrq,R.PA.bwid,2);
    
    [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
    PLV_tvec = XdataT.time{1}(PLV_tvec);
    PLV_tvec_store{N} = PLV_tvec(end)-PLV_tvec(1);
    amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
    lb_amp_sw = cont2slidingwindow(lb_amp,WinSize,floor(R.PA.WinOver*WinSize));
    
    dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,round(R.PA.WinOver*WinSize));
    SW_sampr = max(diff(PLV_tvec));
    tseries = dphi_12_sw; qstable = find(PLV>PLVeps);
    mwid = R.PA.mwid;
    period = (mwid/cohfrq)/SW_sampr;
    % Minimum number of cycles to consider sync
    ampT = [amp lb_amp];
%     ampT = 10.*log10(ampT./[signalEnvAmp;signalEnvAmp]');
    [phi_dist{N} amp_dist{N} seg_ddt1 segL_ddt{N} consecSegs H] = analysestablesegs(qstable,tseries,ampT,period,mwid,fsamp);
end

surr.pA_pli_dist_save = vertcat(phi_dist{:});
surr.amp_pli_dist_save = horzcat(amp_dist{:});
surr.segL_pli_dist_save = [segL_ddt{:}];
surr.PA_time = [PLV_tvec_store{:}];


