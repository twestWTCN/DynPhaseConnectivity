function [PLV_Phi_Amp,Amp_vs_PLV_RP,OVL] = SW_PLV_Analysis_v1(R,amp,phi,dphi_12,XdataT,alt_amp,cohfrq)
% Sliding Window PLV
    fsamp = XdataT.fsample;
    AbsAmp = median(abs(hilbert(XdataT.trial{1}))');
    AbsAmp = [AbsAmp AbsAmp];
    SurrAmp = [R.PA.Ampeps R.PA.Ampeps];
    swsize = R.PA.slidingwindow;
    WinSize = swsize*fsamp;
    % Compute Main Metric
    [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
    PLV_tvec = XdataT.time{1}(PLV_tvec);
    
    % Convert timeseries to SW scale
    amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
    alt_amp_sw = cont2slidingwindow(alt_amp,WinSize,floor(R.PA.WinOver*WinSize));
    dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,floor(R.PA.WinOver*WinSize));
    SW_sampr = 1/max(diff(PLV_tvec));

    % Minimum number of cycles to consider sync
    mwid = R.PA.mwid;
    period = (mwid/cohfrq).*SW_sampr;
    % Amplitude Trace
    ampT = [amp_sw alt_amp_sw];
    ampT = ((ampT-median(ampT))./median(ampT))*100;
%     ampT = ((ampT-AbsAmp)./AbsAmp)*100;
    % Analyse Segments
    % Relative phase vs Amp
    RP = dphi_12_sw; PLV_stable = find(PLV>R.PA.PLVeps);
    [PLV_dist PLV_phi_dist PLV_amp_dist PLV_segL PLV_consecSegs] = analysestablesegs_PLV_vs_Phi_Amp(PLV,PLV_stable,RP,ampT,period,SW_sampr);
    
    PLV_Phi_Amp.PLV_dist = PLV_dist;
    PLV_Phi_Amp.PLV_phi_dist = PLV_phi_dist;
    PLV_Phi_Amp.PLV_amp_dist = PLV_amp_dist;
    PLV_Phi_Amp.PLV_segL = PLV_segL;
    PLV_Phi_Amp.PLV_consecSegs = PLV_consecSegs;
    
    % Supra-Low Beta Amp vs Relative Phase
%     amp_alt_eps = prctile(ampT(:,4),20);
    AMP_supra = find(ampT(:,4)>15); % Greater than 20% of median
    [AMP_dist AMP_phi_dist AMP_plv_dist AMP_segL AMP_consecSegs] =  analysestablesegs_Amp_vs_PLV_RP(AMP_supra',RP,ampT,period,SW_sampr);
    
    Amp_vs_PLV_RP.AMP_dist = AMP_dist;
    Amp_vs_PLV_RP.AMP_phi_dist = AMP_phi_dist;
    Amp_vs_PLV_RP.AMP_plv_dist = AMP_plv_dist;
    Amp_vs_PLV_RP.AMP_segL = AMP_segL;
    
    OVL = ((numel(intersect([PLV_consecSegs{:}],[AMP_consecSegs{:}]))./SW_sampr)./XdataT.time{1}([end]))*100;
