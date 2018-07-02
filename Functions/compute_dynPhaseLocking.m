function [phi_dist amp_dist segL_ddt ampT maxPLV  PLV_tvec Ampeps surr OVL] = compute_dynPhaseLocking(R,Xdata,band,cohfrq,powfrq,optfrq)
if nargin<6
    optfrq = 1;
end
% Get the overal signal amplitudes
sig =  Xdata.trial{1};
signalEnvAmp = median(abs(hilbert(sig)),2);

% Compute optimal frequencies based on PLV
if optfrq ==1
    [cohfrq(band),maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band); %
else
    maxPLV = 1;
end
% [Ampeps,PLVeps] = phase_amp_surrComp(R,Xdata,band,cohfrq,powfrq,R.PA.AmpSurrN);
% disp(PLVeps)
% disp(Ampeps)

PLVeps = [0.75   0.70    0.70];
Ampeps = [0 0   24
          0 0   24];
for band = 1:size(R.bandef)
    SNR_eps_z(:,band) = 10*log10(Ampeps(:,band)./signalEnvAmp);
end
% Compute data transforms (Hilbert)
[amp phi dphi_12 dphi_12_dt XdataT] = comp_instant_angle_phase(Xdata,cohfrq(band),R.PA.bwid,band);
[lb_amp] = comp_instant_angle_phase(Xdata,powfrq(2),R.PA.bwid,2);

% Sliding Window PLV
if R.PA.SType == 1
    fsamp = XdataT.fsample;
    swsize = R.PA.slidingwindow;
%     swsize = floor((10/cohfrq(band))*fsamp);
    WinSize = swsize*fsamp;
    [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
    PLV_tvec = XdataT.time{1}(PLV_tvec);
    amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
    lb_amp_sw = cont2slidingwindow(lb_amp,WinSize,floor(R.PA.WinOver*WinSize));
    
    clear SNR_sw
    SNR_sw(:,1:2) = 10.*log10(amp_sw./signalEnvAmp');
    SNR_sw(:,3:4) = 10.*log10(lb_amp_sw./signalEnvAmp');
%     figure(1)
%     SNR_Inspector(R,PLV_tvec,SNR_sw,SNR_eps_z,band,{'SMA','STN','SMA','STN'})
    dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,floor(R.PA.WinOver*WinSize));
    SW_sampr = max(diff(PLV_tvec));
    PLVeps(band) = prctile(PLV,75);
    tseries = dphi_12_sw; qstable = find(PLV>PLVeps(band));
    mwid = 12; %R.PA.mwid;
    period = (mwid/cohfrq(band))/SW_sampr;
    % Minimum number of cycles to consider sync
    ampT = [amp lb_amp];
%     ampT = 10.*log10(ampT./[signalEnvAmp;signalEnvAmp]');

    [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,ampT,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
    amp_lb_eps = prctile(lb_amp_sw(:,2),75);
    qstable = find(lb_amp_sw(:,2)>amp_lb_eps);
    [phi_lb_dist amp_lb_dist seg_lb_ddt1 segL_lb_ddt consecSegs_lb] = analysestablesegs(qstable',tseries,ampT,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
    
    OVL = ((numel(intersect(seg_lb_ddt1,seg_ddt1)).*SW_sampr)./Xdata.time{1}([end]))*100;
%         surr = phase_amp_surrCompPost(R,Xdata,band,cohfrq(band),powfrq(2),PLVeps(band),100);
    surr= [];
%     any(amp_dist(:)>300)
%     plot_example_phaseanalysis_SW(XdataT,amp,phi,PLV,seg_ddt1,PLVeps(band),lb_amp_sw(:,2),seg_lb_ddt1,amp_lb_eps,PLV_tvec);
elseif R.PA.SType == 2 % Sliding Window PhaseAng. Stability
    clear SNR_sw
    SNR_sw(:,1) = 10.*log10(amp(:,1)./signalEnvAmp(1));
    SNR_sw(:,2) =  10.*log10(amp(:,2)./signalEnvAmp(2));
    SNR_sw(:,3) =  10.*log10(amp(:,3)./signalEnvAmp(2));
    % Plot SNR
    %                         figure(1)
    %                         SNR_Inspector(R,Xdata.time{1},SNR_sw,SNR_eps_z,breg)
    
    %%% Set Length Constraints
    fsamp = Xdata.fsample;
    mwid = R.PA.mwid; % Minimum number of cycles to consider sync
    period = (mwid/maxfrq)*fsamp;
    
    %%% Find segments
    %                 dphi_12_dt_sm = smooth(dphi_12_dt,period/6)';
    qstable = find(abs(dphi_12_dt')<SRPeps); % The points that are below threshold
    tseries = wrapToPi(dphi_12(2:end));
    [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
%     figure(2)
%     plot_example_phaseanalysis_trace(Xdata,amp,phi,dphi_12_dt,seg_ddt1,SRPeps,fsamp);
%     clf(1); clf(2)
end
