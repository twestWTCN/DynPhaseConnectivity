function [phi_dist amp_dist segL_ddt] = compute_dynPhaseLocking(R,Xdata,band,stn_lb_frq)

% Get the overal signal amplitudes
sig =  Xdata.trial{1};
signalEnvAmp = median(abs(hilbert(sig)),2);

% Compute optimal frequencies based on PLV
[maxfrq,maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band,stn_lb_frq); %
[SRPeps Ampeps SNR_eps PLVeps] = phase_amp_surrComp(R,Xdata,band,maxfrq,stn_lb_frq,R.PA.AmpSurrN);
SNR_eps_z(1) = 10*log10(SNR_eps(1)./signalEnvAmp(1));
SNR_eps_z(2) =  10*log10(SNR_eps(2)./signalEnvAmp(2));
SNR_eps_z(3) =  10*log10(SNR_eps(3)./signalEnvAmp(2));
% Compute data transforms (Hilbert)
[amp phi dphi_12 dphi_12_dt Xdata] = comp_instant_angle_phase(Xdata,maxfrq,stn_lb_frq,R.PA.bwid,band);

% Sliding Window PLV
if R.PA.SType == 1
    fsamp = Xdata.fsample;
    %                         swsize = R.PA.slidingwindow;
    swsize = floor((10/maxfrq)*fsamp);
    WinSize = swsize;
    [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver,R.PA.optimalPLFrqMeth);
    PLV_tvec = Xdata.time{1}(PLV_tvec);
    amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
    clear SNR_sw
    SNR_sw(:,1) = 10.*log10(amp_sw(:,1)./signalEnvAmp(1));
    SNR_sw(:,2) =  10.*log10(amp_sw(:,2)./signalEnvAmp(2));
    SNR_sw(:,3) =  10.*log10(amp_sw(:,3)./signalEnvAmp(2));
    %                                             figure(1)
    %                                             SNR_Inspector(R,PLV_tvec,SNR_sw,SNR_eps_z,band)
    dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,round(R.PA.WinOver*WinSize));
    SW_sampr = max(diff(PLV_tvec));
    tseries = dphi_12_sw; qstable = find(PLV>PLVeps);
    mwid = R.PA.mwid;
    period = (mwid/maxfrq)/SW_sampr;
    % Minimum number of cycles to consider sync
    [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
    
    %                         any(amp_dist(:)>300)
                            plot_example_phaseanalysis_SW(Xdata,amp,phi,PLV,seg_ddt1,PLVeps,PLV_tvec);
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
    figure(2)
    plot_example_phaseanalysis_trace(Xdata,amp,phi,dphi_12_dt,seg_ddt1,SRPeps,fsamp);
    clf(1); clf(2)
end
