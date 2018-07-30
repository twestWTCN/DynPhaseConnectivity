function [PLV_Phi_Amp,Amp_vs_PLV_RP,OVL] = compute_dynPhaseLocking_v2(R,Xdata,band,altband,cohfrq,powfrq,optfrq)
if nargin<6
    optfrq = 1;
end
% Compute optimal frequencies based on PLV
if optfrq ==1
    [cohfrq(band),maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band); %
else
    maxPLV = 1;
end

% Compute data transforms (Hilbert)
[amp phi dphi_12 dphi_12_dt XdataT] = comp_instant_angle_phase(Xdata,cohfrq(band),R.PA.bwid(band));
alt_amp = comp_instant_angle_phase(Xdata,powfrq(2),R.PA.bwid(altband));
% Make Amps adjustment over median
[PLV_Phi_Amp,Amp_vs_PLV_RP,OVL] = SW_PLV_Analysis_v1(R,amp,phi,dphi_12,XdataT,alt_amp,cohfrq(band));
 