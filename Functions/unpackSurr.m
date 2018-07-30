function [RP Amp SegL TotTimesurr OVLsurr] = unpackSurr(surr)
for i = 1:numel(surr.PLV_Phi_Amp)
    RP{i} = surr.PLV_Phi_Amp(i).PLV_phi_dist';
    Amp{i} = surr.PLV_Phi_Amp(i).PLV_amp_dist;
    SegL{i} = surr.PLV_Phi_Amp(i).PLV_segL;
end
OVLsurr = surr.OVL;
TotTimesurr = sum(surr.TotTime)
RP = [RP{:}];
Amp = [Amp{:}];
SegL = [SegL{:}];