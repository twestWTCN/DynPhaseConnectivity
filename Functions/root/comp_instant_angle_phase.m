function [amp phi dphi_12 dphi_12_dt Xdata] = comp_instant_angle_phase(Odata,frq,bwid)
fsamp = Odata.fsample;
Xdata.trial{1} = ft_preproc_bandpassfilter(Odata.trial{1}, fsamp, [frq-bwid frq+bwid], [], 'but', 'twopass', 'reduce');
amp(:,1) = abs(hilbert(Xdata.trial{1}(1,:)));
amp(:,2) = abs(hilbert(Xdata.trial{1}(2,:)));
amp = amp(fsamp:end-fsamp,:); % Truncate Filter Artefacts
phi(:,1) = angle(hilbert(Xdata.trial{1}(1,:)));
phi(:,2) = angle(hilbert(Xdata.trial{1}(2,:)));
phi = phi(fsamp:end-fsamp,:); % Truncate Filter Artefacts
dphi_12 = unwrap(diff(phi,1,2));

% % optional amp weighting
% % ampw = unwrap((amp(:,1).*amp(:,2))./(max(amp(:,1))*max(amp(:,2))));  %%
% % dphi_12 = angle(ampw.*exp(i.*(phi(:,1)-phi(:,2))));%%
% % dphi_12 = (dphi_12-(1/sqrt(length(phi(:,1)))))./(1-(1/sqrt(length(phi(:,1))))); %%
dphi_12_dt = diff(dphi_12);
%% Prepare Continuous Filtered Data
Xdata.trial{1} = Xdata.trial{1}(:,fsamp:end-fsamp);
Xdata.time{1} = Odata.time{1}(:,fsamp:end-fsamp);
Xdata.fsample = fsamp;
