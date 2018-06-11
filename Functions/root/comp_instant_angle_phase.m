function [amp phi dphi_12 dphi_12_dt Xdata] = comp_instant_angle_phase(Odata,frq,stn_lb_frq,bwid,band)
fsamp = Odata.fsample;
Xdata.trial{1} = ft_preproc_bandpassfilter(Odata.trial{1}, fsamp, [frq-bwid(band) frq+bwid(band)], [], 'but', 'twopass', 'reduce');
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
 
% % if LowAmpFix == 1
% %     TAmp =  amp(:,1).* amp(:,2);
% %     TAmpNorm = (TAmp./mean(TAmp));
% %     TAmpNormNeg = TAmpNorm;
% %     TAmpNormNeg(TAmpNorm>1) = 1;
% %     dphi_12_dt = dphi_12_dt.*TAmpNormNeg(2:end);
% % end

%% lb_stn_bandpass
betaS.trial{1} = ft_preproc_bandpassfilter(Odata.trial{1}, fsamp, [stn_lb_frq-bwid(2) stn_lb_frq+bwid(2)], [], 'but', 'twopass', 'no');
Xdata.trial{1}(3,:) = betaS.trial{1}(strncmp(Odata.label,'STN',3),:);
az = abs(hilbert(betaS.trial{1}(2,:)));
amp(:,3) = az(fsamp:end-fsamp); % % Truncate Filter Artefacts

%% Prepare Continuous Filtered Data
Xdata.trial{1} = Xdata.trial{1}(:,fsamp:end-fsamp);
Xdata.time{1} = Odata.time{1}(:,fsamp:end-fsamp);
Xdata.fsample = fsamp;
