function SNR_Inspector(R,timev,SNR_sw,SNR_eps_z,band,titular)
R.bandnames = R.bandname;
for i = 1:size(SNR_sw,2)
    cmap = linspecer(4);
    ax(i) = subplot(1,numel(SNR_eps_z),i);
    a(1) = plot(timev,SNR_sw(:,i),'color',cmap(i,:)); hold on;
    a(2) = plot(timev([1 end]),[SNR_eps_z(i) SNR_eps_z(i)],'color',cmap(i,:));
    xlabel('Time'); ylabel('dB'); ylim([-15 5]);% legend({'Signal',[R.PA.SNReps_prctile 'th Percentile']})
    if i == 1
    title([titular{1} ' ' R.bandnames{band} ' band SNR'])
    elseif i == 2
    title([titular{2} ' ' R.bandnames{band} ' band SNR'])
    elseif i == 3
    title([titular{2} ' ' R.bandnames{band-1} ' band SNR'])
        elseif i == 4
    title([titular{2} ' ' R.bandnames{band-1} ' band SNR'])
    end
    grid on 
end
linkaxes(ax,'xy');
set(gcf,'Position',[351         771        1520         327]); shg