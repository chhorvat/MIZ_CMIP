clear 
close all

plot_yrs = 1979:2014;

bands_SIE = 2:3;
bands_MIZ = 2;
mos_GMT = [1:12];
mos_SIE = [1:12]; 
mos_MIZ = [1:5 10:12]; 

plot_preamble

% Plotting parameters
obuse = [1 2 3]; % Which of the SIC obs to use
smooth_window = 10; % Width of smoothing window
smooth_type = 'loess'; % Type of smoothing

past_avg_model = 130:130+smooth_window; % Years in model time to use for anomaly
past_avg_obs = 1:1+smooth_window; % Years in obs time to use for anomaly

%% Parametric Plot oF SIE vs GMT

clearvars -except OBS CMIP CLIVAR bands* obuse smooth* past_avg*

close all

horvat_colors


% Anomaly from 1850:1950 mean
Ax{1} = subplot(221);


[hCMIP,~] = boundedline(CMIP.plotyrs,smoothdata(nanmean(CMIP.SIE_ensmean,1)/1e6,2,smooth_type,smooth_window), ...
    permute(nanstd(CMIP.SIE_ensmean,[],1)/1e6,[2 3 1]),'alpha'); 

hold on

[hCLIV,~] = boundedline(CLIVAR.plotyrs,smoothdata(CLIVAR.SIE_ensmean/1e6,2,smooth_type,smooth_window), ...
    permute(CLIVAR.SIE_ensstd/1e6,[2 3 1]),'alpha'); 

% plot(CMIP.plotyrs,smoothdata(nanmean(CMIP.SIE_ensmean,1),2,smooth_type,smooth_window)/1e6,'--k','linewidth',1)
hObs = boundedline(OBS.plotyrs,smoothdata(nanmean(OBS.SIE(obuse,:),1)/1e6,2,smooth_type,smooth_window), ...
    permute(nanstd(OBS.SIE(obuse,:),[],1)/1e6,[2 3 1]),'cmap',[0 0 0],'alpha');  

% plot(OBS.plotyrs,smoothdata(OBS.SIE(obuse,:),2,smooth_type,smooth_window)/1e6,'k','linewidth',1); 
xlim([1979 2014]);
ylim([0 15]); 
ylabel('10$^6$ sq km','interpreter','latex');
grid on; box on; set(gca,'xminorgrid','on');
title('Sea Ice Extent','interpreter','latex');


legend([hCMIP; hCLIV; hObs],strrep(vertcat(CLIVAR.names(:),'CMIP6','Obs'),'_','-'),'orientation','horizontal', ...
    'location',[.11 .015 .78 .05]); 

%%
Ax{2} = subplot(222);
cla
anom_CLIVAR_SIE = smoothdata(CLIVAR.SIE_ensmean,2,smooth_type,smooth_window);
anom_CLIVAR_SIE = bsxfun(@minus,anom_CLIVAR_SIE,nanmean(anom_CLIVAR_SIE(:,past_avg_model),2))/1e6; 

anom_CMIP_SIE =  smoothdata(CMIP.SIE_ensmean,2,smooth_type,smooth_window); 
anom_CMIP_SIE = bsxfun(@minus,anom_CMIP_SIE,nanmean(anom_CMIP_SIE(:,past_avg_model),2))/1e6; 
std_CMIP_SIE = nanstd(anom_CMIP_SIE,[],1); 
anom_CMIP_SIE = nanmean(anom_CMIP_SIE,1); 

anom_OBS_SIE = smoothdata(OBS.SIE,2,smooth_type,smooth_window);
anom_OBS_SIE = bsxfun(@minus,anom_OBS_SIE,nanmean(anom_OBS_SIE(:,past_avg_obs),2))/1e6; 


% plot(CLIVAR.plotyrs,anom_CLIVAR_SIE,'linewidth',1)
boundedline(CMIP.plotyrs,anom_CMIP_SIE,permute(anom_CMIP_SIE,[2 3 1]),'alpha'); 
hold on
boundedline(CLIVAR.plotyrs,anom_CLIVAR_SIE,permute(CLIVAR.SIE_ensstd/1e6,[2 3 1]),'alpha'); 
% plot(CMIP.plotyrs,anom_CMIP_SIE,'--k','linewidth',1)
boundedline(OBS.plotyrs,nanmean(anom_OBS_SIE(obuse,:),1),nanstd(anom_OBS_SIE,[],1),'cmap',[0 0 0],'alpha'); 
hold off

% plot(CMIP.plotyrs,anom_CMIP_SIE+std_CMIP_SIE,'--r')
% plot(CMIP.plotyrs,anom_CMIP_SIE-std_CMIP_SIE,'--r')
% plot(OBS.plotyrs,anom_OBS_SIE(obuse,:),'k','linewidth',1); 
xlim([1979 2014]);
ylabel('10$^6$ sq km','interpreter','latex');
grid on; box on; set(gca,'xminorgrid','on');
title('Anomaly Sea Ice Extent','interpreter','latex');

%%

Ax{3} = subplot(223);

boundedline(CMIP.plotyrs,smoothdata(nanmean(CMIP.MIZ_ensmean,1)/1e6,2,smooth_type,smooth_window), ...
    permute(nanstd(CMIP.MIZ_ensmean,[],1)/1e6,[2 3 1]),'alpha'); 

hold on

boundedline(CLIVAR.plotyrs,smoothdata(CLIVAR.MIZ_ensmean/1e6,2,smooth_type,smooth_window), ...
    permute(CLIVAR.MIZ_ensstd/1e6,[2 3 1]),'alpha'); 

% plot(CMIP.plotyrs,smoothdata(nanmean(CMIP.MIZ_ensmean,1),2,smooth_type,smooth_window)/1e6,'--k','linewidth',1)
boundedline(OBS.plotyrs,smoothdata(nanmean(OBS.MIZ(obuse,:),1)/1e6,2,smooth_type,smooth_window), ...
    permute(nanstd(OBS.MIZ(obuse,:),[],1)/1e6,[2 3 1]),'cmap',[0 0 0],'alpha'); 

% plot(OBS.plotyrs,smoothdata(OBS.SIE(obuse,:),2,smooth_type,smooth_window)/1e6,'k','linewidth',1); 
xlim([1979 2014]);
ylim([0 2]); 
ylabel('10$^6$ sq km','interpreter','latex');
grid on; box on; set(gca,'xminorgrid','on');
title('MIZ Extent','interpreter','latex');

%%

Ax{4} = subplot(224);


cla
anom_CLIVAR_MIZ = smoothdata(CLIVAR.MIZ_ensmean,2,smooth_type,smooth_window);
anom_CLIVAR_MIZ = bsxfun(@minus,anom_CLIVAR_MIZ,nanmean(anom_CLIVAR_MIZ(:,past_avg_model),2))/1e6; 

anom_CMIP_MIZ =  smoothdata(CMIP.MIZ_ensmean,2,smooth_type,smooth_window); 
anom_CMIP_MIZ = bsxfun(@minus,anom_CMIP_MIZ,nanmean(anom_CMIP_MIZ(:,past_avg_model),2))/1e6; 
std_CMIP_MIZ = nanstd(anom_CMIP_MIZ,[],1); 
anom_CMIP_MIZ = nanmean(anom_CMIP_MIZ,1); 

anom_OBS_MIZ = smoothdata(OBS.MIZ,2,smooth_type,smooth_window);
anom_OBS_MIZ = bsxfun(@minus,anom_OBS_MIZ,nanmean(anom_OBS_MIZ(:,past_avg_obs),2))/1e6; 


% plot(CLIVAR.plotyrs,anom_CLIVAR_MIZ,'linewidth',1)
% boundedline(CMIP.plotyrs,anom_CMIP_MIZ,permute(anom_CMIP_MIZ,[2 3 1]),'alpha'); 
hold on
boundedline(CLIVAR.plotyrs,anom_CLIVAR_MIZ,permute(CLIVAR.MIZ_ensstd/1e6,[2 3 1]),'alpha'); 
% plot(CMIP.plotyrs,anom_CMIP_MIZ,'--k','linewidth',1)
boundedline(OBS.plotyrs,nanmean(anom_OBS_MIZ(obuse,:),1),nanstd(anom_OBS_MIZ,[],1),'cmap',[0 0 0],'alpha'); 
hold off

% plot(CMIP.plotyrs,anom_CMIP_MIZ+std_CMIP_MIZ,'--r')
% plot(CMIP.plotyrs,anom_CMIP_MIZ-std_CMIP_MIZ,'--r')
% plot(OBS.plotyrs,anom_OBS_MIZ(obuse,:),'k','linewidth',1); 
xlim([1979 2014]);
ylabel('10$^6$ sq km','interpreter','latex');
grid on; box on; set(gca,'xminorgrid','on');
title('Anomaly Sea Ice Extent','interpreter','latex');

%%


delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',12,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [16 8];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'timeseries_anomaly.pdf')


