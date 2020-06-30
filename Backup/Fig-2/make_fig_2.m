clear
close all

addpath('../');

bands_SIE = 2:5;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse

mos_SIE = [9];
mos_MIZ = [9]; 

obuse = [1 2 3]; 

horvat_colors;

[~,~,indices] = unique(CLIVAR.namevec);

%%
Ax{1} = subplot(221);

clivarplot = squeeze(mean(CLIVAR.SIE(:,mos_SIE,:),2))/1e6;
obsplot = squeeze(mean(OBS.SIE(:,mos_SIE,:),2))/1e6;

plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

for i = 1:size(clivarplot,1); 
    
    yval = clivarplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_clivar(i) = b(2);
    
end

for i = 1:size(obsplot,1)
    
    yval = obsplot(i,hist_yrs_obs); 
    yval = (yval - yval(1))/yval(1); 
  
    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
end


gvec = vertcat({'Obs','Obs','Obs',CLIVAR.namevec{:}});

boxplot(100*[slope_obs slope_clivar],gvec,'Colorgroup',gvec,'orientation','horizontal','notch','off')

xline(100*mean(slope_obs),'--k','linewidth',1)
title('Sea Ice Extent Trend','interpreter','latex'); 
grid on; box on; 
xlabel('\%/yr','interpreter','latex'); 

Ax{2} = subplot(222);

clivarplot = squeeze(mean(CLIVAR.MIZ(:,mos_MIZ,:),2))/1e6;
obsplot = squeeze(mean(OBS.MIZ(:,mos_MIZ,:),2))/1e6;

plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

for i = 1:size(clivarplot,1); 
    
    yval = clivarplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_clivar(i) = b(2);
    
end

for i = 1:size(obsplot,1)
 
    yval = obsplot(i,hist_yrs_obs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
end


gvec = vertcat({'Obs','Obs','Obs',CLIVAR.namevec{:}});

boxplot(100*[slope_obs slope_clivar],gvec,'Colorgroup',gvec,'orientation','horizontal','notch','off')

xline(100*mean(slope_obs),'--k','linewidth',1)

title('MIZ Extent Trend','interpreter','latex'); 
grid on; box on; 
xlabel('\%/yr','interpreter','latex'); 
set(gca,'yticklabel',{}); 

Ax{3} = subplot(233);

clivarplot = squeeze(mean(CLIVAR.MIZ_F(:,mos_SIE,:),2))/1e6;
obsplot = squeeze(mean(OBS.MIZ_F(:,mos_SIE,:),2))/1e6;

plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

for i = 1:size(clivarplot,1); 
    
    yval = clivarplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_clivar(i) = b(2);
    
end

for i = 1:size(obsplot,1)
    
    yval = obsplot(i,hist_yrs_obs); 
    yval = (yval - yval(1))/yval(1); 
  
    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
end


gvec = vertcat({'Obs','Obs','Obs',CLIVAR.namevec{:}});

boxplot(100*[slope_obs slope_clivar],gvec,'Colorgroup',gvec,'orientation','horizontal','notch','off')

xline(100*mean(slope_obs),'--k','linewidth',1)
title('MIZ Frac Trend','interpreter','latex'); 
grid on; box on; 
xlabel('\%/yr','interpreter','latex');
set(gca,'yticklabel',{}); 

%%
Ax{4} = subplot(234);

clivarplot = squeeze(mean(CLIVAR.SIA(:,mos_MIZ,:),2))/1e6;
cmipplot = squeeze(mean(CMIP.SIA(:,mos_MIZ,:),2))/1e6;
obsplot = squeeze(mean(OBS.SIA(:,mos_MIZ,:),2))/1e6;

plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

for i = 1:size(clivarplot,1); 
    
    yval = clivarplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_clivar(i) = b(2);
    
end

for i = 1:size(cmipplot,1); 
    
    yval = cmipplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_cmip(i) = b(2);
    
end

for i = 1:size(obsplot,1)
 
    yval = obsplot(i,hist_yrs_obs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
    
end

q = repmat({'CMIP'},[size(CMIP.MIZ,1) 1]);
gvec = vertcat({'Obs','Obs','Obs',CLIVAR.namevec{:},q{:}});

boxplot(100*[slope_obs slope_clivar slope_cmip],gvec,'Colorgroup',gvec,'orientation','horizontal','notch','off')

xline(100*mean(slope_obs),'--k','linewidth',1)

title('MIZA Frac Trend','interpreter','latex'); 
grid on; box on; 
xlabel('\%/yr','interpreter','latex'); 
set(gca,'yticklabel',{}); 



Ax{5} = subplot(235);

clivarplot = squeeze(mean(CLIVAR.MIZA(:,mos_MIZ,:),2))/1e6;
cmipplot = squeeze(mean(CMIP.MIZA(:,mos_MIZ,:),2))/1e6;
obsplot = squeeze(mean(OBS.MIZA(:,mos_MIZ,:),2))/1e6;

plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

for i = 1:size(clivarplot,1); 
    
    yval = clivarplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_clivar(i) = b(2);
    
end

for i = 1:size(cmipplot,1); 
    
    yval = cmipplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_cmip(i) = b(2);
    
end

for i = 1:size(obsplot,1)
 
    yval = obsplot(i,hist_yrs_obs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
    
end

q = repmat({'CMIP'},[size(CMIP.MIZ,1) 1]);
gvec = vertcat({'Obs','Obs','Obs',CLIVAR.namevec{:},q{:}});

boxplot(100*[slope_obs slope_clivar slope_cmip],gvec,'Colorgroup',gvec,'orientation','horizontal','notch','off')

xline(100*mean(slope_obs),'--k','linewidth',1)

title('MIZA Trend','interpreter','latex'); 
grid on; box on; 
xlabel('\%/yr','interpreter','latex'); 
set(gca,'yticklabel',{}); 


Ax{6} = subplot(236);

clivarplot = squeeze(mean(CLIVAR.MIZA_F(:,mos_MIZ,:),2))/1e6;
cmipplot = squeeze(mean(CMIP.MIZA_F(:,mos_MIZ,:),2))/1e6;
obsplot = squeeze(mean(OBS.MIZA_F(:,mos_MIZ,:),2))/1e6;

plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

for i = 1:size(clivarplot,1); 
    
    yval = clivarplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_clivar(i) = b(2);
    
end

for i = 1:size(cmipplot,1); 
    
    yval = cmipplot(i,hist_yrs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope_cmip(i) = b(2);
    
end

for i = 1:size(obsplot,1)
 
    yval = obsplot(i,hist_yrs_obs); 
    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
    
end

q = repmat({'CMIP'},[size(CMIP.MIZ,1) 1]);
gvec = vertcat({'Obs','Obs','Obs',CLIVAR.namevec{:},q{:}});

boxplot(100*[slope_obs slope_clivar slope_cmip],gvec,'Colorgroup',gvec,'orientation','horizontal','notch','off')

xline(100*mean(slope_obs),'--k','linewidth',1)

title('MIZA Frac Trend','interpreter','latex'); 
grid on; box on; 
xlabel('\%/yr','interpreter','latex'); 
set(gca,'yticklabel',{}); 

%%
delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [10 5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

% saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-2/Fig-2.pdf')
saveas(gcf,'Fig-2-Winter.pdf'); 