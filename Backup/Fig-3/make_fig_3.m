clear
close all

addpath('../');

bands_SIE = 2:5;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse

mos_SIE = [1 2 12];
mos_MIZ = [1 2 12]; 
mos_GMT = 1:12; 

obuse = [1 2 3]; 

horvat_colors;

[~,~,indices] = unique(CLIVAR.namevec);

%%
Ax{1} = subplot(131);


plot_yrs = 1979:2014;

hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

clear slope* err*

clivarGMT = squeeze(mean(CLIVAR.GMT(:,mos_GMT,:),2));
obsGMT = OBS.GMT; 

% Get GMT trends
for i = 1:size(clivarGMT,1); 
    
    yval = clivarGMT(i,hist_yrs); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    slope_GMT(i) = b(2);
    
end

yval = obsGMT(hist_yrs_obs); 
[b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);

err_GMT_obs = bint(2,2) - b(2); 
slope_GMT_obs = b(2); 

%% Get plot trends
modelplot = squeeze(mean(CLIVAR.SIA(:,mos_SIE,:),2));
% modelplot = cat(1,modelplot,squeeze(mean(CMIP.SIE(:,mos_SIE,:),2)));
obsplot = squeeze(mean(OBS.SIA(:,mos_SIE,:),2));

for i = 1:size(modelplot,1); 
    
    yval = modelplot(i,hist_yrs); 
%    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope(i) = b(2);
   
    
end

for i = 1:size(obsplot,1)
    
    yval = obsplot(i,hist_yrs_obs); 
%    yval = (yval - yval(1))/yval(1); 
  
    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
    
end

gvec =strrep(CLIVAR.namevec,'_','-');
% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');

h = gscatter(slope,slope_GMT,gvec);

ylabel(' \Delta GMT (1/yr)'); 
xlabel(' \Delta SIA (1/yr)'); 
hold on

           
errobs = std(slope_obs,1,2); 

errorbar(mean(slope_obs),slope_GMT_obs,-err_GMT_obs,err_GMT_obs,-errobs,errobs,'k','linewidth',1)

legend('off'); 
grid on; 
box on; 
%%


Ax{2} = subplot(132);

modelplot = squeeze(mean(CLIVAR.MIZA(:,mos_MIZ,:),2));
% modelplot = cat(1,modelplot,squeeze(mean(CMIP.MIZ(:,mos_MIZ,:),2)));
obsplot = squeeze(mean(OBS.MIZA(:,mos_MIZ,:),2));

for i = 1:size(modelplot,1); 
    
    yval = modelplot(i,hist_yrs); 
%    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope(i) = b(2);
   
    
end

for i = 1:size(obsplot,1)
    
    yval = obsplot(i,hist_yrs_obs); 
%    yval = (yval - yval(1))/yval(1); 
  
    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
    
end

gvec =strrep(CLIVAR.namevec,'_','-');
% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');

h = gscatter(slope,slope_GMT,gvec);

ylabel(' \Delta GMT (1/yr)'); 
xlabel(' \Delta MIZA (1/yr)'); 
hold on

           
errobs = std(slope_obs,1,2); 

errorbar(mean(slope_obs),slope_GMT_obs,-err_GMT_obs,err_GMT_obs,-errobs,errobs,'k','linewidth',1)

legend('off'); 
grid on; 
box on; 
%%
Ax{2} = subplot(133);

modelplot = squeeze(mean(CLIVAR.MIZA_F(:,mos_MIZ,:),2));
% modelplot = cat(1,modelplot,squeeze(mean(CMIP.MIZ(:,mos_MIZ,:),2)));
obsplot = squeeze(mean(OBS.MIZA_F(:,mos_MIZ,:),2));

for i = 1:size(modelplot,1); 
    
    yval = modelplot(i,hist_yrs); 
%    yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope(i) = b(2);
   
    
end

for i = 1:size(obsplot,1)
    
    yval = obsplot(i,hist_yrs_obs); 
%    yval = (yval - yval(1))/yval(1); 
  
    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs(i) = (bint(2,2) - b(2));
    slope_obs(i) = b(2);
    
end

gvec =strrep(CLIVAR.namevec,'_','-');
% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');

h = gscatter(slope,slope_GMT,gvec);

ylabel(' \Delta GMT (1/yr)'); 
xlabel(' \Delta MIZ Area Frac (1/yr)'); 
hold on

           
errobs = std(slope_obs,1,2); 

errorbar(mean(slope_obs),slope_GMT_obs,-err_GMT_obs,err_GMT_obs,-errobs,errobs,'k','linewidth',1)

legend('off'); 
grid on; 
box on; 

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
saveas(gcf,'Fig-3-Winter-CLIVAR.pdf'); 