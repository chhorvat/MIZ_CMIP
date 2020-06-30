clear
close all

addpath('../');

bands_SIE = 2:5;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse

mos_SIE = [9];
mos_MIZ = [9]; 
mos_GMT = 1:12; 

horvat_colors;

[~,~,indices] = unique(CLIVAR.namevec);

%%
Ax{1} = subplot(131);


plot_yrs = 2020:2055;

model_yrs = 171:206; % 2020:2055

clear slope* err*

clivarGMT = squeeze(mean(CLIVAR.GMT(:,mos_GMT,:),2));

% Get GMT trends
for i = 1:size(clivarGMT,1); 
    
    yval = clivarGMT(i,model_yrs); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    slope_GMT(i) = b(2);
    
end

%% Get plot trends
modelplot = squeeze(mean(CLIVAR.SIA(:,mos_SIE,:),2));
% modelplot = cat(1,modelplot,squeeze(mean(CMIP.SIE(:,mos_SIE,:),2)));

for i = 1:size(modelplot,1); 
    
    yval = modelplot(i,model_yrs); 
    % yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope(i) = b(2);
   
    
end


gvec =strrep(CLIVAR.namevec,'_','-');
% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');

h = gscatter(slope,slope_GMT,gvec);

ylabel('Slope in GMT'); 
xlabel('Slope in SIE'); 
hold on
grid on; box on; 

          
legend('off'); 

%%


Ax{2} = subplot(132);

modelplot = squeeze(mean(CLIVAR.MIZA(:,mos_MIZ,:),2));
% modelplot = cat(1,modelplot,squeeze(mean(CMIP.MIZ(:,mos_MIZ,:),2)));

for i = 1:size(modelplot,1); 
    
    yval = modelplot(i,model_yrs); 
    % yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope(i) = b(2);
   
    
end

gvec =strrep(CLIVAR.namevec,'_','-');
% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');

h = gscatter(slope,slope_GMT,gvec);

ylabel('Slope in GMT'); 
xlabel('Slope in MIZ'); 
hold on
grid on; box on; 

legend('off'); 

%% 
Ax{3} = subplot(133);

modelplot = squeeze(mean(CLIVAR.MIZA_F(:,mos_MIZ,:),2));
% modelplot = cat(1,modelplot,squeeze(mean(CMIP.MIZ(:,mos_MIZ,:),2)));

for i = 1:size(modelplot,1); 
    
    yval = modelplot(i,model_yrs); 
    % yval = (yval - yval(1))/yval(1); 

    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err(i) = (bint(2,2) - b(2));
    slope(i) = b(2);
   
    
end

gvec =strrep(CLIVAR.namevec,'_','-');
% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');

h = gscatter(slope,slope_GMT,gvec);

ylabel('Slope in GMT'); 
xlabel('Slope in MIZA Fraction'); 
hold on
grid on; box on; 

legend('off'); 


%%
delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [10 4];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

% saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-2/Fig-2.pdf')
saveas(gcf,'Future-Sept-CLIVAR.pdf'); 