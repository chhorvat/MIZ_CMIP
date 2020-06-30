figure(4)

h = scatterhist([CMIP.GMT_slope_hist CLIVAR.GMT_slope_hist ],[CMIP.SIE_slope_hist CLIVAR.SIE_slope_hist ],'group',gvec,'kernel','on');

%%

figure(5)

clear Ax

Ax{1} = subplot(121);
cla
[~,~,nameind] = unique(CLIVAR.namevec); 

bins = linspace(-.05,.05,10); 

    h = histfit([CLIVAR.MIZ_slope_hist CMIP.MIZ_slope_hist],15,'kernel'); 
    hold on
    
    h(1).FaceAlpha = 0.25; 
    h(1).FaceColor = clabs(1,:); 

for i = 1:length(OBS.SIE_slope_hist)
    
    xline(OBS.MIZ_slope_hist(i),'k','linewidth',1);
    
end

%%

Ax{2} = subplot(122);
cla
[~,~,nameind] = unique(CLIVAR.namevec); 

bins = linspace(-.05,.05,10); 

    h = histfit([CLIVAR.MIZ_slope_fut CMIP.MIZ_slope_fut],15,'kernel'); 
    hold on
    
    h(1).FaceAlpha = 0.25; 
    h(1).FaceColor = clabs(1,:); 


