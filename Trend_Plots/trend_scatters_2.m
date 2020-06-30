close
clear Ax
figure(1)

Ax{1} = subplot(121);
cla

gvec = {CLIVAR.namevec(:),repmat({'CMIP'},[size(CMIP.MIZ_ensmean,1) 1])};
gvec = strrep(vertcat(gvec{:}),'_','-');

gscatter([CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist],[CLIVAR.SIE_slope_hist CMIP.SIE_slope_hist],gvec);

ylabel('SIE Trend'); xlabel('GMT Trend');
hold on
grid on; box on;

leg = legend;
set(leg,'AutoUpdate','off')


for i = 1:length(OBS.MIZ_slope_hist)
    
    for j = 1:length(OBS.GMT_slope_hist)
        
        errorbar(OBS.GMT_slope_hist(j),OBS.SIE_slope_hist(i),-OBS.SIE_slope_err(i),OBS.SIE_slope_err(i),-OBS.GMT_slope_err(j),OBS.GMT_slope_err(j),'k','linewidth',1)
        
    end
    
end

gscatter([CLIVAR.GMT_slope_fut CMIP.GMT_slope_fut],[CLIVAR.SIE_slope_fut CMIP.SIE_slope_fut],gvec,'rygcbm','o');


%%

Ax{2} = subplot(122);

gvec = {CLIVAR.namevec(:),repmat({'CMIP'},[size(CMIP.MIZ_ensmean,1) 1])};
gvec = strrep(vertcat(gvec{:}),'_','-');

gscatter([CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist],[CLIVAR.MIZ_slope_hist CMIP.MIZ_slope_hist],gvec);

ylabel('MIZ Trend'); xlabel('GMT Trend');
hold on
grid on; box on;

legend('off');


for i = 1:length(OBS.MIZ_slope_hist)
    
    for j = 1:length(OBS.GMT_slope_hist)
        
        errorbar(OBS.GMT_slope_hist(j),OBS.MIZ_slope_hist(i),-OBS.MIZ_slope_err(i),OBS.MIZ_slope_err(i),-OBS.GMT_slope_err(j),OBS.GMT_slope_err(j),'k','linewidth',1)
        
    end
    
end
title('1979-2014');

gscatter([CLIVAR.GMT_slope_fut CMIP.GMT_slope_fut],[CLIVAR.MIZ_slope_fut CMIP.MIZ_slope_fut],gvec,'rygcbm','o');



delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [10 6];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'trend_scatters_2.pdf')