figure(3)

clear Ax

Ax{1} = subplot(121);

% gvec = {CLIVAR.namevec(:)}; % For labeling box plot
gvec = {repmat({'CMIP'},[size(CMIP.MIZ,1) 1]),CLIVAR.namevec(:)};
gvec = strrep(vertcat(gvec{:}),'_','-');

h = gscatter([ CMIP.MIZ_slope_hist CLIVAR.MIZ_slope_hist],[ CMIP.SIE_slope_hist CLIVAR.SIE_slope_hist],gvec);

for i = 1:length(h)

    set(h(i), 'Color',clabs(i,:),'MarkerFaceColor','none','MarkerSize',20)

end

% gscatter([CLIVAR.GMT_slope_hist],[CLIVAR.SIE_slope_hist],gvec);

ylabel('SIE Trend'); xlabel('MIZ Trend');
hold on
grid on; box on;

leg = legend;
set(leg,'AutoUpdate','off')


for i = 1:length(OBS.MIZ_slope_hist)
    
    
    errorbar(OBS.MIZ_slope_hist(i),OBS.SIE_slope_hist(i),-OBS.SIE_slope_err(i),OBS.SIE_slope_err(i),-OBS.MIZ_slope_err(i),OBS.MIZ_slope_err(i),'k','linewidth',1)
    
    
end
title('1979-2014');

%%

Ax{2} = subplot(122);

h = gscatter([ CMIP.MIZ_slope_fut CLIVAR.MIZ_slope_fut],[ CMIP.SIE_slope_fut CLIVAR.SIE_slope_fut],gvec);

for i = 1:length(h)

    set(h(i), 'Color',clabs(i,:),'MarkerFaceColor','none','MarkerSize',20)

end

ylabel('SIE Trend'); xlabel('MIZ Trend');
hold on
grid on; box on;
title('2020-2055');

legend('off')

set(Ax{1},'xlim',get(Ax{2},'xlim')); 
set(Ax{1},'ylim',get(Ax{2},'ylim')); 

%%



delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [10 4];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'trend_MIZ_SIE.pdf')