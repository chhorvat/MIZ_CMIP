figure(1)

Ax{1} = subplot(221);

% gvec = {CLIVAR.namevec(:)}; % For labeling box plot
gvec = {repmat({'CMIP'},[size(CMIP.MIZ,1) 1]),CLIVAR.namevec(:)};
gvec = strrep(vertcat(gvec{:}),'_','-');

h = gscatter([CMIP.GMT_slope_hist CLIVAR.GMT_slope_hist ],[CMIP.SIE_slope_hist CLIVAR.SIE_slope_hist ],gvec);
% gscatter([CLIVAR.GMT_slope_hist],[CLIVAR.SIE_slope_hist],gvec);

for i = 1:length(h)

    set(h(i), 'Color',clabs(i,:),'MarkerFaceColor','none','MarkerSize',10)

end

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
title('1979-2014');
%%

Ax{2} = subplot(222);

gvec = {repmat({'CMIP'},[size(CMIP.MIZ,1) 1]),CLIVAR.namevec(:)};
gvec = strrep(vertcat(gvec{:}),'_','-');

h = gscatter([CMIP.GMT_slope_hist CLIVAR.GMT_slope_hist ],[CMIP.MIZ_slope_hist CLIVAR.MIZ_slope_hist ],gvec);

for i = 1:length(h)

    set(h(i), 'Color',clabs(i,:),'MarkerFaceColor','none','MarkerSize',10)

end

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

%%

Ax{3} = subplot(223);

h = gscatter([CMIP.GMT_slope_fut CLIVAR.GMT_slope_fut ],[CMIP.SIE_slope_fut CLIVAR.SIE_slope_fut],gvec);

for i = 1:length(h)

    set(h(i), 'Color',clabs(i,:),'MarkerFaceColor','none','MarkerSize',10)

end

ylabel('SIE Trend'); xlabel('GMT Trend');
hold on
grid on; box on;
title('Future');
legend('off');

%%

Ax{4} = subplot(224);

h = gscatter([ CMIP.GMT_slope_fut CLIVAR.GMT_slope_fut] ,[ CMIP.MIZ_slope_fut CLIVAR.MIZ_slope_fut],gvec);

for i = 1:length(h)

    set(h(i), 'Color',clabs(i,:),'MarkerFaceColor','none','MarkerSize',10)

end

ylabel('MIZ Trend'); xlabel('GMT Trend');
hold on
grid on; box on;

title('Future');

legend('off');



% set(Ax{2},'xlim',get(Ax{4},'xlim')); set(Ax{4},'ylim',get(Ax{2},'ylim')); 
% set(Ax{1},'xlim',get(Ax{3},'xlim')); set(Ax{1},'ylim',get(Ax{3},'ylim')); 
% 

%%
delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [10 6];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'trend_scatters.pdf')