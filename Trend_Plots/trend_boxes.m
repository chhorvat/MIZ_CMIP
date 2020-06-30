figure(2)

gvec = vertcat([CLIVAR.namevec(:); repmat({'CMIP6'},[size(CMIP.MIZ,1),1])]); % For labeling box plot
gvec = strrep(gvec,'_','-');

Ax{1} = subplot(221);

boxplot([CLIVAR.SIE_slope_hist CMIP.SIE_slope_hist],gvec,'Colorgroup',gvec,'orientation','horizontal')
hold on
xline(0,'--k','linewidth',1)
grid on
for i = 1:length(OBS.SIE_slope_hist)
    
    xline(OBS.SIE_slope_hist(i),'r','linewidth',1);
    
end

title('Hist SIE');

Ax{2} = subplot(222);

boxplot([CLIVAR.MIZ_slope_hist CMIP.MIZ_slope_hist],gvec,'Colorgroup',gvec,'orientation','horizontal')
hold on
xline(0,'--k','linewidth',1)
grid on
for i = 1:length(OBS.MIZ_slope_hist)
    
    xline(OBS.MIZ_slope_hist(i),'r','linewidth',1);
    
end
title('Hist MIZ');

Ax{3} = subplot(223);

gvec = vertcat([CLIVAR.namevec(:);repmat({'CMIP6'},[size(CMIP.MIZ,1),1])]); % For labeling box plot
gvec = strrep(gvec,'_','-');

boxplot([CLIVAR.SIE_slope_fut CMIP.SIE_slope_fut],gvec,'Colorgroup',gvec,'orientation','horizontal')
hold on
xline(0,'--k','linewidth',1)
grid on; box on;
title('Future SIE');
% xlim([0.5 6.5]);

Ax{4} = subplot(224);


boxplot([CLIVAR.MIZ_slope_fut CMIP.MIZ_slope_fut],gvec,'Colorgroup',gvec,'orientation','horizontal')
hold on
xline(0,'--k','linewidth',1)
grid on; box on;
title('Future MIZ');
% xlim([0.5 6.5]);

set(Ax{2},'xlim',get(Ax{4},'xlim')); 
set(Ax{1},'xlim',get(Ax{3},'xlim')); 

pos = [16 6];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
saveas(gcf,'trend_boxes.pdf');