clear
close all

addpath('../');

bands_SIE = 2:6;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse

mos_SIE = [9];
mos_MIZ = [9];

obuse = [1 2 3];

horvat_colors;

sigmaval = 3;

qval=norminv(1-exp(-sigmaval));
q75=norminv(.75);
whiskval = (qval - q75)/(2*q75);


% Create sub-ensemble from CMIP
q = squeeze(mean(CMIP.MIZA_F_ensmean(:,mos_SIE,165),2));
num_CMIP_models = sum(~isnan(q));


%%

Ax{1} = subplot(321);

addpath('../boundedlin/');

cmipplot = squeeze(mean(CMIP.SIA_ensmean(:,mos_SIE,130:165),2));
cmipplot = bsxfun(@minus,cmipplot,cmipplot(:,1)); 

gmtplot = squeeze(CMIP.GMT_ensmean(:,130:165));
gmtplot = bsxfun(@minus,gmtplot,gmtplot(:,1)); 


plot(CMIP.plotyrs,sensplot); 
xlim([1979 2014])



%%

%%
Ax{3} = subplot(323);

% gvec = {CLIVAR.namevec(:)}; % For labeling box plot
gvec = {CLIVAR.namevec(:),repmat({'CMIP6-EME'},[size(CMIP.MIZ_ensmean,1) 1])};
gvec = strrep(vertcat(gvec{:}),'_','-');

x1 = [mean(CLIVAR.SIA_slope_hist(:,mos_SIE),2); mean(CMIP.SIA_slope_hist(:,mos_SIE),2)];
x2 = [CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist];

scath = gscatter(x1,x2,gvec);

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

[r,p] = corrcoef(x1,x2);

fprintf('p = %d, corrcoeff = %d \n',p(1,2),r(1,2));


legend('off');

for i = 1:length(scath)
    
    set(scath(i), 'MarkerFaceColor',clabs(i,:),'Marker','o','MarkerEdgeColor','none','MarkerSize',3)
    
end

set(scath(end),'MarkerFaceColor','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3,'linewidth',1)


hold on


scath(end+1) = scatter(mean(OBS.SIA_slope_hist(obuse,mos_MIZ),2),repmat(OBS.GMT_slope_hist,[length(obuse) 1]),100,'x','linewidth',2,'markeredgecolor','k');

ylabel('$\Delta$ GMT ($^\circ$/yr)','interpreter','latex');
title('Sep. 1979-2014','interpreter','latex'); 
xlabel(''); 
% xlabel('(\%/yr)','interpreter','latex');
grid on; box on;

%%

Ax{4} = subplot(324);

x1 = [mean(CLIVAR.MIZA_F_slope_hist(:,mos_MIZ),2); mean(CMIP.MIZA_F_slope_hist(:,mos_MIZ),2)];
x2 = [CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist];

scath = gscatter(x1,x2,gvec);

legend('off');

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

[r,p] = corrcoef(x1,x2);

fprintf('p = %d, corrcoeff = %d \n',p(1,2),r(1,2));


for i = 1:length(scath)
    
    set(scath(i), 'MarkerFaceColor',clabs(i,:),'Marker','o','MarkerEdgeColor','none','MarkerSize',3)
    
end

set(scath(end),'MarkerFaceColor','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3,'linewidth',1)

hold on

    %    errorbar(OBS.MIZA_F_slope_hist(i),OBS.GMT_slope_hist,-OBS.GMT_slope_err,OBS.GMT_slope_err,-OBS.MIZA_F_slope_err(i),OBS.MIZA_F_slope_err(i),'k','linewidth',1)
scath(end+1) = scatter(mean(OBS.MIZA_F_slope_hist(obuse,mos_MIZ),2),repmat(OBS.GMT_slope_hist,[length(obuse) 1]),100,'x','linewidth',2,'markeredgecolor','k');

hold off

title('Sep. 1979-2014','interpreter','latex'); 

xlabel('','interpreter','latex'); ylabel(''); %$\Delta$ GMT ($^\circ$/yr)','interpreter','latex');
grid on; box on;
ylim(get(Ax{3},'ylim'));


%%


Ax{5} = subplot(325);

obsout = mean(OBS.SIA_slope_hist(obuse,mos_SIE),2);

modelout = [mean(CLIVAR.SIA_slope_hist(:,mos_SIE),2); mean(CMIP.SIA_slope_hist(:,mos_SIE),2)];
modelnames = {CLIVAR.namevec(:) repmat({'CMIP6-EME'},[size(CMIP.MIZ_ensmean,1) 1])};
modelnames = vertcat(modelnames{:}); 

gvec = {repmat({'OBS'},[length(obuse),1]) modelnames}; 
gvec = strrep(vertcat(gvec{:}),'_','-');

[~,~,indices] = unique(modelnames,'stable');


x = [obsout; modelout];

for i = 1:6
    tval(i+1) = ttest2(modelout(indices==i),obsout,'Vartype','unequal');
    
    if tval(i+1) == 1
        tcols(i+1,:) = [0.8906    0.1016    0.1094]; 
    else
        tcols(i+1,:) = [0.2148    0.4922    0.7188];
    end

end

tcols(1,:) = [0 0 0]; 



boxplot(x,gvec,'orientation','horizontal','whisker',whiskval,'symbol','','colors',tcols);

xline(median(mean(OBS.SIA_slope_hist(:,mos_SIE),2)),'--k','linewidth',1);

grid on; box on;
xlabel('$\Delta$SIA (\%/yr)','interpreter','latex');

set(gca,'TickLabelInterpreter','latex');
% set(gca,'yticklabel',{});

%%
Ax{6} = subplot(326);

obsout = mean(OBS.MIZA_F_slope_hist(obuse,mos_SIE),2);

modelout = [mean(CLIVAR.MIZA_F_slope_hist(:,mos_SIE),2); mean(CMIP.MIZA_F_slope_hist(:,mos_SIE),2)];

x = [obsout; modelout];

for i = 1:6
    tval(i+1) = ttest2(modelout(indices==i),obsout,'Vartype','unequal');
    
    if tval(i+1) == 1
        tcols(i+1,:) = [0.8906    0.1016    0.1094]; 
    else
        tcols(i+1,:) = [0.2148    0.4922    0.7188];
    end

end

tcols(1,:) = [0 0 0]; 
boxplot(x,gvec,'orientation','horizontal','whisker',whiskval,'symbol','','colors',tcols);




xline(median(mean(OBS.MIZA_F_slope_hist(obuse,mos_MIZ),2)),'--k','linewidth',1);

grid on; box on;

xlabel('$\Delta$ MIZF (\%/yr)','interpreter','latex');
% set(gca,'yticklabel',{});
set(gca,'TickLabelInterpreter','latex');
set(gca,'yaxislocation','right')
%%

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',8,'xminortick','on','yminortick','on','TickLabelInterpreter','latex')
    posy = get(Ax{i},'position');
    
    if i > 3
        set(Ax{i},'position',get(Ax{i},'position')+[0 0 -.035 0 ]);
    end
    
end

pos = [8.5 5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

drawnow

%%
pos1 = get(Ax{5},'position');
pos2 = get(Ax{6},'position');

for i = 1:6
    
    posy = get(Ax{i},'position');
    
    if mod(i,2) == 1
        
        set(Ax{i},'position',[pos1(1) posy(2) pos1(3) posy(4)]);
    else
        set(Ax{i},'position',[pos2(1) posy(2) pos2(3) posy(4)]);
    end
    
    
end

posy = get(Ax{4},'position');
posy2 = get(Ax{3},'position');

set(Ax{4},'position',[posy(1) posy2(2) posy(3) posy2(4)]);
set(Ax{4},'ytick',get(Ax{3},'ytick')); 
%%
tightfig;

set(Ax{3},'position',get(Ax{3},'position') + [0 -.045 0 0]); 
set(Ax{4},'position',get(Ax{4},'position') + [0 -.045 0 0]); 

L = legend(scath,{q{:},'CMIP6-EME','Obs'},'location',[.15 .62 .7 .05],'fontsize',8,'orientation','horizontal');
L.ItemTokenSize = [8 8];
saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-1/Fig-1.pdf')
% saveas(gcf,'Extents_Sept.pdf');
