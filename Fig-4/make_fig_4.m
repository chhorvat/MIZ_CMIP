clear
close all

addpath('../');

bands_SIE = 2:6;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse

get_trends;

mos = 9; 

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}; 

obuse = [1 2 3];

horvat_colors;

sigmaval = 3;

qval=norminv(1-exp(-sigmaval));
q75=norminv(.75);
whiskval = (qval - q75)/(2*q75);


% Create sub-ensemble from CMIP
q = squeeze(mean(CMIP.MIZA_F_ensmean(:,mos,200),2));
num_CMIP_models = sum(~isnan(q));

%%
Ax{1} = subplot(221);

% gvec = {CLIVAR.namevec(:)}; % For labeling box plot
gvec = CLIVAR.namevec;
gvec = strrep(gvec,'_','-');

x1 = mean(CLIVAR.MIZA_F_slope_fut(:,mos),2);
x2 = CLIVAR.GMT_slope_fut;

scath = gscatter(x1,x2,gvec);

legend('off');

x2 = x2(~isnan(x1));
x1 = x1(~isnan(x1));

[r,p] = corrcoef(x1,x2);

fprintf('p = %d, corrcoeff = %d \n',p(1,2),r(1,2));


for i = 1:length(scath)
    
    set(scath(i), 'MarkerFaceColor',clabs(i,:),'Marker','o','MarkerEdgeColor','none','MarkerSize',3)
    
end

% xlabel('$\Delta$ MIZ \% (\%/yr)','interpreter','latex'); 
xlabel(''); ylabel(''); %$\Delta$ GMT ($^\circ$/yr)','interpreter','latex');
grid on; box on;
ylim(get(Ax{1},'ylim'));
title([months{mos} '. 2020-2055'],'interpreter','latex'); 
ylabel('$\Delta$ GMT ($^\circ$/yr)','interpreter','latex'); 

%%

Ax{2} = subplot(223);

gvec = CLIVAR.namevec; 
gvec = strrep(gvec,'_','-');

modelout = mean(CLIVAR.MIZA_F_slope_fut(:,mos),2);

boxplot(modelout,gvec,'orientation','horizontal','whisker',whiskval,'symbol','','colors',[0 0 0]);

% title('MIZ \% Trend','interpreter','latex');
grid on; box on;

xlabel('$\Delta$MIZF (\%/yr)','interpreter','latex');
% set(gca,'yticklabel',{});
set(gca,'TickLabelInterpreter','latex');
set(gca,'yaxislocation','left')
%%
Ax{3} = subplot(222);

mos = 12;

% gvec = {CLIVAR.namevec(:)}; % For labeling box plot
gvec = CLIVAR.namevec;
gvec = strrep(gvec,'_','-');

x1 = mean(CLIVAR.MIZA_F_slope_fut(:,mos),2);
x2 = CLIVAR.GMT_slope_fut;

scath = gscatter(x1,x2,gvec);

legend('off');

x2 = x2(~isnan(x1));
x1 = x1(~isnan(x1));

[r,p] = corrcoef(x1,x2);

fprintf('p = %d, corrcoeff = %d \n',p(1,2),r(1,2));


for i = 1:length(scath)
    
    set(scath(i), 'MarkerFaceColor',clabs(i,:),'Marker','o','MarkerEdgeColor','none','MarkerSize',3)
    
end

% xlabel('$\Delta$ MIZ \% (\%/yr)','interpreter','latex');
xlabel(''); ylabel(''); %$\Delta$ GMT ($^\circ$/yr)','interpreter','latex');
grid on; box on;
ylim(get(Ax{1},'ylim'));
title([months{mos} '. 2020-2055'],'interpreter','latex'); 
set(gca,'yticklabel',{})
% ylabel('$\Delta$ GMT ($^\circ$/yr)','interpreter','latex'); 

%%

Ax{4} = subplot(224);

gvec = CLIVAR.namevec; 
gvec = strrep(gvec,'_','-');

modelout = mean(CLIVAR.MIZA_F_slope_fut(:,mos),2);

boxplot(modelout,gvec,'orientation','horizontal','whisker',whiskval,'symbol','','colors',[0 0 0]);

% title('MIZ \% Trend','interpreter','latex');
grid on; box on;

xlabel('$\Delta$MIZF (\%/yr)','interpreter','latex');
set(gca,'yticklabel',{});
set(gca,'TickLabelInterpreter','latex');
set(gca,'yaxislocation','left')

%%
delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',8,'xminortick','on','yminortick','on','TickLabelInterpreter','latex')
    posy1 = get(Ax{i},'position');
    
end

pos = [6.5 3];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

drawnow

%%
posy1 = get(Ax{1},'position');
posy2 = get(Ax{2},'position');
% 
set(Ax{1},'position',[posy2(1) posy1(2) posy2(3) posy1(4)]);
set(Ax{1},'xlim',get(Ax{2},'xlim')); 

posy = get(Ax{3},'position');
posy2 = get(Ax{4},'position');
% 
set(Ax{3},'position',[posy2(1) posy(2) posy1(3) posy(4)]);
set(Ax{4},'position',[posy2(1) posy2(2) posy1(3) posy2(4)]);

set(Ax{4},'xlim',get(Ax{2},'xlim')); 
set(Ax{3},'xlim',get(Ax{1},'xlim')); 
% %%
tightfig;
%%

delete(findall(gcf,'Tag','legtag'))

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(Ax)

    posy = get(Ax{i},'position');
    
    annotation('textbox',[posy(1) posy(2)+posy(4)-.025 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',10,'Tag','legtag');
    
end

% L = legend(scath,{q{:},'CMIP6-EME'},'location',[.15 .63 .7 .05],'fontsize',8,'orientation','horizontal');
% L.ItemTokenSize = [8 8];
saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-4/Fig-4.pdf')
% saveas(gcf,'Extents_Sept.pdf');
