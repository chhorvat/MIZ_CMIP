clear
close all

addpath('../');

bands_SIE = 2:6;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse

get_trends;

gvec = {CLIVAR.namevec(:),repmat({'CMIP6-EME'},[size(CMIP.MIZ_ensmean,1) 1])};
gvec = strrep(vertcat(gvec{:}),'_','-');


for i = 1:12
    
    
    x1 = [mean(CLIVAR.SIA_slope_hist(:,i),2); mean(CMIP.SIA_slope_hist(:,i),2)];
    x2 = [CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist]';
    x1 = x1(~isnan(x2));
    x2 = x2(~isnan(x2));
    
    [r_SIA(i,:,:),p_SIA(i,:,:)] = corrcoef(zscore(x1),zscore(x2));

    x1 = [mean(CLIVAR.SIE_slope_hist(:,i),2); mean(CMIP.SIE_slope_hist(:,i),2)];
    x2 = [CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist]';

    x1 = x1(~isnan(x2));
    x2 = x2(~isnan(x2));
    
    [r_SIE(i,:,:),p_SIE(i,:,:)] = corrcoef(zscore(x1),zscore(x2));
    
    x1 = [mean(CLIVAR.MIZA_F_slope_hist(:,i),2); mean(CMIP.MIZA_F_slope_hist(:,i),2)];
    x2 = [CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist]';
        
    x1 = x1(~isnan(x2));
    x2 = x2(~isnan(x2));
    
    
    [r_MIZ(i,:,:),p_MIZ(i,:,:)] = corrcoef(zscore(x1),zscore(x2));
    
    
end

%%
SIA_ann = [mean(CLIVAR.SIA_slope_hist,2); mean(CMIP.SIA_slope_hist,2)];
SIE_ann = [mean(CLIVAR.SIE_slope_hist,2); mean(CMIP.SIE_slope_hist,2)];
MIZ_ann = [mean(CLIVAR.MIZA_F_slope_hist,2); mean(CMIP.MIZA_F_slope_hist,2)];
GMT_ann = [CLIVAR.GMT_slope_hist CMIP.GMT_slope_hist]';

SIA_ann = SIA_ann(~isnan(SIA_ann)); 
SIE_ann = SIE_ann(~isnan(SIE_ann)); 
MIZ_ann = MIZ_ann(~isnan(MIZ_ann)); 
GMT_ann = GMT_ann(~isnan(GMT_ann)); 

[r,p] = corrcoef(zscore(SIA_ann),zscore(GMT_ann)); 
[r2,p2] = corrcoef(zscore(MIZ_ann),zscore(GMT_ann)); 
[r3,p3] = corrcoef(zscore(SIA_ann),zscore(MIZ_ann)); 
[r4,p4] = corrcoef(zscore(SIE_ann),zscore(GMT_ann)); 

%% 
close all

figure
horvat_colors; 

Ax{1} = subplot(1,1,1); 

mos = 1:12; 
monames = {'J','F','M','A','M','J','J','A','S','O','N','D'}; 

p1 = plot(mos,r_SIA(:,1,2).^2,'color',clabs(1,:)); 
hold on
p2 = plot(mos,r_MIZ(:,1,2).^2,'color',clabs(2,:)); 
% p3 = plot(mos,r_SIE(:,1,2).^2,'--','color',clabs(1,:)); 

sig_SIA = p_SIA(:,1,2) < .01; 
sig_MIZ = p_MIZ(:,1,2) < .01; 


for i = 1:12
    
    if sig_SIA(i)
        scatter(mos(i),r_SIA(i,1,2).^2,50,'filled','markerfacecolor',clabs(1,:),'markeredgecolor',clabs(1,:)); 
    else
        scatter(mos(i),r_SIA(i,1,2).^2,'markerfacecolor','none','markeredgecolor',clabs(1,:)); 
    end
    
    if sig_MIZ(i)
        scatter(mos(i),r_MIZ(i,1,2).^2,50,'filled','markerfacecolor',clabs(2,:),'markeredgecolor',clabs(2,:));
    else
        scatter(mos(i),r_MIZ(i,1,2).^2,'markerfacecolor','none','markeredgecolor',clabs(2,:));
    end
    
end


plot(mos,r_MIZ(:,1,2).^2,'color',clabs(2,:)); 

xlim([1 12]);
ylim([0 1]); 
set(gca,'xtick',mos,'xticklabel',monames); 
set(gca,'ytick',[0 .25 .5 .75 1]); 
grid on; box on;
ylabel('r$^2$','interpreter','latex')

L = legend([p1 p2],{'SIA','MIZ'},'location','best')
L.ItemTokenSize = [8 8];
set(L,'location','best'); 

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',8,'xminortick','on','yminortick','on','TickLabelInterpreter','latex')
    
end

pos = [3 2];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-2/Fig-2.pdf')
