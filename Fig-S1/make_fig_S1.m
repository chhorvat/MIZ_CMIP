clear
close all

addpath('../');

bands_SIE = 2:6;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse


obuse = [1 2 3];

horvat_colors;

[~,~,indices] = unique(CLIVAR.namevec,'stable');


sigmaval = 2;
qval=norminv(1-exp(-sigmaval));
q75=norminv(.75);
whiskval = (qval - q75)/(2*q75);

%%
plot_yrs = 1979:2014;
model_yrs = 130:165; % 1979-2014
obs_yrs = 1:36; % 1979-2014


clear slope* err*

GMT = cat(1,squeeze(mean(CLIVAR.GMT,2)),squeeze(CMIP.GMT_ensmean));

gvec = {CLIVAR.namevec,repmat({'CMIP-EME'},[size(CMIP.MIZ_ensmean,1) 1])};
gvec = vertcat(gvec{:});

obsGMT = OBS.GMT;

slope_GMT = nan(size(GMT,1),1);

% Get GMT trends
for i = 1:size(GMT,1)
    
    yval = GMT(i,model_yrs);
    
    if sum(isnan(yval)) == 0
        [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_GMT(i) = b(2);
    end
    
end

yval = obsGMT(obs_yrs);
[b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);

err_GMT_obs = bint(2,2) - b(2);
slope_GMT_obs = b(2);




%% Get plot trends

model = cat(1,CLIVAR.SIA,CMIP.SIA_ensmean);

obsplot = OBS.SIA(obuse,:,:);

slope = nan(size(model,1),12);
slope_obs = nan(size(obsplot,1),12);

for i = 1:size(model,1)
    for j = 1:12
        
        yval = squeeze(model(i,j,model_yrs));
        % yval = 100*(yval - yval(1))/yval(1);
        
        if sum(isnan(yval)) == 0
            
            [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
            slope(i,j) = b(2);
            
        end
        
    end
end

for i = 1:size(obsplot,1)
    for j = 1:12
        
        yval = squeeze(obsplot(i,j,obs_yrs));
        % yval = 100*(yval - yval(1))/yval(1);
        
        
        if sum(isnan(yval)) == 0
            
            [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
            slope_obs(i,j) = b(2);
            
        end
        
    end
end

sens_SIE = bsxfun(@rdivide,slope,slope_GMT);
sens_SIE_obs = slope_obs/slope_GMT_obs;

% gvec = {gvec,repmat({'CMIP-EME'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');
gnam{end+1} = 'CMIP-EME';

[a,b,indices] = unique(gvec,'stable');

nob = length(obuse);

obvar = stdcorr(sens_SIE_obs,1,nob);
obmean = nanmean(sens_SIE_obs,1);

tval = nan(6,12);

for i = 1:6
    for j = 1:12
        
        tval(i,j) = ttest2(sens_SIE(indices==i,j),sens_SIE_obs(:,j),'Vartype','unequal');
        
    end
end




for i = 1:6
    
    Ax{i} = subplot(2,6,i);
    
    tcols = repmat([0.2148    0.4922    0.7188],[12 1]);
    tcols(tval(i,:) == 1,:) = repmat([0.8906    0.1016    0.1094],[sum(tval(i,:)) 1]);
    
    boxplot(sens_SIE(indices==i,:),'boxstyle','filled','widths',1,'whisker',whiskval,'symbol','','colors',tcols);
    hold on
    %   plot(mean(sens_SIE_obs),'k');
    
    ylim([-5 0]);
    
    title(gnam{i},'interpreter','latex','fontsize',8);
    grid on;
    box on;
    
    %   plot(sens_SIE_obs','k');
    
    plot(obmean,'k');
    plot(obmean + obvar,'--k');
    plot(obmean - obvar,'--k');
    
    set(gca,'xtick',[1:12],'xticklabel',{});
    
    if i > 1
        set(gca,'yticklabel',{});
    else
        ylabel('$10^6$km$^2$/deg K','interpreter','latex');
    end
    
end

obvar1 = obvar; 

%%

model = cat(1,CLIVAR.MIZA_F,CMIP.MIZA_F_ensmean);

obsplot = OBS.MIZA_F(obuse,:,:);

slope = nan(size(model,1),12);
slope_obs = nan(size(obsplot,1),12);


for i = 1:size(model,1)
    for j = 1:12
        
        yval = squeeze(model(i,j,model_yrs));
        % yval = 100*(yval - yval(1))/yval(1);
        if sum(isnan(yval)) == 0
            [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
            slope(i,j) = b(2);
            
        end
    end
end

for i = 1:size(obsplot,1)
    for j = 1:12
        
        yval = squeeze(obsplot(i,j,obs_yrs));
        % yval = 100*(yval - yval(1))/yval(1);
        
        if sum(isnan(yval)) == 0
            
            [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
            slope_obs(i,j) = b(2);
            
        end
        
        
    end
end

sens_MIZ = bsxfun(@rdivide,slope,slope_GMT);
sens_MIZ_obs = slope_obs/slope_GMT_obs;

% gvec = {gvec,repmat({'CMIP-EME'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');
gnam{end+1} = 'CMIP-EME';

[~,~,indices] = unique(gvec,'stable');

obmean =  mean(sens_MIZ_obs);
obvar = stdcorr(sens_MIZ_obs,1,nob);

for i = 1:6
    for j = 1:12

        tval(i,j) = ttest2(sens_MIZ(indices==i,j),sens_MIZ_obs(:,j),'Vartype','unequal');
                
    
    end
end





for i = 1:6
    %%
    Ax{6+i} = subplot(2,6,6+i);
    
    
    tcols = repmat([0.2148    0.4922    0.7188],[12 1]);
    tcols(tval(i,:) == 1,:) = repmat([0.8906    0.1016    0.1094],[sum(tval(i,:)) 1]);
    
    boxplot(sens_MIZ(indices==i,:),'boxstyle','filled','widths',1,'whisker',whiskval,'symbol','','colors',tcols);
    hold on
    %   plot(mean(sens_SIE_obs),'k');
    
    %   title(gnam{i},'interpreter','latex');
    plot(obmean,'k');
    plot(obmean + obvar,'--k');
    plot(obmean - obvar,'--k');
    
    grid on;
    box on;
    
    %   plot(sens_SIE_obs','k');
    
    
    set(gca,'xtick',[1:12],'xticklabel',{'','','M','','','J','','','S','','','D'});
    
    if i > 1
        set(gca,'yticklabel',{});
    else
        ylabel('\%/deg K','interpreter','latex');
    end
    
    ylim([-5 40]);
    
end


%%
delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',8,'xminortick','off','yminortick','on')
    posy = get(Ax{i},'position');
end

pos = [6.5 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');


delx = .01;
dely = .05;

xw = (.9-5*delx)/6;

yw = (.8 - 1*dely)/2;

for i = 1:6
    
    postop = [.075+(i-1)*delx + (i-1)*xw .1 + dely + yw xw yw];
    set(Ax{i},'position',postop);
    
    posbot = [.075+(i-1)*delx + (i-1)*xw .1 xw yw];
    set(Ax{i+6},'position',posbot);
    
    
end

drawnow

saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-S1/Fig-S1.pdf')
saveas(gcf,'Sensitivities-by-model.pdf');

%%
