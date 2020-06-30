clear
close all

addpath('../');

bands_SIE = 2:5;
bands_MIZ = 2:4;

plot_preamble;

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse


obuse = [1 2 3];

horvat_colors;

[~,~,indices] = unique(CLIVAR.namevec);

%%
plot_yrs = 2020:2049;
model_yrs = 171:200;  %
early = (1:5);

clear slope* err*

GMT = cat(1,squeeze(mean(CLIVAR.GMT,2)),CMIP.GMT);

gvec = {CLIVAR.namevec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gvec = vertcat(gvec{:});

obsGMT = OBS.GMT;

% Get GMT trends
for i = 1:size(GMT,1);
    
    yval = GMT(i,model_yrs);
    
    [b, bint] = regress(yval',[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    slope_GMT(i) = b(2);
    
end

%% Get plot trends

model = cat(1,CLIVAR.SIE,CMIP.SIE);

for i = 1:size(model,1)
    for j = 1:12
        
        yval = squeeze(model(i,j,model_yrs));
        % yval = (yval - mean(yval(early)))/mean(yval(early));
        
        [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope(i,j) = b(2);
        
    end
end

sens_SIE = bsxfun(@rdivide,slope,slope_GMT');

% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');
gnam{end+1} = 'CMIP';

[~,~,indices] = unique(gvec);

for i = 1:6
    
    Ax{i} = subplot(2,6,i);
    
    boxplot(sens_SIE(indices==i,:));
    hold on
    
    ylim([-5 0]);
    
    title(gnam{i},'interpreter','latex');
    grid on;
    box on;
    
end

%%

model = cat(1,CLIVAR.MIZ,CMIP.MIZ);


for i = 1:size(model,1)
    for j = 1:12
        
        yval = squeeze(model(i,j,model_yrs));
        % yval = (yval - mean(yval(early)))/mean(yval(early));
        
        [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_MIZ(i,j) = b(2);
        
    end
end

sens_MIZ = bsxfun(@rdivide,slope_MIZ,slope_GMT');

% gvec = {gvec,repmat({'CMIP'},[size(CMIP.MIZ,1) 1])};
gnam = strrep(CLIVAR.names,'_','-');
gnam{end+1} = 'CMIP';

[~,~,indices] = unique(gvec);

for i = 1:6
    
    Ax{6+i} = subplot(2,6,6+i);
    
    boxplot(sens_MIZ(indices==i,:));
    hold on
    
    ylim([-1 1]);
    
    title(gnam{i},'interpreter','latex');
    grid on;
    box on;
    
end


%%
delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
end

pos = [12 5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

% saveas(gcf,'/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/Marginal-Ice-Zone-CMIP/Fig-2/Fig-2.pdf')
saveas(gcf,'Future-Sensitivities-bymodel.pdf');

%%
%%
figure

movec = 'JFMAMJJASOND';

for i = 1:12
    
    Ax{i} = subplot(2,12,i);
    
    boxplot(sens_SIE(:,i),gvec);
    hold on
        
    ylim([-5 0]);
    
    grid on;
    box on;
    set(gca,'xticklabel',{});
    title(movec(i),'interpreter','latex');
end

for i = 1:12
    
    Ax{12+i} = subplot(2,12,12+i);
    
    boxplot(sens_MIZ(:,i),gvec);
    hold on
        
    ylim([-1 1]);
    
    grid on;
    box on;
    set(gca,'xticklabel',{});
    
end

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
end

pos = [24 5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'Future-Sensitivities-bymo.pdf');
