
clear
close all

load('/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/SSMI/GMT_obs'); 

plot_yrs = 1979:2018;

clabs = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    255,255,51
    166,86,40]/256;

bands_SIE = 2:5;
bands_MIZ = 2:4;
mos = [1:12];
mos_GMT = [1:12]; 


%% Start with CLIVAR LENS

load('../ensemble_MIZ.mat');

[slope_MIZ_CLIVAR,slope_SIZ_CLIVAR,slope_PIZ_CLIVAR,slope_SIE_CLIVAR,slope_GMT_CLIVAR,namevec] = deal(cell(1,length(ensemble_names)));

% Look at each ensemble member
for i = 1:length(ensemble_names)
    
    % Number of ensemble members for this grouping
    nens = size(MIZ_band{i},1);
    % For grouping later
    namevec{i} = repmat({ensemble_names{i}},[nens,1]);
    
    % Get average sea ice extent for certain months over time.
    dum = reshape(sum(MIZ_band{i}(:,bands_SIE,:),2),size(MIZ_band{i},1),12,[]);
    SIE = squeeze(mean(dum(:,mos,:),2));
    
    % Get average MIZ extent for certain months over time.
    dum = reshape(sum(MIZ_band{i}(:,bands_MIZ,:),2),size(MIZ_band{i},1),12,[]);
    MIZ_ext = squeeze(mean(dum(:,mos,:),2));

    % Get annual mean GMT
    dum = reshape(Surf_Temp{i},nens,12,[]);
    GMT = squeeze(mean(dum(:,mos_GMT,:),2));
        
    % For each ensemble member in this model set
    for j = 1:size(MIZ_ext,1)
        % Get linear fits over 1979-2018 using linear regression.
        % For MIZ, SIE, SIZ, PIZ.
        
        % Look at change in MIZ
        yval = SIE(j,:)'/1e6;
        yval = (yval - yval(1,:))/yval(1,:);
        [b,~] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_SIE_CLIVAR{i}(j) = slopemult*b(2);
        
        % Look at change in MIZ
        yval = MIZ_ext(j,:)'/1e6;
         yval = (yval - yval(1,:))/yval(1,:);
        [b,~] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_MIZ_CLIVR{i}(j) = slopemult*b(2);
        
        % Look at change in seasonal ice zone
        yval = SIZ_ext{i}(j,:)'/1e6;
        yval = (yval - yval(1,:))/yval(1,:);
        [b,~] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_SIZ_CLIVAR{i}(j) = slopemult*b(2);
        
        % Look at change in perennial ice zone
        yval = PIZ_ext{i}(j,:)'/1e6;
        yval = (yval - yval(1,:))/yval(1,:);
        [b,~] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_PIZ_CLIVAR{i}(j) = slopemult*b(2);
        
        % Look at change in GMT
        yval = GMT(j,:)';
%        yval = (yval - yval(1,:))/yval(1,:);
        [b,~] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
        slope_GMT_CLIVAR{i}(j) = slopemult*b(2);
        
    end
    
    %% Add each ensemble iteratively to plot to the
    subplot(121)
    hold on
    plotens(plot_yrs,SIE/1e6,clabs(i,:));
    
    
    subplot(122)
    hold on
    plotens(plot_yrs,MIZ_ext/1e6,clabs(i,:));
    %       ylim([0 4]);
    
    
    
end

%% Add observations 

for i = 1:length(PIZ_obs)
    
    SIE_obs(:,i) = squeeze(mean(sum(MIZ_band_obs{i}(bands_SIE,mos,:),1),2))/1e6;
    MIZ_obs(:,i) = squeeze(mean(sum(MIZ_band_obs{i}(bands_MIZ,mos,:),1),2))/1e6;

    % Look at SIE over time
    yval = squeeze(sum(mean(MIZ_band_obs{i}(bands_SIE,mos,:),2),1))/1e6;
    yval = (yval - yval(1,:))/yval(1,:);
    [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs_SIE(i,:) = slopemult*(bint(2,:) - b(2));
    slope_obs_SIE(i) = slopemult*b(2);
    
    % Look at SIE over time
    yval = squeeze(sum(mean(MIZ_band_obs{i}(bands_MIZ,mos,:),2),1))/1e6;
    yval = (yval - yval(1,:))/yval(1,:);
    [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs_MIZ(i,:) = slopemult*(bint(2,:) - b(2));
    slope_obs_MIZ(i) = slopemult*b(2);
    
end

for i = 1:length(GMT_obs)
    
    yval = GMT_obs{i}; % Should be years from 1979-2018
    [b, bint] = regress(yval,[ones(numel(plot_yrs),1) plot_yrs(:)-1979 ]);
    err_obs_GMT(i,:) = slopemult*(bint(2,:) - b(2));
    slope_obs_GMT(i) = slopemult*b(2);

end

%%


Ax{1} = subplot(221);

lspecs = {'-','--','-.'};

for i = 1:length(PIZ_obs)
    
    plot(plot_yrs,SIE_obs(:,i),[lspecs{i} 'k'],'linewidth',2,'Tag','Ens_Mean');
    
end
xlim([1979 2018]);
ylabel('10$^6$ sq km','interpreter','latex');
grid on; box on; set(gca,'xminorgrid','on'); 
title('Sea Ice Extent','interpreter','latex');
ylim([0 15]);

% Add ensemble indicators to ens mean plots
p = findall(gca,'Tag','Ens_Mean');

% legend(gca,p(1:length(ensemble_names)+2),strrep(fliplr([ensemble_names,obs_names]),'_','-'),'location','southwest');



Ax{2} = subplot(222);

for i = 1:length(PIZ_obs)
    
    plot(plot_yrs,MIZ_obs(:,i),[lspecs{i} 'k'],'linewidth',2,'Tag','Ens_Mean');
    
end

xlim([1979 2018]);
ylabel('10$^6$ sq km','interpreter','latex');
grid on; box on; set(gca,'xminorgrid','on'); 
title('Marginal Ice Zone','interpreter','latex');
% ylim([0 2.5]);


%% Now compare trends in each ensemble set

Ax{3} = subplot(223);

cla
gvec = vertcat(namevec{:}); % For labeling box plot
gvec = strrep(gvec,'_','-');
SIE_slopes = cell2mat(slope_SIE_CLIVAR);
MIZ_slopes = cell2mat(slope_MIZ_CLIVR);
SIZ_slopes = cell2mat(slope_SIZ_CLIVAR);
PIZ_slopes = cell2mat(slope_PIZ_CLIVAR);
GMT_slopes = cell2mat(slope_GMT_CLIVAR);

subplot(223)

gscatter(GMT_slopes,SIE_slopes,gvec);
grid on; box on; 
ylabel('SIE Trend'); xlabel('GMT Trend'); 
hold on

leg = legend;

leg.String{length(ensemble_names)+1} = GMT_names{1};

set(leg,'AutoUpdate','off')
legend('off'); 

for i = 1:length(PIZ_obs)
    
    for j = 1:length(GMT_obs)
        
    errorbar(slope_obs_GMT(j),slope_obs_SIE(i),err_obs_SIE(i,1),err_obs_SIE(i,2),err_obs_GMT(j,1),err_obs_GMT(j,2),'k','linewidth',1)
    
    end
    
end

%%
Ax{4} = subplot(224);

gscatter(GMT_slopes,MIZ_slopes,gvec);
grid on; box on; 
ylabel('MIZ Trend'); xlabel('GMT Trend'); 
hold on

leg = legend;

leg.String{length(ensemble_names)+1} = GMT_names{1};

set(leg,'AutoUpdate','off')
legend('off'); 

for i = 1:length(PIZ_obs)
    
    for j = 1:length(GMT_obs)
        
    errorbar(slope_obs_GMT(j),slope_obs_MIZ(i),err_obs_MIZ(i,1),err_obs_MIZ(i,2),err_obs_GMT(j,1),err_obs_GMT(j,2),'k','linewidth',1)
    
    end
    
end

%%

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',10,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    
end

pos = [10 6];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');



