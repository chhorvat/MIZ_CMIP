
addpath('/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/DATA_TOOLS');


if ~exist('bands_SIE')
    
    bands_SIE = 2:6;
    bands_MIZ = 2:5;
    
end
%% Start with CLIVAR LENS

load('CLIVAR_MIZ.mat');

CLIVAR = struct();
CLIVAR.SIE = [];
CLIVAR.MIZ = [];
CLIVAR.GMT = [];
CLIVAR.SIA = [];
CLIVAR.CIA = [];
CLIVAR.MIZA = [];
CLIVAR.MIZA_F = [];
CLIVAR.MIZ_F = [];

CLIVAR.namevec = {};
CLIVAR.plotyrs = 1850:2100; %
CLIVAR.names = ensemble_names;

% Look at each ensemble member
for i = 1:length(ensemble_names)
    
    if ~isempty(MIZ_band{i})
        
        % Number of ensemble members for this grouping
        nens = size(MIZ_band{i},1);
        
        % Total years in the grouping
        
        totyears = 251;
        SIE = nan(nens,12,totyears);
        SIA = nan(nens,12,totyears);
        MIZA = nan(nens,12,totyears);
        MIZ = nan(nens,12,totyears);
        GMT = nan(nens,12,totyears);
        CIA = nan(nens,12,totyears);
        
        nyears = size(MIZ_band{i},3)/12;
        
        yearback = min(totyears,nyears);
        mosback = 12*yearback;
        
        % For grouping later
        CLIVAR.namevec = cat(1,CLIVAR.namevec,repmat({ensemble_names{i}},[nens,1]));
        
        % Get average sea ice extent for certain months over time.
        CIA(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(:,end,end-mosback+1:end),2),nens,12,[])/1e6;
        CIA = reshape(CIA,nens,12,totyears);
        
        % Get average sea ice extent for certain months over time.
        SIE(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(:,bands_SIE,end-mosback+1:end),2),nens,12,[])/1e6;
        SIE = reshape(SIE,nens,12,totyears);
        
        SIA(:,:,end-yearback+1:end) = reshape(sum(MIZA_band{i}(:,[1 bands_SIE],end-mosback+1:end),2),nens,12,[])/1e6;
        SIA = reshape(SIA,nens,12,totyears);
        
        % Get average sea ice extent for certain months over time.
        CIA(:,:,end-yearback+1:end) = reshape(sum(MIZA_band{i}(:,end,end-mosback+1:end),2),nens,12,[])/1e6;
        CIA = reshape(CIA,nens,12,totyears);
        
        % Get average sea ice extent for certain months over time.
        MIZ(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(:,bands_MIZ,end-mosback+1:end),2),nens,12,[])/1e6;
        MIZ = reshape(MIZ,nens,12,totyears);
        
        MIZ_F = 100*MIZ./SIE;
        
        % Get average sea ice extent for certain months over time.
        MIZA(:,:,end-yearback+1:end) = reshape(sum(MIZA_band{i}(:,[1 bands_MIZ],end-mosback+1:end),2),nens,12,[])/1e6;
        MIZA = reshape(MIZA,nens,12,totyears);
        
        MIZA_F = 100*MIZA./SIA;
        
        % Add along first dimension
        CLIVAR.SIA = cat(1,CLIVAR.SIA,SIA);
        CLIVAR.SIA_ensmean(i,:,:) = nanmean(SIA,1);
        CLIVAR.SIA_ensstd(i,:,:) = stdcorr(SIA,1,nens);
        
        % Add along first dimension
        CLIVAR.CIA = cat(1,CLIVAR.CIA,CIA);
        CLIVAR.CIA_ensmean(i,:,:) = nanmean(CIA,1);
        CLIVAR.CIA_ensstd(i,:,:) = stdcorr(CIA,1,nens);
        
        
        % Add along first dimension
        CLIVAR.SIE = cat(1,CLIVAR.SIE,SIE);
        CLIVAR.SIE_ensmean(i,:,:) = nanmean(SIE,1);
        CLIVAR.SIE_ensstd(i,:,:) = stdcorr(SIE,1,nens);
        
        
        % Add along first dimension
        CLIVAR.MIZ = cat(1,CLIVAR.MIZ,MIZ);
        CLIVAR.MIZ_ensmean(i,:,:) = nanmean(MIZ,1);
        CLIVAR.MIZ_ensstd(i,:,:) = stdcorr(MIZ,1,nens);
        
        % Add along first dimension
        CLIVAR.MIZA = cat(1,CLIVAR.MIZA,MIZA);
        CLIVAR.MIZA_ensmean(i,:,:) = nanmean(MIZA,1);
        CLIVAR.MIZA_ensstd(i,:,:) = stdcorr(MIZA,1,nens);
        
        % Add along first dimension
        CLIVAR.MIZ_F = cat(1,CLIVAR.MIZ_F,MIZ_F);
        CLIVAR.MIZ_F_ensmean(i,:,:) = nanmean(MIZ_F,1);
        CLIVAR.MIZ_F_ensstd(i,:,:) = stdcorr(MIZ_F,1,nens);
        
        % Add along first dimension
        CLIVAR.MIZA_F = cat(1,CLIVAR.MIZA_F,MIZA_F);
        CLIVAR.MIZA_F_ensmean(i,:,:) = nanmean(MIZA_F,1);
        CLIVAR.MIZA_F_ensstd(i,:,:) = stdcorr(MIZA_F,1,nens);
        
        % Get Global mean temperature
        GMT(:,:,end-yearback+1:end) = reshape(Surf_Temp{i}(:,end-mosback+1:end),nens,12,[]);
        GMT = reshape(GMT,nens,12,totyears);
        
        CLIVAR.GMT = cat(1,CLIVAR.GMT,GMT);
        CLIVAR.GMT_ensmean(i,:,:) = nanmean(GMT,1);
        CLIVAR.GMT_ensstd(i,:,:) = stdcorr(GMT,1,nens);
        
        
    end
    
end

%% Now add CMIP6 - these only go up to 2014

clearvars -except *CLIVAR mos* bands*

load('CMIP_MIZ.mat');

CMIP = struct();
CMIP.SIE = [];
CMIP.SIA = [];
CMIP.CIA = [];

CMIP.MIZA = [];
CMIP.MIZA_F = [];
CMIP.MIZ = [];
CMIP.MIZ_F = [];
CMIP.GMT = [];
[CMIP.GMT_ensmean,CMIP.GMT_ensstd] = deal(nan(length(ensemble_names),1,251)); 

[CMIP.SIE_ensmean,CMIP.SIE_ensstd,CMIP.SIA_ensmean,CMIP.SIA_ensstd,CMIP.CIA_ensmean,CMIP.CIA_ensstd, ...
    CMIP.MIZA_ensmean,CMIP.MIZA_ensstd,CMIP.MIZ_ensmean,CMIP.MIZ_ensstd, ...
    CMIP.MIZA_F_ensmean,CMIP.MIZA_F_ensstd,CMIP.MIZ_F_ensmean,CMIP.MIZ_F_ensstd] ... 
    = deal(nan(length(ensemble_names),12,251)); 

CMIP.namevec = {};
CMIP.plotyrs = 1850:2100; %
CMIP.names = ensemble_names;

% Look at each ensemble member
for i = 1:length(ensemble_names)
    
    if ~isempty(MIZ_band{i})
        
        usable = sum(Surf_Temp{i}(:,165),2) > 0;
        
        % Number of ensemble members for this grouping
        nens = sum(usable);
        
        if nens > 9
            enough_ens = 1;
        else
            enough_ens = nan;
        end
        
        % Total years in the grouping
        
        totyears = 251;
        SIE = nan(nens,12,totyears);
        SIA = nan(nens,12,totyears);
        CIA = nan(nens,12,totyears);
        MIZA = nan(nens,12,totyears);
        MIZ = nan(nens,12,totyears);
        GMT = nan(nens,totyears);
        
        nyears = size(MIZ_band{i},3)/12;
        
        yearback = min(totyears,nyears);
        mosback = 12*yearback;
        
        % Get average sea ice extent for certain months over time.
        SIE(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(usable,bands_SIE,end-mosback+1:end),2),nens,12,[])/1e6;
        SIE = reshape(SIE,nens,12,totyears);
        SIE(SIE == 0) = nan;
        
        % Get average sea ice extent for certain months over time.
        SIA(:,:,end-yearback+1:end) = reshape(sum(MIZA_band{i}(usable,[1 bands_SIE],end-mosback+1:end),2),nens,12,[])/1e6;
        SIA = reshape(SIA,nens,12,totyears);
        SIA(SIA == 0) = nan;
        
        % Get average sea ice extent for certain months over time.
        CIA(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(usable,end,end-mosback+1:end),2),nens,12,[])/1e6;
        CIA = reshape(CIA,nens,12,totyears);
        
        % Get average sea ice extent for certain months over time.
        MIZA(:,:,end-yearback+1:end) = reshape(sum(MIZA_band{i}(usable,[1 bands_MIZ],end-mosback+1:end),2),nens,12,[])/1e6;
        MIZA = reshape(MIZA,nens,12,totyears);
        MIZA(MIZA == 0) = nan;
        
        MIZA_F = 100*MIZA ./ SIA;
        
        % Get average sea ice extent for certain months over time.
        MIZ(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(usable,bands_MIZ,end-mosback+1:end),2),nens,12,[])/1e6;
        MIZ = reshape(MIZ,nens,12,totyears);
        MIZ(MIZ == 0) = nan;
        
        MIZ_F = 100*MIZ ./ SIE;
        
        % Get Global mean temperature
        GMT(:,end-yearback+1:end) = reshape(Surf_Temp{i}(usable,end-yearback+1:end),nens,[]);
        GMT(GMT==0) = nan;
        
        
        %% Add to existing record
        ens_yrly = size(Surf_Temp{i},1) - sum(isnan(Surf_Temp{i}),1);
        
        % For grouping later
        CMIP.namevec = cat(1,CMIP.namevec,repmat({ensemble_names{i}},[nens,1]));
        
        % Add along first dimension
        CMIP.MIZ = cat(1,CMIP.MIZ,MIZ);
        CMIP.MIZ_ensmean(i,:,:) = enough_ens*nanmean(MIZ,1);
        CMIP.MIZ_ensstd(i,:,:) = stdcorr(MIZ,1,ens_yrly);
        
        % Add along first dimension
        CMIP.SIA = cat(1,CMIP.SIA,SIA);
        CMIP.SIA_ensmean(i,:,:) = enough_ens*nanmean(SIA,1);
        CMIP.SIA_ensstd(i,:,:) = stdcorr(SIA,1,ens_yrly);
        
        % Add along first dimension
        CMIP.CIA = cat(1,CMIP.CIA,CIA);
        CMIP.CIA_ensmean(i,:,:) = enough_ens*nanmean(CIA,1);
        CMIP.CIA_ensstd(i,:,:) = stdcorr(CIA,1,ens_yrly);
        
        % Add along first dimension
        CMIP.SIE = cat(1,CMIP.SIE,SIE);
        CMIP.SIE_ensmean(i,:,:) = enough_ens*nanmean(SIE,1);
        CMIP.SIE_ensstd(i,:,:) = stdcorr(SIE,1,ens_yrly);
        
        % Add along first dimension
        CMIP.MIZ_F = cat(1,CMIP.MIZ_F,MIZ_F);
        CMIP.MIZ_F_ensmean(i,:,:) = enough_ens*nanmean(MIZ_F,1);
        CMIP.MIZ_F_ensstd(i,:,:) = stdcorr(MIZ_F,1,ens_yrly);
        
        % Add along first dimension
        CMIP.MIZA = cat(1,CMIP.MIZA,MIZA);
        CMIP.MIZA_ensmean(i,:,:) = enough_ens*nanmean(MIZA,1);
        CMIP.MIZA_ensstd(i,:,:) = stdcorr(MIZA,1,ens_yrly);
        
        % Add along first dimension
        CMIP.MIZA_F = cat(1,CMIP.MIZA_F,MIZA_F);
        CMIP.MIZA_F_ensmean(i,:,:) = enough_ens*nanmean(MIZA_F,1);
        CMIP.MIZA_F_ensstd(i,:,:) = stdcorr(MIZA_F,1,ens_yrly);
        
        CMIP.GMT = cat(1,CMIP.GMT,GMT);
        CMIP.GMT_ensmean(i,:,:) = enough_ens*nanmean(GMT,1);
        CMIP.GMT_ensstd(i,:,:) = stdcorr(GMT,1,ens_yrly);
        
        
        
    end
    
end




%% Now add observations

load('/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/SSMI/GMT_obs');
load('OBS_MIZ.mat');

OBS = struct();
OBS.SIE = [];
OBS.SIA = [];
OBS.MIZ = [];
OBS.MIZ_F = [];
OBS.MIZA = [];
OBS.MIZA_F = [];
OBS.GMT = [];
OBS.CIA = [];

OBS.names = obs_names;

OBS.plotyrs_SI = 1979:2018;
OBS.plotyrs_GMT = 1880:2019;
OBS.plotyrs = 1979:2018;

[~,common_yrs,~] = intersect(OBS.plotyrs_GMT,OBS.plotyrs_SI);

for i = 1:length(PIZ_obs)
    
    OBS.SIA(i,:,:) = squeeze(sum(MIZA_band_obs{i}([1 bands_SIE],:,:),1))/1e6;
    OBS.MIZA(i,:,:) = squeeze(sum(MIZA_band_obs{i}([1 bands_MIZ],:,:),1))/1e6;
    
    OBS.SIE(i,:,:) = squeeze(sum(MIZ_band_obs{i}(bands_SIE,:,:),1))/1e6;
    OBS.MIZ(i,:,:) = squeeze(sum(MIZ_band_obs{i}(bands_MIZ,:,:),1))/1e6;
    
    OBS.CIA(i,:,:) = squeeze(sum(MIZ_band_obs{i}(end,:,:),1))/1e6;
    
end

OBS.MIZA_F = 100*OBS.MIZA./OBS.SIA;
OBS.MIZ_F = 100*OBS.MIZ./OBS.SIE;

OBS.GMT = GMT_obs{1}(common_yrs);

%%
%% Now add FSD simulations

% load('/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/SSMI/GMT_obs');
% load('FSD_MIZ.mat');
% 
% FSD = struct();
% FSD.SIE = [];
% FSD.SIA = [];
% FSD.MIZ = [];
% FSD.MIZ_F = [];
% FSD.MIZA = [];
% FSD.MIZA_F = [];
% FSD.GMT = [];
% FSD.CIA = [];
% 
% FSD.names = names;
% 
% FSD.plotyrs = 1850:2100; %
% FSD.plotyrs_GMT = 1880:2019;
% 
% for i = 1:length(MIZ_band)
%     
%     FSD.SIA(i,:,:) = reshape(squeeze(sum(MIZA_band{i}([1 bands_SIE],:,:),1))/1e6,12,[]);
%     FSD.MIZA(i,:,:) = reshape(squeeze(sum(MIZA_band{i}([1 bands_MIZ],:,:),1))/1e6,12,[]);
%     
%     FSD.SIE(i,:,:) = reshape(squeeze(sum(MIZ_band{i}(bands_SIE,:,:),1))/1e6,12,[]);
%     FSD.MIZ(i,:,:) = reshape(squeeze(sum(MIZ_band{i}(bands_MIZ,:,:),1))/1e6,12,[]);
%     
%     FSD.CIA(i,:,:) = reshape(squeeze(sum(MIZ_band{i}(end,:,:),1))/1e6,12,[]);
%     
% end
% 
% FSD.MIZA_F = 100*FSD.MIZA./FSD.SIA;
% FSD.MIZ_F = 100*FSD.MIZ./FSD.SIE;
% 
% FSD.GMT = nan(length(FSD.plotyrs),1); 
% FSD.GMT(31:170) = GMT_obs{1};
% 
% 
