clearvars -except *bands

low_bands = [0 .15 .3 .7 .8 .9];
high_bands = [.15 .3 .7 .8 .9 1];

fileloc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/SSMI/';

obs_names = {'Bootstrap','NASATEAM','OSISAF'};

% Each dataset contains years from 1979-2019. 
% OSISAF 450 only covers through 2015. 
% OSISAF 450-b covers 2016-2019. 

[PIZ_obs,SIZ_obs,MIZ_band_obs,MIZA_band_obs] = deal(cell(length(obs_names),1));

for sind = 1:length(obs_names)
    
    load([fileloc obs_names{sind} '_monthly.mat'],'*_NH')
    
    
    conc = SSMI_monthly_NH;
    
    % Cover data gaps
    % Take average of prior and later year for each month
    
    % Cover for the lack of data on Dec 1987 and Jan 1988 in
    % Boostrap/NASAteam
    if sind < 3
    
        conc(:,:,12,9) = mean(conc(:,:,12,[8 10]),4);
        conc(:,:,1,10) = mean(conc(:,:,1,[9 11]),4);
    
    end
    
    if sind == 3
    % OSISAF loses data in April-May 1986
        conc(:,:,4:5,8) = mean(conc(:,:,4:5,[7 9]),4); 
    end
    
    MIZ = zeros(length(low_bands),size(SSMI_monthly_NH,3),size(SSMI_monthly_NH,4));
    lat = lat_NH;
    lon = lon_NH;
    grid_area = area_NH;
    
    
    for i = 1:length(low_bands)
        
        usable = (lat > 0).*(conc > low_bands(i)).*(conc <= high_bands(i));
        usable(usable == 0) = nan;

        MIZ(i,:,:) = squeeze(nansum(nansum(bsxfun(@times,usable,grid_area),1),2));
        MIZA(i,:,:) = squeeze(nansum(nansum(bsxfun(@times,bsxfun(@times,usable,conc),grid_area),1),2));
        
        hasice = reshape(conc,size(conc,1),size(conc,2),12,[]) > .15;
        
        someice = squeeze(1 * (lat > 0) .* (sum(hasice,3) > 0));
        allice = squeeze(1 * (lat > 0) .* (sum(hasice,3) == 12));
        
        SIZ_ext = squeeze(sum(sum(bsxfun(@times,someice,grid_area),1),2));
        PIZ_ext = squeeze(sum(sum(bsxfun(@times,allice,grid_area),1),2));
        
    end
    
    MIZA_band_obs{sind} = MIZA;
    MIZ_band_obs{sind} = MIZ;
    SIZ_obs{sind} = SIZ_ext;
    PIZ_obs{sind} = PIZ_ext;
    
end

save('OBS_MIZ.mat','*_obs','obs_names')
