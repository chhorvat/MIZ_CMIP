clear

%% First do this for the climate model ensembles

ensembles_loc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/Model_Conc_Data/CLIVAR_LE/';


ensemble_names = {'GFDL_ESM2M','CANESM2','GFDL_CM3','CESM','CSIRO_MK36'};

low_bands = [0 .15 .3 .7 .8 .9];
high_bands = [.15 .3 .7 .8 .9 1];


%%
[MIZ_band,MIZA_band,PIZ_ext,SIZ_ext,Surf_Temp] = deal(cell(length(ensemble_names),1));

%%
for ensind = 1:length(ensemble_names)
    
    disp(['Ensemble ' ensemble_names{ensind}]);
    
    ensfold = [ensembles_loc ensemble_names{ensind} '/'];
    
    load([ensfold ensemble_names{ensind} '_grid.mat'],'lat','lon','grid_area','grid_area_atmos');
    
    conclist = dir([ensfold 'SIC/*.nc']);
    
    tslist = dir([ensfold 'TS/*.nc']);
    
    %%
    
    for i = 1:length(conclist)
        %%
        conc = ncread([ensfold 'SIC/' conclist(i).name],'sic');
        
        % Take the maximum backwards number of years possible for each.
        % Recognizing that they all end in 2100
        
        ntick = size(MIZ_band{ensind},3);
        
        if ntick > 1
            conc = conc(:,:,end-ntick+1:end);
        end
        
        if nanmax(conc(:)) > 2
            conc = conc / 100;
        end
        
        for j = 1:length(low_bands)
            
            usable = (lat > 0).*(conc > low_bands(j)).*(conc <= high_bands(j));
            usable(usable == 0) = nan;
            
            MIZ_band{ensind}(i,j,:) = squeeze(nansum(nansum(bsxfun(@times,usable,grid_area),1),2));
            MIZA_band{ensind}(i,j,:) = squeeze(nansum(nansum(bsxfun(@times,bsxfun(@times,usable,conc),grid_area),1),2));
            
        end
        
        hasice = reshape(conc,size(conc,1),size(conc,2),12,[]) > .15;
        
        someice = squeeze(1 * (lat > 0) .* (sum(hasice,3) > 0));
        allice = squeeze(1 * (lat > 0) .* (sum(hasice,3) == 12));
        
        SIZ_ext{ensind}(i,:) = squeeze(nansum(nansum(bsxfun(@times,someice,grid_area),1),2));
        PIZ_ext{ensind}(i,:) = squeeze(nansum(nansum(bsxfun(@times,allice,grid_area),1),2));
        
        
    end
    
    for i = 1:length(tslist)
        
        ts = ncread([ensfold 'TS/' tslist(i).name],'ts');
        ts = ts(:,:,end-ntick+1:end);
        
        Surf_Temp{ensind}(i,:) = sum(sum(bsxfun(@times,ts,grid_area_atmos),1),2) / sum(sum(grid_area_atmos));
        
    end
    
    
end

%%
save('CLIVAR_MIZ.mat','ensemble_names','PIZ*','MIZ*','SIZ*','*_bands','Surf_Temp');




