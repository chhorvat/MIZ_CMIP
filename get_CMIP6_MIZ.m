clear

ensembles_loc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/Model_Conc_Data/CMIP6/';

load SIMIP_data

% ensemble_names = {'BCC-CSM2-MR','BCC-ESM1','CanESM5','CanESM5-CanOE', ...
%     'CESM2','CESM2-FV2','CESM2-WACCM','CESM2-WACCM-FV2','CIESM','INM-CM4-8', ...
%     };

% Issues with or missing runs from
% CanESM
% E3SM-1-0 r1p1f1 1975-1980

try
    
    load('CMIP_MIZ.mat')
    load('list_of_CMIP6_runs.mat');
    
    
catch err
    
    load('list_of_CMIP6_runs','ensemble_names');
    
    [MIZ_band,MIZA_band,PIZ_ext,SIZ_ext,Surf_Temp] = deal(cell(length(ensemble_names),1));
    
    ensemble_finished = 0*zeros(length(ensemble_names),1); % Whether all runs from a given model are done
    runlist = {}; % Names of each run we have done
    run_finished = 0 * length(SIMIP_names_hist); % Whether a given run is done
    
end

% E3SM-1-0 is all busted
% HadGEM3-GC31-MM also but because it wants too much memory

low_bands = [0 .15 .3 .7 .8 .9];
high_bands = [.15 .3 .7 .8 .9 1];

unfinished_ensembles = ensemble_names(~ensemble_finished);

%%
for todo_ind = 1:length(unfinished_ensembles)
    %%
    % Which model are we using out of unfinished ones?
    ens_ind = find(strcmp(ensemble_names,unfinished_ensembles{todo_ind}));
     
    %
        
    locs = find(contains(SIMIP_names_hist,[ensemble_names{ens_ind} '_']));
    
    nens = length(locs);
    % If we have actual ensemble members to look at in CMIP6 dataset
    
    if nens > 0
        
        doneflag = 1; % Reverts to zero if there is an error in loading
        
        disp(['Ensemble ' ensemble_names{ens_ind}]);
        
        memb_ct = 1;
                
        for member_ind = 1:nens
            
            % Folder containing all runs
            ensfold = [ensembles_loc ensemble_names{ens_ind} '/'];
            
            % Identify the run we want to consider first
            run_id = SIMIP_names_hist{locs(member_ind)}(length(ensemble_names{ens_ind})+2:end);
            
            % See whether we have the rcp8.5 run for GMT or not
            
            try
                
                % Load in grid
                load([ensfold ensemble_names{ens_ind} '_grid.mat'],'lat','lon','grid_area');
                
                % Find the right files
                conclist = dir([ensfold 'SIC/*' run_id '*.nc']);
                
                disp(['Run ID is ' run_id]);
                
                %% Now read all of the relevant concentrations in
                % The concentration file may be over several tiles
                
                conc = nan(size(grid_area,1),size(grid_area,2),3012);
                T = nan(3012,1);
                
                for run_sub = 1:length(conclist)
                    
                    yr1 = str2num(conclist(run_sub).name(end-15:end-12)); 
                    mo1 = str2num(conclist(run_sub).name(end-11:end-10)); 
                    yr2 = str2num(conclist(run_sub).name(end-8:end-5)); 
                    mo2 = str2num(conclist(run_sub).name(end-4:end-3)); 
                    
                    ind1 = (yr1-1850)*12 + mo1; 
                    ind2 = (yr2-1850)*12 + mo2; 
                    
                    if ind1 < 3012
                        
                        conc(:,:,ind1:ind2) = ncread([ensfold 'SIC/' conclist(run_sub).name],'siconc');
                        T(ind1:ind2) = ncread([ensfold 'SIC/' conclist(run_sub).name],'time');
                        
                    end
                    
                end
                
                if nanmax(conc(:)) > 2
                    conc = conc / 100;
                end
                
                
                %% Now process the full concentration record
                
                
                % Evaluate concentration bands in the obs
                for j = 1:length(low_bands)
                    
                    usable = (lat > 0).*(conc > low_bands(j)).*(conc <= high_bands(j));
                    usable(usable == 0) = nan;
                    
                    MIZ_band{ens_ind}(memb_ct,j,:) = squeeze(nansum(nansum(bsxfun(@times,usable,grid_area),1),2));
                    MIZA_band{ens_ind}(memb_ct,j,:) = squeeze(nansum(nansum(bsxfun(@times,bsxfun(@times,usable,conc),grid_area),1),2));
                    
                end
                
                % hasice is all that have ice at some point during the year
                hasice = reshape(conc,size(conc,1),size(conc,2),12,[]) > .15;
                
                % These are those that have some ice
                someice = squeeze(1 * (lat > 0) .* (sum(hasice,3) > 0));
                % These are those that have ice in all 12 months
                allice = squeeze(1 * (lat > 0) .* (sum(hasice,3) == 12));
                
                % Create SIZ and PIZ for this ensemble member and run
                SIZ_ext{ens_ind}(memb_ct,:) = squeeze(nansum(nansum(bsxfun(@times,someice,grid_area),1),2));
                PIZ_ext{ens_ind}(memb_ct,:) = squeeze(nansum(nansum(bsxfun(@times,allice,grid_area),1),2));
                
                %% Now pull in the CMIP6 GMT data
                
                % See if this particular run string is in the rcp8p5 set
                loc2 = find(strcmp(SIMIP_names_hist{locs(member_ind)},SIMIP_names_8p5));
                
                if ~isempty(loc2) 
                    
                    % If it is, take the gmst as the rcp8p5 data
                    Surf_Temp{ens_ind}(memb_ct,1:250) = SIMIP_GMT_8p5(:,loc2);
                    Surf_Temp{ens_ind}(memb_ct,251) = nan;
                    
                else
                    
                    % If not, set all years beyond 2014 as nan
                    Surf_Temp{ens_ind}(memb_ct,1:165) = SIMIP_GMT_hist(:,locs(member_ind));
                    Surf_Temp{ens_ind}(memb_ct,166:251) = nan;
                    
                end
                
                run_finished(locs(memb_ct)) = 1;
                
                runlist{end+1} = [ensemble_names{ens_ind} '_' run_id];
                
                
                
            catch err
                
                
                doneflag = 0;
                
                disp(['Failure on ' ensfold 'SIC/*' run_id '*.nc']);
                
                
            end
            
            memb_ct = memb_ct + 1;
            
        end
        
        % If there were no errors, add to the list of finished runs
        if doneflag
            
            % Isn't finished if there are no ensemble members!
            ensemble_finished(ens_ind) = 1;
            
        end
        
    end
    
    save('CMIP_MIZ.mat','ensemble_names*','ensemble_finished','runlist','run_finished','SIMIP_names_hist','Surf_Temp','PIZ*','MIZ*','SIZ*','*_bands');
    
end



