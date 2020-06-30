
load('CLIVAR_MIZ.mat');

OUT.SIE = [];
OUT.MIZ = [];
OUT.GMT = [];
OUT.namevec = {};
OUT.plotyrs = 1850:2100; %
OUT.names = ensemble_names;

for months = 1:12
    
    % Look at each ensemble member
    for i = 1:length(ensemble_names)
        
        if ~isempty(MIZ_band{i})
            
            % Number of ensemble members for this grouping
            nens = size(MIZ_band{i},1);
            % Total years in the grouping
            
            totyears = 251;
            SIE = nan(nens,12,totyears);
            MIZ = nan(nens,12,totyears);
            GMT = nan(nens,12,totyears);
            
            nyears = size(MIZ_band{i},3)/12;
            
            yearback = min(totyears,nyears);
            mosback = 12*yearback;
            
            
            % Get average sea ice extent for certain months over time.
            MIZ(:,:,end-yearback+1:end) = reshape(sum(MIZ_band{i}(:,bands_MIZ,end-mosback+1:end),2),nens,12,[]);
            MIZ = reshape(mean(MIZ(:,months,:),2),nens,totyears);
            
            % Add along first dimension
            OUT.MIZ = cat(1,OUT.MIZ,MIZ);
            
        end
        
    end
    
end

MIZ = OUT.MIZ;
MIZ = reshape(MIZ,[],12,251);
MIZ = MIZ(:,:,130:165);

for i = 1:size(MIZ,1)
    for j = 1:size(MIZ,2)
        dum = smoothdata(squeeze(MIZ(i,j,:)),'loess',10);
        dum = (dum - dum(1))/dum(1);
        
        [b,~] = regress(dum,[ones(36,1) [1:36]']);
        slope(i,j) = b(2);
    end
end

MIZ = 

