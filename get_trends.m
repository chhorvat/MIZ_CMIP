%% Get Trends in CLIVAR

clearvars -except OBS CMIP CLIVAR bands* smooth* obuse


hist_yrs = 130:165; % 1979-2014
hist_yrs_obs = 1:36; % 1979-2014
fut_yrs = 171:206; % 2020-2055

slopemult = 1; % length(hist_yrs); % To per decade


yval = CLIVAR.SIA/1e6;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs(1)))/dum(hist_yrs(1));
        
        [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
        CLIVAR.SIA_slope_hist(i,j) = slopemult*b(2);
        
        [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
        CLIVAR.SIA_slope_fut(i,j) = slopemult*b(2);
        
    end
end

yval = CLIVAR.SIE/1e6;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs(1)))/dum(hist_yrs(1));
        
        [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
        CLIVAR.SIE_slope_hist(i,j) = slopemult*b(2);
        
        [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
        CLIVAR.SIE_slope_fut(i,j) = slopemult*b(2);
        
    end
end

yval = CLIVAR.MIZA_F;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs(1)))/dum(hist_yrs(1));
        
        [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
        CLIVAR.MIZA_F_slope_hist(i,j) = slopemult*b(2);
        
        [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
        CLIVAR.MIZA_F_slope_fut(i,j) = slopemult*b(2);
        
    end
end

yval = squeeze(mean(CLIVAR.GMT,2));

for i = 1:size(yval,1)
        
        dum = squeeze(yval(i,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        
        [b,~] = regress(dum(hist_yrs)',[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
        CLIVAR.GMT_slope_hist(i) = slopemult*b(2);
        
        [b,~] = regress(dum(fut_yrs)',[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
        CLIVAR.GMT_slope_fut(i) = slopemult*b(2);
        
end
%% Get Trends in CMIP

yval = CMIP.SIA_ensmean/1e6;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs(1)))/dum(hist_yrs(1));
        
        if ~isnan(sum(dum(hist_yrs)))
            
            [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
            CMIP.SIA_slope_hist(i,j) = slopemult*b(2);
            
        else
            
            CMIP.SIA_slope_hist(i,j) = nan;
            
        end
        
        if ~isnan(sum(dum(fut_yrs)))
            
            [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
            
            CMIP.SIA_slope_fut(i,j) = slopemult*b(2);
            
        else
            CMIP.SIA_slope_fut(i,j) = nan;
        end
        
    end
end

yval = CMIP.SIE_ensmean/1e6;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs(1)))/dum(hist_yrs(1));
        
        if ~isnan(sum(dum(hist_yrs)))
            
            [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
            CMIP.SIE_slope_hist(i,j) = slopemult*b(2);
            
        else
            
            CMIP.SIE_slope_hist(i,j) = nan;
            
        end
        
        if ~isnan(sum(dum(fut_yrs)))
            
            [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
            
            CMIP.SIE_slope_fut(i,j) = slopemult*b(2);
            
        else
            CMIP.SIE_slope_fut(i,j) = nan;
        end
        
    end
end

yval = CMIP.GMT_ensmean;

for i = 1:size(yval,1)
    
    dum = squeeze(yval(i,:))'; % smoothdata(yval(i,:),smooth_type,smooth_window)';
    
    if ~isnan(sum(dum(hist_yrs)))
        
        [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
        CMIP.GMT_slope_hist(i) = slopemult*b(2);
        
    else
        
        CMIP.GMT_slope_hist(i) = nan;
        
    end
    
    if ~isnan(sum(dum(fut_yrs)))
        
        [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
        
        CMIP.GMT_slope_fut(i) = slopemult*b(2);
        
    else
        CMIP.GMT_slope_fut(i) = nan;
    end
    
end

yval = CMIP.MIZA_F_ensmean;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs(1)))/dum(hist_yrs(1));
        
        if ~isnan(sum(dum(hist_yrs)))
            
            [b,~] = regress(dum(hist_yrs),[ones(length(hist_yrs),1) [hist_yrs-hist_yrs(1)]']);
            CMIP.MIZA_F_slope_hist(i,j) = slopemult*b(2);
            
        else
            
            CMIP.MIZA_F_slope_hist(i,j) = nan;
            
        end
        
        if ~isnan(sum(dum(fut_yrs)))
            
            [b,~] = regress(dum(fut_yrs),[ones(length(fut_yrs),1) [fut_yrs - fut_yrs(1)]']);
            
            CMIP.MIZA_F_slope_fut(i,j) = slopemult*b(2);
            
        else
            
            CMIP.MIZA_F_slope_fut(i,j) = nan;
            
        end
        
    end
    
end
%% Get Trends in OBS

yval = OBS.SIA/1e6;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs_obs(1)))/dum(hist_yrs_obs(1));
        
        [b,bint] = regress(dum(hist_yrs_obs),[ones(length(hist_yrs_obs),1) [hist_yrs_obs-hist_yrs_obs(1)]']);
        OBS.SIA_slope_hist(i,j) = slopemult*b(2);
        OBS.SIA_slope_err(i,j) = slopemult*(bint(2,1) - b(2));
        
    end
end

yval = OBS.SIE/1e6;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs_obs(1)))/dum(hist_yrs_obs(1));
        
        [b,bint] = regress(dum(hist_yrs_obs),[ones(length(hist_yrs_obs),1) [hist_yrs_obs-hist_yrs_obs(1)]']);
        OBS.SIE_slope_hist(i,j) = slopemult*b(2);
        OBS.SIE_slope_err(i,j) = slopemult*(bint(2,1) - b(2));
        
    end
end


yval = OBS.MIZA_F;

for i = 1:size(yval,1)
    for j = 1:size(yval,2)
        
        dum = squeeze(yval(i,j,:)); % smoothdata(yval(i,j,:),smooth_type,smooth_window)';
        dum = 100*(dum - dum(hist_yrs_obs(1)))/dum(hist_yrs_obs(1));
        
        [b,bint] = regress(dum(hist_yrs_obs),[ones(length(hist_yrs_obs),1) [hist_yrs_obs-hist_yrs_obs(1)]']);
        OBS.MIZA_F_slope_hist(i,j) = slopemult*b(2);
        OBS.MIZA_F_slope_err(i,j) = slopemult*(bint(2,1) - b(2));
        
    end
end

yval = OBS.GMT';

for i = 1:size(yval,1)
    
    dum = squeeze(yval(i,:))'; % smoothdata(yval(i,:),smooth_type,smooth_window)';
    
    [b,bint] = regress(dum(hist_yrs_obs),[ones(length(hist_yrs_obs),1) [hist_yrs_obs-hist_yrs_obs(1)]']);
    OBS.GMT_slope_hist(i) = slopemult*b(2);
    OBS.GMT_slope_err(i) = slopemult*(bint(2,1) - b(2));
    
end
