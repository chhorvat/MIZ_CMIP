clear
close all

plot_yrs = 1979:2018;

bands_SIE = 2:3;
bands_MIZ = 2;
mos_GMT = [1:12];
mos_SIE = [9];
mos_MIZ = [1:5 10:12];

addpath('../'); 

plot_preamble;

smooth_window = 1;
smooth_type = 'loess';
obuse = [1 2 3]; % Which of the SIC obs to use

%%

get_trends; 
%%
horvat_colors; 

trend_scatters; 

%%
trend_boxes; 

%%
trend_phase; 

%% 


