function plotens(xax,field,color)

ens_mean = mean(field,1); 

ens_var = std(field,0,1); 

plot(xax,field,'linewidth',.25,'color',[217,217,217]/256); 

hold on

plot(xax,ens_mean,'linewidth',1,'color',color,'Tag','Ens_Mean','linewidth',1);

% plot(xax,ens_mean + ens_var,'--','color',color,'linewidth',1);
% plot(xax,ens_mean - ens_var,'--','color',color,'linewidth',1);





