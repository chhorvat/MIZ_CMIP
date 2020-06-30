figure

movec = 'JFMAMJJASOND'; 

for i = 1:12
    
    Ax{i} = subplot(2,12,i);
    
    boxplot(sens_SIE(:,i),gvec); 
    hold on
    
    yline(mean(sens_SIE_obs(:,i)),'k'); 
    
    ylim([-1 1]); 

    grid on; 
    box on; 
        set(gca,'xticklabel',{}); 
    title(movec(i),'interpreter','latex'); 
end

for i = 1:12
    
    Ax{12+i} = subplot(2,12,12+i);
    
    boxplot(sens_MIZ(:,i),gvec); 
    hold on
    
    yline(mean(sens_MIZ_obs(:,i)),'k'); 
    
    ylim([-25 25]); 

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

saveas(gcf,'Sensitivities-bymo.pdf');
