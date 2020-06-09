function fh = Plot_CompareBehaviorStates_2States(motif_activity,yvals)

rng('default');
y=[];
for i = 1:size(motif_activity,1) %normalize and average
   temp = cat(1,motif_activity{i,:})';      
   temp = temp/nanmean(temp(:)); 
   y(:,i) = nanmean(temp);
end

figure('Position',[0 0 1000 5000]); hold on; 
s1 = subplot(312,'Units','centimeters','Position',[8 6 10 2.0]); hold on

axes(s1);
imagesc(y,yvals)
colormap(gca,flipud(redgreencmap(256,'Interpolation','quadratic')));
c = colorbar;
set(c,'YTick',[yvals(1), 1, yvals(2)]);
ylabel(c,{'Normalized Motif';'Intensity'},'FontSize',16,'Fontweight','normal','FontName','Arial');
set(c,'units','centimeters','position',[18.25 6 0.5 2])

for x_grid = 0.5:1:size(y,2)+0.5
    line([x_grid,x_grid],[0.5,size(y,1)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid = 0.5:1:size(y,1)+0.5
    line([0.5,size(y,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
xlabel('Basis Motifs')

xlim([0.5 size(y,2)+.5])
ylim([0.5 size(y,1)+0.5])
set(gca,'YColor','k')
setFigureDefaults

set(gcf,'Position',[680   150  875  650]);

fh = gcf;
end

