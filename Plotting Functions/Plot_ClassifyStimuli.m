function stats = Plot_ClassifyStimuli(pixel_wise,motif_wise,pos)

num_mice = numel(pixel_wise);
rng('default');
figure('position',[ 680   258   440   720]); hold on;
x = ((rand(1,num_mice)/4)-0.125); %x placement randomly distributed around 1

%line at chance
line([-1 2],[0.5 0.5],'color',[0.6, 0.6 0.6],'linestyle','--','linewidth',2);

%plot lines between them
for i = 1:numel(x)
	plot([x(i),x(i)+1],[pixel_wise(i),motif_wise(i)],'color',[0.8 0.8 0.8],'linewidth',1.5);
end

%plot pixelwise auc
plot(x,pixel_wise,'color',[0.5 0.5 0.5],'linestyle','none','marker','.','markersize',20);
line([-0.1,0.1],[mean(pixel_wise),mean(pixel_wise)],'color','k','Linewidth',2);
line([0,0],[mean(pixel_wise)-std(pixel_wise)/sqrt(num_mice),mean(pixel_wise)+std(pixel_wise)/sqrt(num_mice)],'color','k','Linewidth',2)

%plot motifwise auc
plot(x+1,motif_wise,'color',[0.5 0.5 0.5],'linestyle','none','marker','.','markersize',20);
line([0.9,1.1],[mean(motif_wise),mean(motif_wise)],'color','k','Linewidth',2);
line([1,1],[mean(motif_wise)-std(motif_wise)/sqrt(num_mice),mean(motif_wise)+std(motif_wise)/sqrt(num_mice)],'color','k','Linewidth',2)


%significant testing
[h,p] = ttest(pixel_wise,motif_wise,'Tail','right');
AddSig(h,p,[0,1,0.35,0.35],2,0.02,1)
stats.pval_compare = p;

[h,p] = ttest(pixel_wise,0.5,'Tail','right');
AddSig(h,p,[0,0,0.7,0.7],2,0,1)
[h,p] = ttest([motif_wise],0.5,'Tail','right');
AddSig(h,p,[1,1,0.7,0.7],2,0,1)

ylim([0.2 0.7]);
ylabel({'Pixelwise Classification';'Accuracy (AUC)'})
set(gca,'Xlim',[-.5 1.5],'Xtick',[0,1],'XTickLabel',{'Original','Reconstruction'},'XTickLabelRotation',45,...
    'Ytick',(0.2:0.25:0.7),'clipping','on','units','centimeters')
set(gca,'position',pos)
title({'Original Data, But Not';'Motif Reconstruction,';'Distinguish Visual Stimuli'},'Fontsize',16,...
    'units','centimeters','Fontweight','normal','position',[2.05,9],'FontName','Arial')
setFigureDefaults;

end

