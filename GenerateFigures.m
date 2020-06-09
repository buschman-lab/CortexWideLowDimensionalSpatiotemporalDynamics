function GenerateFigures()

%Generates figures from:
%Camden J. MacDowell, Timothy J. Buschman,
%Low-Dimensional Spatiotemporal Dynamics Underlie Cortex-wide Neural Activity,
%Current Biology, 2020, ISSN 0960-9822,https://doi.org/10.1016/j.cub.2020.04.090.
%Contact: camdenm@princeton.edu for questions

%Add subdirectories
addpath(genpath(pwd));

%Load Figure Formating Parameters
fp = figure_parameters;

%% FIGURE 2
data = load('Fig_2_data.mat'); %load data
%data is a structure of measured variable x number of epochs (144)

%Plot Figure 2A, Number of Discovered Motifs
figure('position',[680   382   974   596]); hold on; 
histogram(data.number_motifs,'FaceColor',fp.c_discovery,'EdgeColor',fp.c_discovery-0.1,'linewidth',1,'BinWidth',1);    % Plot with errorbars
ylabel('Number of Epochs')
xlabel({'Number of';'Discovered Motifs'});
line([nanmedian(data.number_motifs),nanmedian(data.number_motifs)],[0 300],'color','k','linestyle','--','linewidth',1.5)
text(nanmedian(data.number_motifs)-12,19,['\mu_{1/2}=', num2str(round(nanmedian(data.number_motifs),2))],'FontSize',16,'FontWeight','normal','FontName','Arial');
xlim([0 30])
ylim([0 20])
title('Motifs Per Epoch','FontName','Arial','FontWeight','normal','FontSize',16);
setFigureDefaults
set(gca,'position',[2 3 3.5 8.5])
legend off

%Plot Figure 2B, Motif Percent Explained Variance
figure('position',[680   382   974   596]); hold on; 
histogram(data.PEV,'FaceColor',fp.c_discovery,'EdgeColor',fp.c_discovery-0.1,'linewidth',1,'BinWidth',2);    % Plot with errorbars
ylabel('Number of Epochs')
xlabel({'Percent Explained';'Variance'});
line([nanmedian(data.PEV),nanmedian(data.PEV)],[0 300],'color','k','linestyle','--','linewidth',1.5)
text(nanmedian(data.PEV)-35,25,['\mu_{1/2}=', num2str(round(nanmedian(data.PEV),2))],'FontSize',16,'FontWeight','normal','FontName','Arial');
xlim([50 100])
ylim([0 41])
title('Motif PEV','FontName','Arial','FontWeight','normal','FontSize',16);
setFigureDefaults
set(gca,'position',[2 3 3.5 8.5])
legend off

%Plot Figure 2C, Motif Occurence
occurance = [data.occurance{:}]/2; %convert to motifs/minute and combine across epochs 
figure('position',[680   382   974   596]); hold on; 
histogram(occurance,'FaceColor',fp.c_discovery,'EdgeColor',fp.c_discovery-0.1,'linewidth',1,'BinWidth',1);    % Plot with errorbars
ylabel('Number of Motifs')
xlabel({'Motif Frequency';'(occurence/minute)'});
line(nanmean(occurance),nanmean(occurance),[0 1000],'color','k','linestyle','--','linewidth',1.5)
text(nanmean(occurance)+3,300,['\mu=', num2str(round(nanmean(occurance),2))],'FontSize',16,'FontWeight','normal','FontName','Arial');
xlim([-1 20])
ylim([0 800])
title('Motif Frequency','FontName','Arial','FontWeight','normal','FontSize',16);
setFigureDefaults
set(gca,'position',[2 3 3.5 8.5])

%Plot Figure 2D, Motif Loadings
loadings = cumsum(data.loadings); %convert to cummulative sum
ci = NaN(size(loadings,1),2);
for i = 1:size(loadings,1)
    ci(i,:) = bootci(1000,@nanmedian,loadings(i,:));
end
figure('position',[1092 252 336 586]); hold on
errorbar(1:1:size(loadings,1),nanmedian(loadings,2),(nanmedian(loadings,2)-ci(:,1)),(ci(:,2)-nanmedian(loadings,2)),'LineWidth',2,'Color',[0.4 0.4 0.4])
title({'All Motifs Are Used to';'Explain Neural Activity'},'FontWeight','normal','units','centimeters','position',[2.5,9.25]);
xlabel({'Motif';'(Ordered By Decreasing';'Variance Per Epoch)'})
ylabel({'Relative Percent Explained';'Variance (Cummulative)'});
ylim([0 100]);
xticks((0:5:size(loadings,1)))
xlim([0 size(loadings,1)+0.5])
setFigureDefaults();
set(gca,'ytick',(0:20:100));
set(gca,'position',[3 4 5 8.5]);

%% FIGURE 3
data = load('Fig_3_data.mat'); %load data
%data is a structure of either measured variable x number of epochs (144)
%or number of motifs (2622) x measured variable

%Figure 3A
%compute the half life of the mean 
x = (1:1:(size(data.motif_autocorr,2))*1)';
y = nanmean(data.motif_autocorr,1)';
f = fit(x,y,'exp1','StartPoint',[0,0]);
figure('position',[542   420   336   555]); hold on;
shadedErrorBar(1:1:size(data.motif_autocorr,2),nanmean(data.motif_autocorr,1),(nanstd(data.motif_autocorr,[],1)),'lineprops',{'-','color',...
[0.4 0.4 0.4],'linewidth',2},'transparent',1,'patchSaturation',0.2);
p1 = plot(f,x,y,'k');
set(p1,'linewidth',2);
legend off;
title({'Correlation of Spatial';'Patterns in Motifs';'Decays Over Time'},...
    'FontWeight','normal','units','centimeters','position',[2.25,9]);
xlabel({'Temporal Offset';'Within Motif (ms)'})
ylabel('Correlation');
ylim([-0.2 0.62]);
xlim([0.5 size(data.motif_autocorr,2)+0.5])
ylim([-0.15 0.6])
setFigureDefaults();
set(gca,'XTickLabel',(1:2:13)*75,'XTick',1:2:13,'XTickLabelRotation',45)
set(gca,'Position',[2.25 3.05 5 8.5],'YColor','k')

%Figure 3B
figure('position',[680   514   242   464]); hold on
COL = {[0.4 0.4 0.4]};
CompareViolins(data.dissimilarity(:,2)',fp,'col',COL,'label','');
ylabel({'Dissimilarity (1-\rho)'})
setFigureDefaults();
set(gca,'position',[3 3 2 8.5],'ylim',[0 1])

%Figure 3D
figure('position',[680  200  242   600]); hold on
CompareViolins(data.static_vs_motifs,fp,'label',{'Motif','Static Network'},'col',{[0 0 0.75],[0.5 0.5 0.5]});
set(gca,'XTickLabelRotation',45)
ylim([0 100]);
ylabel({'Percent Explained Variance'})
setFigureDefaults;
set(gca,'position',[2,4,4,8.5])
[p, h] = signrank(data.static_vs_motifs(1,:),data.static_vs_motifs(2,:));
AddSig(h,p,[1,2,100,100],2,5,1)

%% Figure 4
data = load('Fig_4_data.mat'); %load data
%data is a structure of measured variable x number of epochs (144)

figure('position',[ 150   150   735   678]); hold on;
label = {'Within Mice','Between Mice','Basis Motifs','Static Networks'};
COL = {[0 0 0.75],[0.25 0.75 0],[0.85,0.37,0.001],[0.5 0.5 0.5]};
CompareViolins(data.PEV,fp,'label',label,'col',COL);
ylabel({'Percent Explained Variance'});
ylim([0 100])
setFigureDefaults();
set(gca,'XTickLabelRotation',45,'ytick',(0:20:100),'Clipping','off','position',[3 4 10 10])
[pval, h] = signrank(data.PEV(1,:),data.PEV(3,:));
AddSig(h,pval,[1 3 96.5],3,5,1); 
[pval, h] = signrank(data.PEV(1,:),data.PEV(2,:));
AddSig(h,pval,[1 2 24],3,6,1); 
[pval, h] = signrank(data.PEV(2,:),data.PEV(3,:));
AddSig(h,pval,[2 3 38],3,6,1); 
[pval, h] = signrank(data.PEV(3,:),data.PEV(4,:));
AddSig(h,pval,[3 4 12],3,6,1); 

%% Figure 5

%Fig 5B
data = load('Fig_5_data.mat');
motif_contribution = nanmean(data.motif_contribution);
ci = bootci(100,@nanmean,data.motif_contribution);
motif_id = 1:14;
[~,idx] = sort(motif_contribution,'ascend');
motif_id = motif_id(idx);
motif_contribution = motif_contribution(idx);
ci =ci(:,idx); 

figure; hold on;
line([0,0],[0 14],'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2)
plot(motif_contribution,1:numel(motif_contribution),'o','color','k','linewidth',2)
errorbar(motif_contribution,1:numel(motif_contribution),[],[],motif_contribution-ci(1,:),ci(2,:)-motif_contribution,'color','k','linewidth',1.5,'linestyle','none')
xlabel({'Contibution';'(AUC)'})
ylabel('Basis Motif (Sorted)');
set(gca,'ytick',1:numel(motif_contribution),'yticklabel',motif_id);
setFigureDefaults
% set(gca,'Units','centimeters','Position',[2 3 3 6]); 
set(gca,'Units','centimeters','Position',[3 3 4.5 5]);
set(gcf,'position',[680   400  722   600]);
xlim([-0.1 2])
ylim([0.5 14.5])
set(gca,'xtick',[0 0.5 1 1.5 2])
figure; hold on;
line([0,0],[0 14],'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2)
plot(motif_contribution,1:numel(motif_contribution),'o','color','k','linewidth',2)
errorbar(motif_contribution,1:numel(motif_contribution),[],[],motif_contribution-ci(1,:),ci(2,:)-motif_contribution,'color','k','linewidth',1.5,'linestyle','none')
set(gca,'ytick',1:numel(motif_contribution),'yticklabel',motif_id);
setFigureDefaults
xlim([10,13])
ylim([0.5 14.5])
set(gca,'Units','centimeters','Position',[2 3 1.2 5]); 
set(gcf,'position',[680   400  722   600]);



%Fig 5F
col = {[0.7 0.1 0],[0.5 0.5 0.5],[0.3 0.7 1],[0.5 0.5 0.5]};
%figure
figure('Position',[504   286   540  692]); hold on;
CompareViolins(data.PEV,fp,'col',col,'label',{'Basis: Visual','Avg: Visual','Basis: Tactile','Avg: Tactile'});
title({'Basis Motifs Capture Both Stimulus-';'Related and Ongoing Neural Activity'},'Fontsize',16,'Fontweight','normal','position',[2.5 108 0]);
setFigureDefaults;
set(gca,'XTickLabelRotation',45,'ytick',get(gca,'YTick'))
xlabel('Fit Type: Stimulus Type')
set(gca,'Clipping','off','LineWidth',2,'Fontsize',16,'Fontweight','normal');
set(gca,'Units','centimeters','position',[3 5 7 6]);
ylim([0 100])
[pval, h] = signrank(data.PEV(1,:),data.PEV(2,:));
AddSig(h,pval,[1 2 97],4,8,1); 
[pval, h] = signrank(data.PEV(3,:),data.PEV(4,:));
AddSig(h,pval,[3 4 97],4,8,1); 

%% Figure 6;
data = load('Fig_6_data.mat');

%Figure 6B
Plot_ClassifyStimuli(data.AUC(1,:),data.AUC(2,:),[3,4,5,6.5]);

%Figure 6C
figure('position',[ 680   258   440   720]); hold on;
x = ((rand(1,size(data.PEV,2))/4)-0.125)+(1:3)'; %x placement randomly distributed around 1
for i = 1:size(data.PEV,2) %loop through mice and plot
    plot(x(:,i),data.PEV(:,i),'color',[0.5 0.5 0.5],'linestyle','none','marker','.','markersize',20);
end
ylabel({'Percent Explainable Variance';'(During Stimulus Period)'})
z = std(data.PEV,[],2)/sqrt(size(data.PEV,2));
line([1:3; 1:3],[nanmean(data.PEV,2)-z,nanmean(data.PEV,2)+z]','color','k','linewidth',2);
line([(1:3)-0.1; (1:3)+0.1],[nanmean(data.PEV,2),nanmean(data.PEV,2)]','color','k','linewidth',2);
set(gca,'xlim',[0.5 3.5],'Ylim',[0 100],'FontSize',18,'FontWeight','normal','FontName','Arial',...
    'ytick',(0:25:100),'units','centimeters','position',[2.75 7 6 6.5],...
    'xticklabels',{'Non-Specific Motifs',sprintf('Stim-Specific Motif (#%d)',10),...
    'Pixelwise Stim-Specific Activity'},'XTickLabelRotation',45,'Linewidth',2)
setFigureDefaults;


%% Figure 7
data = load('Fig_7_data.mat');

%figure 7B
Plot_CompareBehaviorStates_2States(data.motif_activity,[0.9 1.1]);

%figure 7E
figure; hold on
COL = {[0.85,0.37,0.001],[0.85,0.37,0.001],[.14 0.4 0.6431],[0.14 0.4 0.6431]};
label = {'Solo','Social'};
CompareViolins(data.PEV,fp,'label',label,'col',COL);
title({'Basis Motifs Generalize Across';'Social and Solo Environments'},'FontWeight','normal','FontName','Arial',...
    'units','centimeters','Position',[1.75 9.25]);
ylabel({'Percent Explained Variance'});
ylim([0 100]);
setFigureDefaults;
set(gca,'XTickLabelRotation',45,'ytick',(0:20:100))
[pval, h] = ranksum(data.PEV(1,:),data.PEV(2,:)); %independent samples solo vs social environment
AddSig(h,pval,[1 2 99],2,5,1); 
set(gca,'Clipping','off','LineWidth',2,'Fontsize',16,'Fontweight','normal');
set(gca,'Units','centimeters','position',[4 5 4 8.5]);
set(gcf,'Position',[524   290   500  600]);

%figure 7F
figure(); hold on
%Plot the median of both
plot(0:.1:1,0:.1:1,'LineWidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--')
scatter(nanmedian(data.Relative_PEV{1},2),nanmedian(data.Relative_PEV{2},2),20,'filled','k')

%Add circles around the significant ones, were significance color
for i = 1:size(data.Relative_PEV{1},1)
    [pval, h] = ranksum(data.Relative_PEV{1}(i,:),data.Relative_PEV{2}(i,:));
    if h == 1                
        if pval<0.0001
            scatter(nanmedian(data.Relative_PEV{1}(i,:)),nanmedian(data.Relative_PEV{2}(i,:)),100,[1 0 0],'filled')
        elseif pval<0.05/size(data.Relative_PEV{1},1) %bonferonni corrected
            scatter(nanmedian(data.Relative_PEV{1}(i,:)),nanmedian(data.Relative_PEV{2}(i,:)),50,[1 0 0],'filled')          
        end
    end
    text([nanmedian(data.Relative_PEV{1}(i,:))-.005,nanmedian(data.Relative_PEV{1}(i,:))-.005],...
    [nanmedian(data.Relative_PEV{2}(i,:)),nanmedian(data.Relative_PEV{2}(i,:))],...
    sprintf('%d',i),'Fontsize',16,'FontWeight','normal',...
    'Color',[0.3 0.3 0.3],'HorizontalAlignment','center') 
    
end

title({'Basis Motifs Are Expressed Differently';'in Solo and Social Environments'},'FontWeight','normal');
xlabel({'Relative Percent Explained';'Variance Solo Environment'});
ylabel({'Relative Percent Explained';'Variance Social Environment'});
ylim([0 0.26]);
xlim([0 0.26]);
setFigureDefaults;
set(gca,'xtick',get(gca,'xtick'),'xticklabels',get(gca,'XTick')*100,'ytick',get(gca,'YTick'),'yticklabels', get(gca,'YTick')*100);
set(gca,'Position',[3 3 8.5 8.5]);
set(gcf,'Position',[652   272   502   518])



end %function end



































































    