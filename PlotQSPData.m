% PlotQSPData: script for plotting the results obtained by running 
%   QSPClusterAnalysis.py and ParseQSPData.m
% Author: Arkadz Kirshtein, https://sites.google.com/site/akirshtein/
% (c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
%
% If using this or related code please cite 
% Kirshtein, A.; Akbarinejad, S.; Hao, W.; Le, T.; Aronow, R.A.; Shahriyari, L. 
%     Data driven mathematical model of colon cancer progression. 
%     (Manuscript submitted for publication).

%Plot sensitivity
paralabels={'\lambda_{T_hD}';  '\lambda_{T_hM}';  '\lambda_{T_h\mu_1}'; '\lambda_{T_CT_h}';'\lambda_{T_CD}'; '\lambda_{T_rT_h}'; '\lambda_{T_r\mu_2}';'\lambda_{T_rG_\beta}';   ...
      '\lambda_{DH}';   '\lambda_{DC}';     '\lambda_{M\mu_2}';   '\lambda_{MI_\gamma}'; '\lambda_{MT_h}'; '\lambda_{C}';    '\lambda_{C\mu_1}'; '\alpha_{NC}';   ...
      '\lambda_{\mu_1T_h}';'\lambda_{\mu_1M}';   '\lambda_{\mu_1D}';   '\lambda_{\mu_2M}';'\lambda_{\mu_2D}';'\lambda_{\mu_2T_r}';  ...
      '\lambda_{HN}';   '\lambda_{HM}';     '\lambda_{HT_h}';    '\lambda_{HT_C}'; '\lambda_{HT_r}'; '\lambda_{I_\gammaT_h}'; '\lambda_{I_\gammaT_C}'; '\lambda_{I_\gammaM}';   ...
      '\lambda_{G_\betaM}';  '\lambda_{G_\betaT_r}';   '\delta_{T_N}';      '\delta_{T_h\mu_2}';'\delta_{T_hT_r}'; '\delta_{T_h}';    '\delta_{T_C\mu_2}'; '\delta_{T_CT_r}'; ...
      '\delta_{T_C}';    '\delta_{T_r\mu_1}';   '\delta_{T_r}';      '\delta_{DH}';   '\delta_{DC}';   '\delta_{D}';     '\delta_{M}';     '\delta_{CG_\beta}';'\delta_{CI_\gamma}';'\delta_{CT_C}'; ...
      '\delta_{C}';     '\delta_{N}';       '\delta_{\mu_1}';     '\delta_{\mu_2}';  '\delta_{H}';     '\delta_{I_\gamma}';    '\delta_{G_\beta}'; 'A_{T_N}';'A_{Dn}';'M_0';'C_0'; '\alpha_{T_NT_h}'; '\alpha_{T_NT_C}'; '\alpha_{T_NT_r}'; '\alpha_{D_ND}'};
grid_level=3; %grid level for sensitivity
  
  
% Steady state primary
figure;
sensvals=cell(5,3);
for cluster=1:5
    Sensdat=dlmread(['Data/Sensitivity/V63-grid' num2str(grid_level) '-cluster-' num2str(cluster) '-of-5-results-sensitivity_steady.csv']);
    for i=1:3
        subplot(5,3,3*(cluster-1)+i)
        threshold=1;
        [~,index]=sort(abs(Sensdat(:,i)),'descend');
        index=index(1:4);
        sensvals{cluster,i}=Sensdat(index,i);
        bar(Sensdat(index,i))
        set(gca,'fontsize',9)
        xlim([0 length(index)+1])
        xticks(1:length(index))
        xticklabels(paralabels(index))
    end
end

subplot(5,3,1)
title('Sensitivity of Cancer')
subplot(5,3,2)
title('Sensitivity of Total Cell Count')
subplot(5,3,3)
title('Sensitivity of minimal eigenvalue')
for i=1:5
    subplot(5,3,3*i-2)
    ylabel(['Cluster ' num2str(i)])
    subplot(5,3,3*i-1)
    ylabel(['Cluster ' num2str(i)])
    subplot(5,3,3*i)
    ylabel(['Cluster ' num2str(i)])
end
sgtitle('Most sensitive parameters')
  
% Steady state secondary
figure;
subindex=setdiff(1:63,[14:16, 46:50, 59]);
sensvals=cell(5,2);
for cluster=1:5
    Sensdat=dlmread(['Data/Sensitivity/V63-grid' num2str(grid_level) '-cluster-' num2str(cluster) '-of-5-results-sensitivity_steady.csv']);
    for i=1:2
        subplot(5,2,2*(cluster-1)+i)
        threshold=1;
        [~,index]=sort(abs(Sensdat(subindex,i)),'descend');
        index=index(1:4);
        sensvals{cluster,i}=Sensdat(subindex(index),i);
        bar(Sensdat(subindex(index),i))
        set(gca,'fontsize',9)
        xlim([0 length(index)+1])
        xticks(1:length(index))
        xticklabels(paralabels(subindex(index)))
    end
end

subplot(5,2,1)
title('Sensitivity of Cancer')
subplot(5,2,2)
title('Sensitivity of Total Cell Count')
for i=1:5
    subplot(5,2,2*i-1)
    ylabel(['Cluster ' num2str(i)])
    subplot(5,2,2*i)
    ylabel(['Cluster ' num2str(i)])
end
sgtitle('Most sensitive immune parameters')

clear;

%Plot dynamics

% variability shaded region
varshade=true;

% reading parsed data
load('Data\Dynamic\V63-results-all-data.mat')

figure;
cmap=lines(5);
for i=1:9
    subplot(2,5,i);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,i+1)*cells(cluster,i),'color',cmap(cluster,:));
    end
    if varshade
        for cluster=1:5
                xf=[data{cluster}(:,1); flip(data{cluster}(:,1))];
                yf=[mindata{cluster}(:,i+1)*cells(cluster,i);
                    flip(maxdata{cluster}(:,i+1)*cells(cluster,i))];
                fill(xf,yf,cmap(cluster,:),'FaceAlpha',0.1,'EdgeAlpha',0);
        end
    end
    hold off
    xlabel('time (days)')
    ylabel(vars{i})
    xlim([0 T])
end
subplot(2,5,10);
hold on
for cluster=1:5
    plot(data{cluster}(:,1),data{cluster}(:,2:10)*cells(cluster,lmod)','color',cmap(cluster,:));
end
if varshade
    for cluster=1:5
        xf=[data{cluster}(:,1); flip(data{cluster}(:,1))];
        yf=[mindata{cluster}(:,2:10)*cells(cluster,lmod)';
            flip(maxdata{cluster}(:,2:10)*cells(cluster,lmod)')];
        fill(xf,yf,cmap(cluster,:),'FaceAlpha',0.1,'EdgeAlpha',0);
    end
end
hold off
xlabel('time (days)')
ylabel(vars{15})
xlim([0 T])
legend({'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5'})

%optional visible scale adjustment
subplot(2,5,2);ylim([0 7000])
subplot(2,5,3);ylim([0 7000])
subplot(2,5,6);ylim([0 700])

sgtitle('Cell dynamics')

%cytokines
figure;
cmap=lines(5);
for i=10:14
    subplot(1,5,i-9);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,i+1)*cells(cluster,i),'color',cmap(cluster,:));
    end
    if varshade
        for cluster=1:5
                xf=[data{cluster}(:,1); flip(data{cluster}(:,1))];
                yf=[mindata{cluster}(:,i+1)*cells(cluster,i);
                    flip(maxdata{cluster}(:,i+1)*cells(cluster,i))];
                fill(xf,yf,cmap(cluster,:),'FaceAlpha',0.1,'EdgeAlpha',0);
        end
    end
    hold off
    xlabel('time (days)')
    ylabel(vars{i})
    xlim([0 T])
end
legend({'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5'})

%optional visible scale adjustment
subplot(1,5,2);ylim([0 1600])
subplot(1,5,3);ylim([0 11000])
subplot(1,5,4);ylim([0 12])
subplot(1,5,5);ylim([0 35000])
sgtitle('Cytokine dynamics')

clear;
% Plot variable initial conditions
vars={'Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TGF-\beta', 'Total cells'};
tab=readtable('input/Large_Tumor_variable_scaling.csv');
cells=tab{:,:};
lmod=1:9;
plot_id=[lmod, 11:15];

load('Data\InitialDifference\V63-patients-all-data.mat');

T=3500;

clusters=5;

cellscaling=cells(:,lmod);

cmap=cell(1,clusters);
figure_handles=cell(1,clusters);
for c=1:clusters
    figure;
    figure_handles{c}=gcf;
    for k=1:15
        subplot(3,5,k)
        hold on
    end
    cmap{c}=lines(sum(cluster_id==c));
end


for p=1:patients
    figure(figure_handles{cluster_id(p)});
    pcolor=sum(cluster_id(1:p)==cluster_id(p));
    for k=1:14
        subplot(3,5,plot_id(k))
        plot(data{p}(:,1),data{p}(:,k+1)*cells(cluster_id(p),k),'color',cmap{cluster_id(p)}(pcolor,:));
    end
    subplot(3,5,10)
    plot(data{p}(:,1),data{p}(:,2:10)*cellscaling(cluster_id(p),:)','color',cmap{cluster_id(p)}(pcolor,:));
end

for c=1:clusters
    figure(figure_handles{c});
    for k=1:14
        subplot(3,5,plot_id(k))
        xlim([0 T])
        xlabel('time (days)')
        ylabel(vars{k})
        hold off
    end
    subplot(3,5,10)
    xlim([0 T])
    xlabel('time (days)')
    ylabel(vars{15})
    sgtitle(['Dynamics for different patients in cluster ' num2str(c)])
    hold off
end


%optional visible scale adjustment
figure(figure_handles{1})
for i=[2 3 6 7 11 13 14]
    subplot(3,5,i);ylim([0 1.5]*cells(1,i-(i>10)))
end
subplot(3,5,5);ylim([0 5]*cells(1,5))
subplot(3,5,12);ylim([0.5 1.5]*cells(1,11))
subplot(3,5,15);ylim([0.5 1.5]*cells(1,14))

figure(figure_handles{2})
for i=[2 3 6 7 11 13 14]
    subplot(3,5,i);ylim([0 1.5]*cells(2,i-(i>10)))
end
subplot(3,5,5);ylim([0 8]*cells(2,5))
subplot(3,5,12);ylim([0.5 1.5]*cells(2,11))
subplot(3,5,15);ylim([0.5 1.5]*cells(2,14))

figure(figure_handles{3})
for i=[1:4 7 11 12 14 15]
    subplot(3,5,i);ylim([0.5 1.5]*cells(3,i-(i>10)))
end
subplot(3,5,5);ylim([0 70]*cells(3,5))
subplot(3,5,6);ylim([0.5 1.5]*cells(3,6))
subplot(3,5,13);ylim([0 1.5]*cells(3,12))

figure(figure_handles{4})
for i=[2 3 6]
    subplot(3,5,i);ylim([0 1.5]*cells(4,i-(i>10)))
end
subplot(3,5,5);ylim([0 8]*cells(4,5))
subplot(3,5,7);ylim([0.5 1.2]*cells(4,7))
subplot(3,5,11);ylim([0.5 1.2]*cells(4,10))
subplot(3,5,12);ylim([0.8 1.2]*cells(4,11))
subplot(3,5,13);ylim([0 1.2]*cells(4,12))
subplot(3,5,14);ylim([0.5 1.2]*cells(4,13))
subplot(3,5,15);ylim([0.8 1.2]*cells(4,14))

figure(figure_handles{5})
for i=[2 3 7 11 13 14]
    subplot(3,5,i);ylim([0 1.5]*cells(5,i-(i>10)))
end
subplot(3,5,4);ylim([0.8 4]*cells(5,4))
subplot(3,5,5);ylim([0.8 5]*cells(5,5))
subplot(3,5,6);ylim([0 1.1]*cells(5,6))
subplot(3,5,11);ylim([0 1.5]*cells(5,10))
subplot(3,5,12);ylim([0.8 1.1]*cells(5,11))
subplot(3,5,15);ylim([0.8 1.1]*cells(5,14))