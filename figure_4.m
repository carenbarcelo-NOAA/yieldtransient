% Barceló et al. ICES JMS FIGURE 4

% YPR plots for all 16 species

function [Species ma b95 c95] = figure_4

%Multipanel figure to show different YPR shapes

%List of species to plot
Species = {'Blue rockfish','Kelp rockfish','Black rockfish','Gopher rockfish',...
'Lingcod','Copper rockfish','California scorpionfish','Brown rockfish','Vermilion rockfish',...
'Yellowtail rockfish','Cabezon','China rockfish','Kelp greenling',...
'Kelp bass','Olive rockfish','Black and yellow rockfish'}
Colors = hsv(length(Species));

F=[0.17,0.17,0.05,0.17,0.23,0.08,0.19,0.13,0.14,0.06,0.12,0.09,0.17,0.12,0.07,0.17]; % F vector for each species
T = 50;
MPA_Frac = 0.2;

c95=zeros(16);
b95=zeros(16);

figure;
for i = 1:length(Species)
    
[Ages, Surv, Ratio, SurvW, C2, MeanAge, R95, C95, Params, Q75idx] = convolution_model(Species{i},F(i));
%[M,I] = max(change);

ma(i)=MeanAge;
c95(i)=C2;
b95(i)=R95;

subplot(4,4,i);
area(Ages,SurvW,'FaceColor',Colors(i,:),'FaceAlpha',0.5);
hold on
plot([MeanAge MeanAge], [0 0.3],'k--','color',Colors(i,:))
    ylabel('Surv*W');
    xlabel('Age (y)');
    ylim([0 0.3])
    xlim([0 50])
    set(gca,'FontSize', 15)
    set(gcf,'Color','w')
    grid on
    xticks([0 10 20 30 40 50])
    xticklabels({'0','20','30', '40', '50'})
    title(Species{i}) 
    box on
end
end

