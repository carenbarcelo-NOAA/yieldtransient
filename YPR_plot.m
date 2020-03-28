function YPR_plot
%Multipanel figure to show different YPR shapes
%Loop over multiple species to do this
%List of species to plot
Species = {'Blue rockfish','Kelp rockfish','Black rockfish','Gopher rockfish',...
  'Lingcod','Copper rockfish','California scorpionfish','Brown rockfish'}
% Species ={'Vermilion rockfish',...
%   'Yellowtail rockfish','Cabezon','China rockfish','Kelp greenling',...
%   'Kelp bass','Olive rockfish','Black and yellow rockfish'}
Colors = hsv(length(Species));
%F=0.2
F=[0.17,0.17,0.05,0.17,0.23,0.08,0.19,0.13];%,
F=[0.14,0.06,0.12,0.09,0.17,0.12,0.07,0.17]; % F vector for each species
T = 50;
MPA_Frac = 0.2;

figure;
for i = 1:length(Species)
    
[Ages, Surv, Ratio, SurvW, C2, MeanAge, R95, C95, Params, Q75idx] = convolution_model(Species{i},0.17);
%[M,I] = max(change);

subplot(4,2,i);
%plot(Ages,SurvW,'-','color',Colors(i,:));
area(Ages,SurvW,'FaceColor',Colors(i,:),'FaceAlpha',0.5);
hold on
plot([MeanAge MeanAge], [0 0.3],'k--','color',Colors(i,:))
%plot([Q75idx Q75idx], [0 0.3], 'k--')
%annotation('textbox',[.9 .5 .1 .2],'String','Text outside the axes','EdgeColor','none')
    ylabel('Surv*W');
    xlabel('Age (y)');
    ylim([0 0.3])
    xlim([0 50])
    set(gca,'xgrid','on')
    set(gcf,'Color','w')
    title(Species{i}) 
end
end

