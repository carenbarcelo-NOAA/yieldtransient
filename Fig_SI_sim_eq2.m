%Barcelo et al. ICES JMS Supplementary materials:
%Comparing Eq 2 with simulations: 
%Plots of trajectories for biomass and yield

function [YE95a C95 B95a B95 Sp] = Fig_SI_sim_eq2
linS = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
Sp = {'Kelp rockfish','Blue rockfish','Black rockfish','Gopher rockfish',...
  'Lingcod','Copper rockfish','California scorpionfish','Brown rockfish','Vermilion rockfish',...
'Yellowtail rockfish','Cabezon','China rockfish','Kelp greenling',...
'Kelp bass','Olive rockfish','Black and yellow rockfish'}
Fs=[0.17,0.17,0.05,0.17,0.23,0.08,0.19,0.13,0.14,0.06,0.12,0.09,0.17,0.12,0.07,0.17];%,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.17];
mycolours=flipud(lines(16));

T=1:100;
figure;

for i=1:length(Sp)
    fs=Fs(i);
    plotHandles = zeros(1,i);
    plotLabels = cell(1,i); 
    
    %% Trajectory of Eq 2 (ORIGINAL version with Ac)
    [B(i,:), Ages, SurvW, C2a, MeanAge, B95(i,:), C95(i,:), YPRint, B50, C50,ratio, maxY, Params]=convoplots(fs,Sp{i},linS{i});
    
    sC2a =sum(C2a);
    
    subplot(4,4,i)
    plot([-10 -0.0001],[sC2a(1),sC2a(1)],'Color','r','LineWidth',1,'LineStyle','-');
    plotHandles(i)=plot(T(1:100),sC2a(1:100)/sC2a(end),'Color',mycolours(i,:),'LineWidth',1, 'LineStyle','-');
    hold on
    plot([C95(i) C95(i)],[0 5],'Color',mycolours(i,:),'LineStyle','-') %[0 sum(C2a(:,end)).*2]
    xlim([0 70])
    ylim([0 2])
    ylabel('Added Yield', 'fontsize', 14)
    xlabel('Time (y)','fontsize',14)
    hold on
    box on
    set(gcf,'color','w')
    title(Sp(i));

    %% Trajectory from simulations
    [YE95a(i,:), Y50(i,:), B95a(i,:), B50(i,:), Btot(i,:), YE, Bmpa(i,:), mfm(i,:),ratioB(i,:)]=porous_model(Sp(i), Fs(i))

    Yaddnorm(i,:)=YE(2,:);
    plotHandles(i)=plot(Yaddnorm(i,:)./Yaddnorm(i,end),'Color','k','LineWidth',1, 'LineStyle','-');
    plot([YE95a(i) YE95a(i)],[0 5],'Color','k','LineWidth',1, 'LineStyle','-');

   

end