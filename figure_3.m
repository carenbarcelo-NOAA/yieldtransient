% Barceló et al. FIGURE 3
% Biomass, convolution, YPR and convolution/AUC-YPR
function figure_3
linS = {'--','-','-.',':'};
Sp={'Kelp rockfish','Kelp greenling','Gopher rockfish','Blue rockfish'};
Fs=[0.17,0.17,0.17,0.17];
mycolours=flipud(lines(4));
Ac_Am_scenario ='Ac';

T=1:50;
figure;

for i=1:length(Sp)
    fs=Fs(i);
    plotHandles = zeros(1,i);
    plotLabels = cell(1,i); 

    [B, Ages, SurvW, C2a, MeanAge,B95, C95, YPRint, B50, C50, ratio, maxY, Params]=convoplots(fs,Sp{i},linS{i},Ac_Am_scenario);

    subplot(4,1,1)
    axis tight manual % this ensures that getframe() returns a consistent size
    hold on
    plot([-10 -0.0001],[B(1),B(1)],'Color','r','LineWidth',2,'LineStyle','-');
    plotHandles(i)=plot(T(1:50)-1,B(1:50),'Color',mycolours(i,:),'LineWidth',2, 'LineStyle',linS{i});
    plotLabels{i}=[Sp(i)];
    xlim([-5 45])
    ylim([0 5])
    ylabel('Recruitment index', 'fontsize', 14)
    xlabel('Time (y)','fontsize',14)
    plot([0 0],[0 0],'k^','MarkerFaceColor','k') %[0 sum(C2a(:,end)).*2]
    plot([B95 B95],[0 5],'Color','k','LineStyle','-') %[0 sum(C2a(:,end)).*2]
    box on
    set(gcf,'color','w')

    subplot(4,1,3)
    axis tight manual % this ensures that getframe() returns a consistent size
    plot([-10 -0.00001],[sum(C2a(1)),sum(C2a(1))],'Color', [0.5,0,0.5],'LineWidth',2,'LineStyle','-');
    hold on
    plotHandles(i)=plot(T-1, sum(C2a(:,1:50))/max(sum(C2a(:,1:50))),'LineWidth',2,'Color', mycolours(i,:),'LineWidth',2,'LineStyle',linS{i});
    plotLabels{i}=Sp(i);
    xlim([-5 45])
    ylim([0 1.1])
    xlabel('Time (y)','fontsize',14)
    ylabel({'Yield index'},'fontsize',14)
    plot([0 0],[0 3.5e5],'k^','MarkerFaceColor','k') %[0 sum(C2a(:,end)).*2]
    plot([C95 C95],[0 3.5e5],'Color','k','LineStyle','-'); %[0 sum(C2a(:,end)).*2]     

    subplot(4,1,4)
    axis tight manual % this ensures that getframe() returns a consistent size
    plot([-10 -0.00001],[sum(C2a(1)),sum(C2a(1))],'Color', [0.5,0,0.5],'LineWidth',2,'LineStyle','-');
    hold on
    plotHandles(i)=plot(T-1, sum(C2a(:,1:50))/max(sum(C2a(:,1:50))),'LineWidth',2,'Color', mycolours(i,:),'LineWidth',2,'LineStyle',linS{i});
    plotLabels{i}=Sp(i);
    xlim([-1 10])
    ylim([0 0.8])
    xlabel('Time (y)','fontsize',14)
    ylabel({'Yield index'},'fontsize',14)
    plot([0 0],[0 3.5e5],'k^','MarkerFaceColor','k') %[0 sum(C2a(:,end)).*2
    
    subplot(4,1,2)
    axis tight manual % this ensures that getframe() returns a consistent size
    hold on
    plotHandles(i)=area(Ages,SurvW,'FaceColor',mycolours(i,:),'FaceAlpha',0.5,'EdgeColor',mycolours(i,:),'LineStyle',linS{i},'LineWidth',2); %
    plotLabels{i}=Sp(i);
    xlim([-2 30])
    ylim([0 .3])
    ylabel('Yield per recruit','fontsize',14)
    xlabel('Ages','fontsize',14)
    box on
    set(gcf,'color','w');
    hold on
end

    % mean age lines
    plot([6.9130 6.9130],[0 3.5e5],'Color','k','LineStyle','-'); %[0 sum(C2a(:,end)).*2]
    plot([6.5627 6.5627],[0 3.5e5],'Color','k','LineStyle','-'); %[0 sum(C2a(:,end)).*2]
    plot([8.8975 8.8975],[0 3.5e5],'Color','k','LineStyle','-'); %[0 sum(C2a(:,end)).*2]
    plot([8.7472 8.7472],[0 3.5e5],'Color','k','LineStyle','-'); %[0 sum(C2a(:,end)).*2]

    subplot(4,1,1)
    ht = text(0,1.2,'MPA begins');
    set(ht,'Rotation',90)
    set(ht,'FontSize',9)
    
    subplot(4,1,3)
    ht = text(0,0.18,'MPA begins');
    set(ht,'Rotation',90)
    set(ht,'FontSize',9)
    
    subplot(4,1,4)
    ht = text(0,0.1,'MPA begins');
    set(ht,'Rotation',90)
    set(ht,'FontSize',9)
        
    subplot(4,1,3)
    ht = text(0,0.18,'MPA begins');
    set(ht,'Rotation',90)
    set(ht,'FontSize',9)
    
    subplot(4,1,4)
    ht = text(0,0.1,'MPA begins');
    set(ht,'Rotation',90)
    set(ht,'FontSize',9)
    
    set(gca,'xcolor','k','ycolor','k')
    set(gca,'tickdir','out','ticklength',[0.015 0.015])
    
end
