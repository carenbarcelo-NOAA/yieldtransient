% Barceló et al. FIGURE 2

% Trajectories of Biomass (or Recruitment index) {using equation from Kaplan et al. 2019},...
% YPR, Convolution for yield
function figure_2
T = 0:50; 
Ts = repmat(T(:),[1,4]);

% Biomass Function:
M = 0.14;
F = 0.17;
k = 0.17;
Linf = 38.15;
Ac = 4;

Is = repmat(0:3,[length(T),1]);
Ks = repmat([1 -3 3 -1],[length(T),1]);

%Biomass from Kaplan et al.(2019)
Bnum = Ks.*(exp(-(M+Is.*k).*Ac).*(1-exp(-(M+Is.*k).*Ts)))./(M+Is.*k) ...
    + Ks.*exp(-(M+Is.*k).*(Ac+Ts))./(M+F+Is.*k);

Bdenom = (Ks.*(exp(-(M+Is.*k).*Ac))./(M+F+Is.*k));

cU=0.25;
B = (cU*sum(Bnum,2)./sum(Bdenom,2));

% Influence Function (YPR):
maxA = 21;
Ages = 1:maxA;
isFished = Ages>=Ac;
L = Linf*(1-exp(-k*Ages));
W = L.^3;

Surv = cumprod(ones(length(Ages(Ages>=Ac)),1)*exp(-(M+F)));
Surv = Surv./max(Surv);
Surv = Surv.*exp(-M*Ac);

Surv = [zeros(Ac-1,1);Surv];

SurvW = Surv(:).*W(:);

% Set up the convolution
SurvWnf = SurvW(:)';  % orient it properly in time

SurvW = fliplr(SurvW(:)');  % orient it properly in time
SurvW = SurvW(:)';  % orient it properly in time

SurvW0 = [SurvW, zeros(1,length(T)-maxA)]; % pad with zeros at end so it matches the size of Ratio
SurvW2 = repmat(SurvW0(:)',[length(T),1]);
C2a = repmat(B(:),[1,length(B)]).*SurvW2; % each point on biomass trajectory * influence function

% Now use a loop to time-shift the convolution prior to summing
SurvW2a = SurvW2;
for t = 1:length(B)
   C2a(t,:) = [zeros(1,t-1),C2a(t,1:(length(SurvW0)-t+1))]; 
   SurvW2a(t,:) = [zeros(1,t-1),SurvW2a(t,1:(length(SurvW0)-t+1))]; 
end

B = B/max(B)*3;

% Rescale the functions to same scale
SurvWnf = SurvWnf/max(SurvWnf);
SurvW = SurvW/max(SurvW);

%Plots (later modified in Adobe Illustrator for final versions).
figure;
T=1:50;
subplot(2,2,1)
    hold on
    plot([-50 0],[B(1),B(1)],'r-','LineWidth',3)
    plot(T-1,B(1:50),'r-','LineWidth',3)
    xlim([-10 50])
    ylim([-.2 3.5])
    ylabel('Recruitment index','FontSize',14)
    xlabel('Time (y)','fontsize',14)
    set(gca,'linewidth',1.5);
    set(gca,'YTickLabel',[])
    plot([0 0],[0 3.5],'k--') %[0 sum(C2a(:,end)).*2]
    box on
subplot(2,2,2)
    plot(Ages-1,SurvWnf*0.5-.05,'b-','LineWidth',3)
    xlim([-10 40])
    ylim([-0.05 .5])
    set(gcf,'color','white')
    set(gca,'linewidth',1.5);
    set(gca,'YTickLabel',[])
    box on
    xlabel('Ages (y)','fontsize',14)
    ylabel('YPR','FontSize',14)
subplot(2,2,3)
         plot(Ages-1,SurvWnf*0.5-.05,'b-','LineWidth',3)
    xlim([-10 40])
    ylim([-0.05 .5])
    set(gcf,'color','white')
    set(gca,'linewidth',1.5);
    set(gca,'YTickLabel',[])
    box on
    xlabel('Ages (y)','fontsize',14)
    ylabel('YPR','FontSize',14)  
subplot(2,2,4) 
    T=1:40;
    sC2a=sum(C2a(:,1:t));
    axis tight manual % this ensures that getframe() returns a consistent size
    plot([-10 0],[sum(C2a(1))+50000,sum(C2a(1))+50000],'Color', [0.5,0,0.5],'LineWidth',2)
    hold on
    plot(T-1,sC2a(1:40)+50000,'LineWidth',2,'Color', [0.5,0,0.5],'LineWidth',2)
    xlim([-10 50])
    ylim([40000 sum(C2a(:,end))+50500])%.*2
    xlabel('Time (y)','fontsize',14)
    ylabel('Relative Yield','fontsize',14)
    set(gca,'YTickLabel',[])
    hold on
    yyaxis right
    plot([-10 0],[B(1),B(1)],'r-','LineWidth',3)
    set(gca,'linewidth',1.5);    
    ylim([0 1]);
    ax2 = gca; % current axes
    ax2.YColor = 'k';
    set(gca,'YTickLabel',[])
end
