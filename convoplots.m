% Biomass, Convolution and YPR figure 
function [B, Ages, SurvW, C2a, MeanAge,B95, C95, YPRint, B50, C50,ratio, maxY, Params]=convoplots(F,Species,linS)
linS=linS;
Params = define_Params(Species);

% Time vector for biomass:
T = 0:1000;
Ts = repmat(T(:),[1,4]); % expand for the vonBert polynomial expansion below

% Biomass Function:
M = Params.M; % natural mortality rate
F = F; % fishing mortality rate
k = Params.k; % kappa
Linf = Params.Linf; % Linf
Ac = Params.Af; % Age at first capture
Am =Params.Amat; % Age at maturity

Is = repmat(0:3,[length(T),1]);
Ks = repmat([1 -3 3 -1],[length(T),1]);

Ax=Ac;
      
% Equation 2 - recruitment index
    Bnum_x = Ks.*(exp(-(M+Is.*k).*Ax).*(1-exp(-(M+Is.*k).*Ts)))./(M+Is.*k) ...
    + Ks.*exp(-(M+Is.*k).*(Ax+Ts))./(M+F+Is.*k);
    
    Bdenom_x = (Ks.*(exp(-(M+Is.*k).*Ax))./(M+F+Is.*k));
   
    B = (sum(Bnum_x,2)./sum(Bdenom_x,2));
   
% Influence Function: YPR
maxA = Params.A;
Ages = 1:maxA;
isFished = Ages>=Ac;
L = Params.Linf*(1-exp(-Params.k*Ages));
W = L.^(3);

Surv = cumprod(ones(length(Ages(Ages>=Ac)),1)*exp(-(Params.M+F)));
Surv = Surv./max(Surv);
Surv = Surv.*exp(-Params.M*Ac);

Surv = [zeros(Ac-1,1);Surv];
SurvW = Surv(:).*W(:);

% Set up the convolution
SurvW = SurvW(:)';  % orient it properly in time
SurvW0 = [SurvW, zeros(1,length(T)-maxA)]; % pad with zeros at end so it matches the size of Ratio
SurvW2 = repmat(SurvW0(:)',[length(T),1]);
C2a = repmat(B(:),[1,length(SurvW2)]).*SurvW2; % each point on biomass trajectory * influence function

% Now use a loop to time-shift the convolution prior to summing
SurvW2a = SurvW2;
for t = 1:length(B)
   C2a(t,:) = [zeros(1,t-1),C2a(t,1:(length(SurvW0)-t+1))]; 
   SurvW2a(t,:) = [zeros(1,t-1),SurvW2a(t,1:(length(SurvW0)-t+1))]; 
end

% rescale the functions to same scale
SurvW2a = SurvW2a/max(SurvW2a(:))*.2;
SurvW=SurvW./sum(SurvW); %Rescaling to have same area under curve

%Find 95% times (use find() as index to T, which runs from 0:maxT)
%B95=find(abs((B-B(end))./B(end))<=0.05,1,'first')-1;
B95=T(find(B >= 0.95*max(B),1));
B50=T(find(B >= 0.50*max(B),1));

C95=T(find(sum(C2a) >= 0.95*(max(sum(C2a))),1));
C50=T(find(sum(C2a) >= 0.50*(max(sum(C2a))),1));

MeanAge = sum(SurvW(:).*Ages(:)./sum(SurvW));
YPRint=trapz(Ages(:),SurvW(:));% Area under the curve for YPR function
ratio = max(B)/B(1);
maxY =max(sum(C2a(:,1:50)));

end
