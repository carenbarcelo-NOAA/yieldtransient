function [Ages, Surv, Ratio, SurvW, C2, MeanAge, R95, C95, Params, Q75idx] = convolution_model(Species,F)
Params=define_Params(Species);
Ages = 1:Params.A;

% comment this out in the unfished case
isFished = Ages >= Params.Af;

L = Params.Linf*(1 - exp(-Params.k.*Ages));
W = L.^3;

% fished case
Surv = cumprod(ones(length(Ages(Ages>=Params.Af)),1)*exp(-(Params.M+F))); %fished
Surv = Surv.*exp(-Params.M)*Params.Af; %fished case
Surv = [zeros(Params.Af-1,1);Surv]; %fished case
SurvW = Surv(:).*W(:); 

MeanAge = sum(SurvW(:).*Ages(:)./sum(SurvW));

% Katie's solution:
Tmax = 100;
T = 0:(Tmax-1);
Tmat = repmat(T(:),[1,4]);
Is = 0:3;
Is = repmat(Is,[Tmax,1]);
Ks = [1 -3 3 -1];
Ks = repmat(Ks,[Tmax,1]);

Numerator = Ks.*(exp(-(Params.M+Is.*Params.k).*Params.Af).*(1-exp(-(Params.M+Is.*Params.k).*Tmat))./...
            (Params.M+Is.*Params.k)) + ...
            Ks.*(exp(-(Params.M+Is.*Params.k).*(Params.Af+Tmax))./(Params.M+F+Is.*Params.k));
        
Denominator = Ks.*exp(-(Params.M+Is.*Params.k).*Params.Af)./(Params.M+F+Is.*Params.k);

Ratio = sum(Numerator,2)./sum(Denominator,2);

r1 = diff(Ratio);
%R95= min(find(r1<0.05));

Ratio = Ratio + 1;
T0 = 50;

Ratio = [ones(T0,1); Ratio(:)];

SurvW = flipud(SurvW);
SurvW=SurvW./sum(SurvW); %Rescaling to have same area under curve

% Finding 75% quantile of YPR
Q=prctile(SurvW(:),75);
Q=quantile(SurvW(:),0.75);
Q75idx=max(find(SurvW(:)>Q));

SurvW2 = repmat(SurvW(:)',[length(Ratio),1]);
SurvW2(1:T0,:) = 0;
for i = 1:length(SurvW)
   SurvW2(T0+i,1:(Params.A-i))=0; 
end

C2 = sum(repmat(Ratio,[1,Params.A]).*SurvW2,2);

c1 = diff(C2);
%C95=find(c1<0.05);

%%%%%%
Ratio = Ratio((T0+1):end); % just take the post-MPA part
SurvW = fliplr(SurvW(:)');  % orient it properly in time
SurvW0 = [SurvW, zeros(1,length(Ratio)-Params.A)]; % pad with zeros at end so it matches the size of Ratio
SurvW2 = repmat(SurvW0(:)',[length(Ratio),1]);
C2a = repmat(Ratio(:),[1,length(Ratio)]).*SurvW2; % each point on biomass trajectory * influence function

% Now use a loop to time-shift the convolution prior to summing
for t = 1:length(Ratio)
   C2a(t,:) = [zeros(1,t-1),C2a(t,1:(length(SurvW0)-t+1))]; 
end

C2 = sum(C2a);

Species
% find 95% points
 R95=find(Ratio >= 0.95*max(Ratio),1)
 %Convolution 95%
 C95=find(C2 >= 0.95*max(C2),1)

end







