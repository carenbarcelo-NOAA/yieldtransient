function [YE95a, Y50, B95a, B50, Btot, YE, Bmpa, mfm,ratioB]=porous_model(Species, Fs)

T=200;
mycolours=flip(hsv(16));
mycolours = lines(4);
MPA_frac = 0.2;
Lambda_target = 1;

YE95a=zeros(length(Species),1);

for f=1:length(Species)
    
    fs=Fs;

    [N0, B0, Y0, C0, N, B, Y, C, YE, Params]=MPA_sims(2,MPA_frac,fs,T,'Porous','None',Species{f},Lambda_target,[0,0]);

Yadd=YE(2,:);

Bmpa = squeeze(sum(B(Params.Af:end,1,:)));
ratioB=max(Bmpa)/Bmpa(1);

Bsum=squeeze(sum(B((Params.Af:end),2,:)));
B0sum=squeeze(sum(B0((Params.Af:end),2,:)));
Btot=Bsum';
mfm(f)=(Params.M+Fs)/Params.M;

Csum=squeeze(sum(C));
C0sum=squeeze(sum(C0));
Ctot=horzcat(C0sum,Csum);

% Calculate time to asymptote for biomass
BB = Btot;
dBB1=abs((BB-BB(end))./BB(end));
B95 = find(dBB1<=0.05,1,'first');
B95a(f) = find(abs((BB-BB(end))./BB(end))<=0.05,1,'first')-1;
B50(f) = find(abs((BB-BB(end))./BB(end))<=0.50,1,'first')-1;

% Calculate time to asymptote for yield
YY = Yadd(1:end);
dYY1=abs((YY-YY(end))./YY(end));
Y95 = find(dYY1<=0.05,1,'first')-1;

YE95a(f) = find(abs((YY-YY(end))./YY(end))<=0.05,1,'first')-1;
Y50(f) = find(abs((YY-YY(end))./YY(end))<=0.5,1,'first')-1;

end
end
