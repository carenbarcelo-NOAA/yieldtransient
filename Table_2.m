% Barceló et al. Table 2
% Biomass, yield timescales

function Table_2

Sp = {'Kelp rockfish','Blue rockfish','Black rockfish','Gopher rockfish',...
  'Lingcod','Copper rockfish','California scorpionfish','Brown rockfish','Vermilion rockfish',...
'Yellowtail rockfish','Cabezon','China rockfish','Kelp greenling',...
'Kelp bass','Olive rockfish','Black and yellow rockfish'}

genus = {'Sebastes mystinus','Sebastes atrovirens','Sebastes melanops','Sebastes carnatus'...
    'Ophiodon elongatus','Sebastes caurinus','Scorpaena guttata','Sebastes auriculatus','Sebastes miniatus',...
    'Sebastes flavidus','Scorpaenichthys marmoratus','Sebastes nebulosus','Hexagrammos decagrammus',...
    'Paralabrax clathratus','Sebastes serranoides','Sebastes chrysomelas'}
    
Fs=[0.17,0.17,0.05,0.17,0.23,0.08,0.19,0.13,0.14,0.06,0.12,0.09,0.17,0.12,0.07,0.17]; % F vector for each species
linS = {'--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':'};

T=1:50;

MeanAge=zeros(length(Sp),1);
B95=zeros(length(Sp),1);
C95=zeros(length(Sp),1);
B50=zeros(length(Sp),1);
C50=zeros(length(Sp),1);

M=zeros(length(Sp),1);
Ac=zeros(length(Sp),1);
Linf=zeros(length(Sp),1);
Amat=zeros(length(Sp),1);
k =zeros(length(Sp),1);
MF_M=zeros(length(Sp),1);
maxY=zeros(length(Sp),1);
Species=zeros(length(Sp),1);
Ratio=zeros(length(Sp),1);

output_table2=zeros(length(Sp),11);


for i=1:length(Sp)
    fs=Fs(i);

    [B, Ages, SurvW, C2a, MeanAge(i),B95(i), C95(i), YPRint(i), B50(i), C50(i), ratio(i), maxY(i), Params]=convoplots(fs,Sp{i},linS{i},'Ac');

    M(i,:)=Params.M;
    Ac(i,:)=Params.Af;
    Amax(i,:)=Params.A;
    Amat(i,:)=Params.Amat;
    Linf(i,:)=Params.Linf;
    k(i,:) =Params.k;
    MF_M(i,:)=(M(i)+fs)/fs;
    maxY(i,:)=maxY(i);
    Ratio(i,:) =ratio(i);
    
end

Species = regexp(Sp, '\s/', 'split');
Genus = regexp(genus, '\s/', 'split');

colnames = {'Species' 'Genus' 'F' 'M' 'Ac' 'Amax' 'Amat' 'Linf' 'k' '(M+F)/M' 'YPR-MeanAge' 'B95' 'Y95' 'Ratio' }; %'B50' 'Y50'

output_table2=table(Species(:), Genus(:), Fs(:), M(:) ,Ac(:) ,Amax(:), Amat(:) ,Linf(:) ,k(:), MF_M(:), MeanAge(:), B95(:), C95(:), Ratio(:), 'VariableNames',colnames); %B50(:), C50(:)
writetable(output_table2,'Table2_revision.csv')

keyboard

end




