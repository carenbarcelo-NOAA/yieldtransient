
function [N0, B0, Y0, C0, N, B, Y, C, YE, Params]=MPA_sims(npatches,MPA_frac,F,T,Conn_scenario,DD_scenario,Species,Lambda_target,Noise)
% Blueprint for MPA simulations
% Define function with input arguments: species (make this a string), how
% many patches, possibly F or FLEP, connectivity scenario, etc.
if ~any(strcmp(Conn_scenario,{'Open','Closed','Porous'}))
            error('Conn_scenario must be Open, Closed or Porous')
end
if ~any(strcmp(DD_scenario,{'None','BH'}))
            error('Conn_scenario must be None or BH')
end
% At present this code only really works for 2 patches. Would have to
% revise to work for more
if npatches ~= 2; error('Only works for 2 patches right now!'); end

% Get parameters
Params = define_Params(Species);
Params.Lambda_target = Lambda_target;

% Assemble patch system
Cmat = zeros(npatches); % square connectivity matrix
Cmat(1:2) = [MPA_frac, 1-MPA_frac];
Cmat(3:4) = fliplr(Cmat(1:2)); 
Params.Cmat = Cmat;

% Start at unfished equilibrium, then simulate initial fishing.
N00 = repmat(Params.SAD(:),[1,2]); % unfished age structure in both patches
C00=repmat(0.0001,[1,2]);
Y00=C00; B00 = N00;
% assemble the Leslie matrices

[Leslies0(:,:,1),Lesliefactor] = get_Leslie(Params,F,Conn_scenario,NaN,Params.Lambda_target); % to-be MPA patch
Leslies0(:,:,2) = Leslies0(:,:,1); %get_Leslie(Params,F,Conn_scenario,Target_Lambda); % fished patch

max(eig(Leslies0(:,:,1)))

[N0 B0 Y0 ~,~,C0] = iterate_model(Params,Leslies0,[F,F],N00,Y00,C00,B00,T,Conn_scenario,DD_scenario, MPA_frac,Noise); % run the model
% The N0 generated here will be the starting point for post-MPA dynamics;

% Now make runs with an MPA
% Redistribute F
newF = F*(1+MPA_frac);

Leslies(:,:,1) = get_Leslie(Params,0,Conn_scenario,Lesliefactor); % to-be MPA patch
Leslies(:,:,2) = get_Leslie(Params,newF,Conn_scenario,Lesliefactor); % fished patch
[N, B, Y, ~,~,C,YE] = iterate_model(Params,Leslies,[0,newF],N0(:,:,end),Y0(:,end),C0(:,end),B0(:,:,end),T,Conn_scenario,DD_scenario, MPA_frac,Noise); % run the model


Ysum=sum(Y);
Y0sum=sum(Y0);
Ytot=horzcat(Y0sum,Ysum);
logYtot=log(Ytot+1);

Csum=sum(C);

