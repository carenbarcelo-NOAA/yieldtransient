function [N, B, Y, E, EP, C, Yextra] = iterate_model(Params,L,F,N0,Y0,C0,B0,T,Conn_scenario,DD_scenario, MPA_frac,Noise)

% Iterate a Leslie-matrix-based model

%Add in noise to recruitment;
if ~exist('Noise','var')
    Noise = ones(T,1);
else
    if length(Noise) < T 
         Noise = exp(normrnd(Noise(1),Noise(2),T,1)); 
    else
         Noise=exp(Noise); 
    end
end


% Initialize
N = zeros(size(N0,1),size(N0,2),T);
B = N;
E=zeros(size(N0,2),T);
EP = zeros(size(N0,1),T);
D_BM = N(1:end-1,:,:); % has one fewer dimension
Y = zeros(size(N0,2),T);
Y(:,1) = Y0;
C = zeros(size(N0,2),T);
C(:,1) = C0;
N(:,:,1)= N0; % initialize density
B(:,:,1) = B0;

% Add in extra state variable for 'extra' recruits in fished patch
Nextra = zeros(size(N));
Bextra = zeros(size(B));
Yextra = zeros(size(Y));

switch Conn_scenario
    case {'Open'}
            for t = 2:T
                
                R = Params.R * Noise(t); % noise in recruitment
                 
                for i = 1:size(N0,2) % number of patches
                N(:,i,t) = L(:,:,i)*N(:,i,t-1);                
                N(1,i,t)=R; %N(1,i,t)=R*MPA_frac; 
                B(:,i,t)= N(:,i,t).* Params.BiomassAge(:); 
                
                if F(i) > 0
                Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i))); 
                C(i,t)=Y(i,t)./F(i); % CPUE
                end                    
                end % end loop over patches
                
                if Noise(t) < 1

                end
                
            end% end loop over time    
    
    case {'Porous'}
        
                        

        for t = 2:T % time
            
            R = Params.R * Noise(t); % noise in recruitment
            
                for i = 1:size(N0,2) % number of patches
                    N(:,i,t) = L(:,:,i)*N(:,i,t-1);
                    Nextra(:,i,t) = L(:,:,i)*Nextra(:,i,t-1);
                    if i==1 % MPA
                    E(i,t)=N(1,i,t); %Eggs produced by mpa
                                        
                    N(1,i,t)=(R*MPA_frac)./N(1,1,1); 
                    N(1,i,t)=(R*MPA_frac); 

                    B(:,i,t)= N(:,i,t).* Params.BiomassAge(:); 
                    end
                
                    if i==2 % Fished patch
                        if F(1) == 0 
                            Etmp = E(1,t)-E(1,1); 
                    E_obase=R*(1-MPA_frac)+Etmp*(1-MPA_frac); 
                        else
                            Etmp = 0;
                    E_obase=R*(1-MPA_frac);
                        end
                    N(1,i,t) = E_obase;
                    
                    Nextra(1,i,t) = Etmp*(1-MPA_frac);
                
                    end
                    B(:,i,t)= N(:,i,t).* Params.BiomassAge(:); 
                    Bextra(:,i,t) = Nextra(:,i,t).* Params.BiomassAge(:); 
                    if F(i) > 0
                    Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i))); 
                    C(i,t)=Y(i,t)./F(i); % CPUE
                    
                    Yextra(i,t) = sum((Bextra(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i))); 
                    
                    end
                end
        end

    case {'Closed'}
                    
        if size(N0,2) == 1
            dispmat = 1;
        elseif size(N0,2) == 2
        dispmat=[MPA_frac,MPA_frac;(1-MPA_frac),(1-MPA_frac)];
        else
            error('Can only have 1 or 2 patches')
        end

        for t = 2:T %this was 2:T
            for i = 1:size(N0,2) % number of patches
                N(:,i,t)=L(:,:,i)*N(:,i,t-1);%closed pop has no external recruits
                E(i,t)=N(1,i,t)*Noise(t);% %Eggs produced in that patch
                B(:,i,t)= N(:,i,t).* Params.BiomassAge(:);
                EP(:,t) = Params.EP0(:)*Noise(t);
                if F(i) > 0
                Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i))); %
                C(i,t)=Y(i,t)./F(i); % CPUE
                end
             end 
           
            % Larval dispersal & density-dependence 
            Rdis=dispmat*E(:,t);
                 
     switch DD_scenario
         case {'BH'}
            Rdis1=(1/(0.25*Params.LEP0)).*Rdis./(1+((1/(0.25*Params.LEP0))/1).*Rdis);
         otherwise
             Rdis1 = Rdis;
     end % end switch
     
            N(1,:,t)=Rdis1;
            
            for i = 1:size(N0,2)
            B(:,i,t)= N(:,i,t) .* Params.BiomassAge(:); %(N(1:end,t-1)-N(2:end)).*(Ba(1:end-1)+Ba(2:end))/2;
            end 
            
        end 
end

         
