function [rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral] = TeTauxDeReaction_TRG(gaz,Te,ChoixHautePression,ChoixAutoabs,exposant)

%fonction qui sort les Taux de réactions pour Le pas de Te sélectionné,
%pour le gaz utilisé, 
%pour High ou Low Ar pressure
% et qui applique les Fudges Factors

%% ============== Préchargement de variables ==============
   %Taux de réactions pour Te=0.05eV jusquà 10eV par step de 0.001eV
   % #gaz 3=Ar 2=Ne, 4=Kr, 5=Xe

    global K_1s_2p K_gs_1s K_quench_1s K_neutral K_gs_2p

%% ============== Boucle d'initialisation sur Te pour la sélection des bons taux de réaction ==============
        %Préallocation d'espace
        rateGround_1s=zeros(5,length(Te));
        rateGround_2p=zeros(10,length(Te));          
        rate1s_2p=zeros(5,10,length(Te));
        rateQuenching=zeros(5,7,length(Te));
        rateNeutral = zeros(6,5);    
        sig_rateGround_1s=zeros(5,length(Te));
        sig_rateGround_2p=zeros(10,length(Te));          
        sig_rate1s_2p=zeros(5,10,length(Te));
        sig_rateQuenching=zeros(5,7,length(Te));
        sig_rateNeutral = zeros(6,5);



        
        %%
         %si on considère les résonants ou non
          if ChoixAutoabs==1
              pas=1; % on prend tout
          else
              pas=2; % on prend que le 1s3 et 1s5 (sans autoabs/résonant)
          end
        
    for j=1:length(Te)       %Boucle sur Te
        %D'abord on trouve la position du Te correspondant dans Allrates
        position=round((Te(j)/0.01)-4);    %-4pour retomber à une position 1 pour Te=0.05 et matcher le pas de TeTauxdeReaction pas de 0.01 de 0.05 à 50 eV     
        e = round((exposant/0.01)-9);           %   position de l'exposant                

      %  position=round((Te(j)/0.1)-4);        %pour ceux de TEOES , enlever fudge                         
        %Puis on sélectionne les taux du fondamental aux 1s via impact électronique

        for l=1:pas:5
            rateGround_1s(l,j)=K_gs_1s(e,gaz,l,position);
        end 

            

        %...les taux du fondamental aux 2p via impact électronique(10 niveaux=1:10)
        for l=1:10
            rateGround_2p(l,j)=K_gs_2p(e,gaz,l,position);
        end 

        %...les taux des 1s aux 2p via impact électronique
        for k=1:pas:5         %Boucle sur les 1s
            for l=1:10    %Boucle sur les 2p
                rate1s_2p(k,l,j)=K_1s_2p(e,gaz,k,l,position);
            end
        end
        
         %...et les taux de quenching des 1s par impact électronique (mixing, collision superélastique et ionisation)
        for k=1:pas:5                       %Boucle sur les 1s
            for l=1:7                  %Boucle sur les 1s; 6=coll superélastique, 7=ionisation
                rateQuenching(k,l,j)=K_quench_1s(e,gaz,k,l,position);
            end
        end
    end
    %Constant pour tout Te, quenching des 1s par les autres gaz présents (mat 5X5)
    rateNeutral=K_neutral(:,gaz);
    
%% Correction par les Fudge Factor de Donnelly Teoes 2007 (Les fudges sont seulement pour les transitions ground -> 2P et 1s -> 2P)
[FUDGEground_2P,FUDGE1s_2P] = TeFudgefactor_TRG(gaz);
%Boucle sur les 1s   pour les fudges
for j=1:5 
    for i=1:10
    rate1s_2p(j,i,:) =  rate1s_2p(j,i,:).*FUDGE1s_2P(j,i); %rate1s_2p(1s,2p,Te) FUDGE1s_2P(1s,2p)
    end
end
for i=1:10
    rateGround_2p(i,:) = rateGround_2p(i,:).*FUDGEground_2P(i); %rateGround_2p(2p,Te) FUDGEground_2P(2p)
end
    clear j i
    
    %% Ajout des incertitudes sur les Taux de réactions.
    %selon Donnelly Teoes 2007
    % Toutes les Gas  K_neutral            = +/- 25%
    % Ne, Ar, Kr, Xe  Gs -> 2P                 = +/- 15% x
    % Ne, Ar, Kr      1s -> 2P                 = +/- 25% x
    % Xe              1s -> 2P                 = +/- 45% x
    % Ne, Ar, Kr      Gs -> 1s                 = +/- 15% x
    % Xe              Gs -> 1s                 = +/- 30% x 
    % Ne, Ar, Kr, Xe  ionization               = +/- 25%
    % Ne, Ar, Kr, Xe  quenching by electrons   = +/- 40%
    % ou le mixing est dans le quenching by electron.
    
    if gaz==5 % différent pour le Xenon
        sig_rateGround_1s = 0.30*rateGround_1s;
        sig_rate1s_2p     = 0.45*rate1s_2p;
    else
        sig_rateGround_1s = 0.15*rateGround_1s; 
        sig_rate1s_2p     = 0.25*rate1s_2p;
    end 
    
    sig_rateGround_2p = 0.15*rateGround_2p;
    % ou les fudges ne sont pas appliqués ici puisqu'ils le sont sur rate1s_2p
    % rateGround_1s avant le calcul des sig
    
    sig_rateQuenching(1:5,1:6,:) = 0.40*rateQuenching(1:5,1:6,:); %mixing
    sig_rateQuenching(1:5,7,:) = 0.25*rateQuenching(1:5,7,:); %ionisation
    
    sig_rateNeutral = 0.25*rateNeutral;

    
    
 clear K_1s_2p K_gs_2p K_gs_1s K_quench_1s K_neutral position
   