function [Erreur,neOptimal,TeOptimal,MecanismeOptimal,MoyenneOptimal,ErreurOptimal,FitCorrige] = TeCalculErreur_TRG(type,ne,Te,Quench,mecanismes,Temp_I_exp,I_theo,FitCorrige,P,sig_I_theo,Temp_sig_I_exp)

%% ======= Initialisation des erreurs et incertitudes =======
Erreur=zeros(length(ne),length(Te),9);
Moyenne=zeros(length(ne),length(Te),9);
IncertitudeTe=zeros(length(ne),length(Te),9);
IncertitudeNm=zeros(length(ne),length(Te),9);
       index(1)=2; %position des Raies dans E
       index(2)=2; %position des Raies dans E
       index(3)=21; %position des Raies dans E
       index(4)=33; %position des Raies dans E
       index(5)=41; %position des Raies dans E

  

   
    
    %% Scan sur tout l'espace des paramètres
    for r=1:length(Quench)
        for j=1:length(ne)    
            for i=1:length(Te)
                %% Extraction de l'intensité théorique à comparer et mise à 0 des raies concernées
                I_simul=zeros(1,41);
                I_exp=zeros(1,41);
                sig_I_exp=zeros(1,41);
                sig_I_simul=zeros(1,41);
                for m=1:41
                    I_simul(m)=I_theo(j,i,r,m)*FitCorrige(m);
                    I_exp(m) = Temp_I_exp(m)*FitCorrige(m);
                    sig_I_exp(m) = Temp_sig_I_exp(m)*FitCorrige(m);
                    sig_I_simul(m)=sig_I_theo(j,i,r,m)*FitCorrige(m);
                end

                if type ==1 % si c'est le w= 1-sig
                    %rejete les raies avec une incertitude plus grande que la
                    %valeur de l'intensité (donne des poids < 0 ...caca)
                    for m=1:length(I_simul)
                        if 1- sqrt( (sig_I_simul(m)/I_simul(m))^2 + (sig_I_exp(m)/I_exp(m))^2 ) < 0
                           I_simul(m)        = 0;
                           sig_I_simul(m)    = 0;
                           I_exp(m)          = 0;
                           sig_I_exp(m)      = 0;
                           %disp('intensité théorique plus petite que l incertitude')
                        end
                    end
                end
                    I_exp=I_exp(I_exp~=0);
                    sig_I_exp=sig_I_exp(sig_I_exp~=0);
                    I_simul=I_simul(I_simul~=0);
                    sig_I_simul=sig_I_simul(sig_I_simul~=0);


                %% poids statistiques pour la STD
                w=zeros(1,size(I_exp,2));
                if     type == 1
                    w = 1 - sqrt((sig_I_simul./I_simul).^2 + (sig_I_exp./I_exp).^2);  %Poid pour la STD
                elseif type == 2
                    w = 1 ./ sqrt((sig_I_simul./I_simul).^2 + (sig_I_exp./I_exp).^2);  %Poid pour la STD
                elseif type ==3
                    w = ones(1,size(I_exp,2));
                end


                %% Initialisation des ratios pour chaque nouveau couple n1s2/Te
                ratio=zeros(1,size(I_exp,2));

                %% Obtention des ratios pour toutes les raies

                for m=1:size(I_exp,2)
                    ratio(m)=I_exp(m)/I_simul(m);
                end
                %% Calcul de l'erreur standard relative
                Moyenne(j,i,r)=mean(ratio);
                Erreur(j,i,r)=(100*std(ratio,w))/Moyenne(j,i,r);


    %             Output(j,i,:) = Erreur; %prend les Ratios en output pis les save apres sous All_ratio.mat

            end
        end
    end
    %   save('All_ratio.mat','Erreur')
    %   save('ne.mat','ne')
    %   save('Te.mat','Te')

        %% 1) Te à la valeur de Ne input
       PosNeMin1 =  round(size(Erreur,1)/2);
       PosQuenchMin1 =  round(size(Erreur,3)/2);
       PosTeMin1= find(Erreur(PosNeMin1,:,PosQuenchMin1 )==min(Erreur(PosNeMin1,:,PosQuenchMin1 )));

        %% 2) Te et Ne en faisant varier Ne mais pas Quench      
       [PosNeMin2,PosTeMin2]=find(Erreur(:,:,PosQuenchMin1)==min(min(Erreur(:,:,PosQuenchMin1))));
           
        %% 3) Te et Ne en faisant varier Quench, le taux de quenching de 2^n : n= -4 à 4
        for r=1:length(Quench)
           MinQuench1(r) = min(Erreur(PosNeMin1,:,r));      
        end
        [val PosQuenchMin3] = min(MinQuench1);
        [val PosTeMin3]=min(Erreur(PosNeMin1,:,PosQuenchMin3));
 
        %% 4) Te et Ne en faisant varier Ne et Quench
        for r=1:length(Quench)
           MinQuench2(r) = min(min(Erreur(:,:,r)));
           
        end
        [val PosQuenchMin4] = min(MinQuench2);
        [PosNeMin4,PosTeMin4]=find(Erreur(:,:,PosQuenchMin4(1))==min(min(Erreur(:,:,PosQuenchMin4(1)))));

       %% Évaluation de Te Ne (la moyenne des 4!) et l'erreur et l'erreur
       
       Te_opt =   [Te(PosTeMin1);
                   Te(PosTeMin2);
                   Te(PosTeMin3);
                   Te(PosTeMin4)];
          
       err=std(Te_opt)/mean(Te_opt);

       TeOptimal(1) = mean(Te_opt);
       TeOptimal(2) = mean(Te_opt)*(1-err);%min
       TeOptimal(3) = mean(Te_opt)*(1+err);%max
       [c PosMeanTeMin] = min(abs(Te-TeOptimal(1)));%position la plus proche
       TeOptimal(4) = PosMeanTeMin(1);

       
       neOptimal(1)=10^((log10(ne(PosNeMin2))+log10(ne(PosNeMin4)))/2);  %moyenne en log  
       
       if ne(PosNeMin2)<ne(PosNeMin4)
           neOptimal(2)=ne(PosNeMin2);%min
           neOptimal(3)=ne(PosNeMin4);%max
       elseif ne(PosNeMin2)>ne(PosNeMin4)
           neOptimal(2)=ne(PosNeMin4);%min
           neOptimal(3)=ne(PosNeMin2);%max      
       else
           neOptimal(2)=neOptimal(1);%min
           neOptimal(3)=neOptimal(1);%max      
       end
       [c PosMeanNeMin] = min(abs(ne-neOptimal(1)));%position la plus proche
       neOptimal(4)=PosMeanNeMin(1);
       %% autre data
       MecanismeOptimal=0;  %on s'en fou
       %on prend le Ne moyen calculé, le Te moyen calculé et 0 Quench
       ErreurOptimal = Erreur(neOptimal(4),TeOptimal(4),PosQuenchMin1);
       
       MoyenneOptimal=Moyenne(neOptimal(4),TeOptimal(4),PosQuenchMin1);
    end

