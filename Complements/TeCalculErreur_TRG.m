function [Erreur,neOptimal,TeOptimal,MoyenneOptimal,ErreurOptimal,FitCorrige] = TeCalculErreur_TRG(type,ne,Te,Temp_I_exp,I_theo,FitCorrige,P,sig_I_theo,Temp_sig_I_exp)

%% fonction qui calcul l'erreur standard relative entre les raies experimentales et théotique
%% Donne Te et Ne à l'optimal

%% ======= Initialisation des erreurs et incertitudes =======
Erreur=zeros(length(ne),length(Te));
Moyenne=zeros(length(ne),length(Te));
IncertitudeTe=zeros(length(ne),length(Te));
IncertitudeNm=zeros(length(ne),length(Te));
       index(1)=2; %position des Raies dans E
       index(2)=2; %position des Raies dans E
       index(3)=21; %position des Raies dans E
       index(4)=33; %position des Raies dans E
       index(5)=41; %position des Raies dans E

 
%% ======= Calcul d'erreur =======
    
    %% Scan sur tout l'espace des paramètres
    for j=1:length(ne)    
        for i=1:length(Te)
            %% Extraction de l'intensité théorique à comparer et mise à 0 des raies concernées
            %initialisation
            I_simul=zeros(1,41);
            I_exp=zeros(1,41);
            sig_I_exp=zeros(1,41);
            sig_I_simul=zeros(1,41);
            %enleve les raies qui auraient été mises à zero par le fit.
            for m=1:41
                I_simul(m)=I_theo(j,i,m)*FitCorrige(m);
                I_exp(m) = Temp_I_exp(m)*FitCorrige(m);
                sig_I_exp(m) = Temp_sig_I_exp(m)*FitCorrige(m);
                sig_I_simul(m)=sig_I_theo(j,i,m)*FitCorrige(m);
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
                %enleve les zeros, on veut des vecteurs de rapport non nuls.
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
           
            %% Obtention des ratios pour toutes les raies
            
            % Initialisation des ratios pour chaque nouveau couple n1s2/Te
            ratio=zeros(1,size(I_exp,2));
            %calcul des ratios
            for m=1:size(I_exp,2)
                ratio(m)=I_exp(m)/I_simul(m);
            end
            %% Calcul de l'erreur standard relative
            Moyenne(j,i)=mean(ratio);
            Erreur(j,i)=(100*std(ratio,w))/Moyenne(j,i);          
        end
    end
    % on save les valeurs de la STD et des vecteurs Ne Te
%    save('All_ratio.mat','Erreur')
%    save('ne.mat','ne')
%    save('Te.mat','Te')
    %% Obtention de la position où se trouve la plus petite erreur

    [PosNeMin, PosTeMin]=find(Erreur==min(min(Erreur)));
    neOptimal(4)=PosNeMin(1);
    TeOptimal(4)=PosTeMin(1);
    MoyenneOptimal=Moyenne(PosNeMin(1),PosTeMin(1));
    ErreurOptimal = Erreur(PosNeMin(1),PosTeMin(1));
    
    %% Obtention de la température et de la densité théorique qui optimisent le spectre expérimental
    if length(PosNeMin) == 1               %S'il n'y a qu'un seul choix à l'optimum
        neOptimal(1)=ne(PosNeMin);
        TeOptimal(1)=Te(PosTeMin);
        disp ('Une seule paire n1s2/Te trouvée') 
    else                                    %S'il y a plusieurs choix à l'optimum
        neOptimal(1)=ne(PosNeMin(1));
        TeOptimal(1)=Te(PosTeMin(1));
    end

    %% Estimation de l'incertitude sur Te de deux méthode
    %premièrement avec un mapping 3D en variant Ne et Te
    for j=1:length(ne)    
        for i=1:length(Te)
            %Si pour un couple Nm,Te en particulier l'erreur est moins de 5% l'erreur minimale...
            if Erreur(j,i)<1.05*min(min(Erreur))
                IncertitudeTe1(j,i)=Te(i);       %Alors on note cette valeur de Te
                IncertitudeNm1(j,i)=ne(j);     %et de Nm car c'est de l'inceritude

            %Si par contre l'erreur pour ce couple est au delà de 5%
            else
                IncertitudeTe1(j,i)=TeOptimal(1);   %Alors on se fou de cette valeur de Te
                IncertitudeNm1(j,i)=neOptimal(1);   %et de Nm et donc on leur attribue l'optimal pour ne pas qu'ils apparaissent dans le calcul
            end
        end
    end
    %Deuxièmement en Variant seulement Te sur le Ne optimal.
        for i=1:length(Te)
            %si sur la courbe de STD en fct de Te lerreur est plus petite que que 5% du min
            if Erreur(PosNeMin,i)<1.05*min(min(Erreur))
                IncertitudeTe2(i)=Te(i); %on note la valeur
            %Si par contre l'erreur pour ce couple est au delà de 5%
            else
                IncertitudeTe2(i)=TeOptimal(1);  %on s'en fou
            end
        end
%% export Te Ne ainsi que les incertitude dans des Vecteurs.
%note la plus petite valeur calculé entre les deux méthode
    if min(min(IncertitudeTe1)) > min(IncertitudeTe2)
        TeOptimal(2)=min(min(IncertitudeTe1));
    else
        TeOptimal(2)=min(IncertitudeTe2);
    end
%note la plus grande valeur calculé entre les deux méthode    
    if max(max(IncertitudeTe1)) < max(IncertitudeTe2)
        TeOptimal(3)=max(max(IncertitudeTe1));
    else
        TeOptimal(3)=max(IncertitudeTe2);
    end

    neOptimal(2)=min(min(IncertitudeNm1)); %on les garde quand meme (ils ne sont plus en output)
    neOptimal(3)=max(max(IncertitudeNm1)); %on les garde quand meme

end

