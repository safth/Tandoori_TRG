function Calcul_Te_TRG_ViaInterface(Overwrite,InfoGraphes,InfoTe,InfoNe,filetype,Tg,P_tot,longueur,FctRep,ChoixErreur,ChoixHautePression,ChoixAutoabs,Commentaire,flow,exposant,Dimension,sig_longueur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CECI EST LE CODE PRINCIPAL POUR CALCULER TE À PARTIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  D'UN SPECTRE CONTENANT UNE TRACE DE GAZ RARES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
warning('off') % ne pas affichier les avertissements lors du calcul des 1s (1s1 qui nexiste pas fait un NAN)
%% ---------- Notes d'utilisation ----------
%====== Les sous fonctions utilisées sont: ======
% Load les raies à fitter ?nearkrxe.txt?
% Fit les raie des spectres avec la fonction IntensiteGaussMax_TRG.m
% Calcul les flow effectifs avec TePumpingSpeed_TRG.m
% Loop Gaz
%   Selectionne les taux de réaction appropriés (précalculé avec CalculTauxDeReactions_TRG).
%   Loop Ne
%       Loop Te
%           Calcul des intensités théoriques pour tout Ne Te (TeIntTheo_TRG.m)
%       End
%   End 
% End
% 
% TeCalculErreur_TRG calcul une %STD pondéré par des poids (100-%erreur) trouvé par un calcul de propagation d?erreur et minimise I_theo/I_exp pour trouver Te et Ne


% Tous ces fichiers ainsi que les sections efficaces sont dans le sous-dossier "Complements"
%------------------------------------------
%%  ###########################################################################################
%%  ######################################### INPUTS ##########################################
%%  ###########################################################################################

% Range des températures électroniques scannées (eV) (pas plus bas que 0.07 eV)
    Te =InfoTe(1):InfoTe(3):InfoTe(2);   
    
% Range de densité électronique scannée (m^-3)
   %5mtorr
   p=1; %changer p quand je change la pression
                                %50 20   10   5     2
debut                       = [5.2 5.2   5.2  5.1   5.2]; %20 10 5 2mtorr
fin                         = [73 72.7  82  78.8   78]; %fin colonne cm %78.8
nb_spectres                 = [70 72    34   31    35]; %35
                                 %50 20      10       5        2
Val_Pente                   = [1e15 6.3e14  4.969e14 3.358e14 1.98e14  ]; 
n_c                         = 2.12e16;

   Z=linspace(debut(p),fin(p),nb_spectres(p)); 
   Ne=  -Val_Pente(p)*Z+Val_Pente(p)*fin(p)+n_c; %4.7209e+15
   
    
% Choix d'affichage des différents graphes (=0 n'affice pas, =1 affiche)
    GraphExp=InfoGraphes(1);         %Graphe de l'extraction des intensité expérimentales des raies ET des fits gaussiens sur les raies
    TeGraph1=InfoGraphes(2);         %Graphe de Iobs/Iexp ou %fond en fonction de Te ou Ne
    TeGraph2=InfoGraphes(3);         %graphe des mécanismes à l'optimum
    TeGraph3=InfoGraphes(5);         %Graphes du fit des raies, 3D de l'erreur et STD
    TeGrapheFinaux=InfoGraphes(4);   %Graphes de l'évolution de NeOptimal et TeOptimal en fonction des fichiers

%% ###########################################################################################
%% ###################################### DÉBUT DU CODE ######################################
%% ###########################################################################################
   %Valeurs de beta et de l'intégrale de convolution trouvé à l'aide de TeBeta.m
   load AllIntegral_new.mat

%% ============== Raies d'émission utilisées dans le calcul ==============
     if flow(3) > 10*flow(4) % si on a 10X plus dargon que de Krypton
         nearkrxe = 'nearkrxeAr.txt'; %met des doublet pis des triplet à des endroits différents!!
     else
         nearkrxe = 'nearkrxe.txt';
     end
    tempfile=fopen(nearkrxe); %Longueurs d'onde et "Doublet"
    data = textscan(tempfile,'%s %s','delimiter',' ','headerlines',6);
    E=str2double(data{1,1});  % Variable avec les 41 raies à fiter en nm
    Doublet=str2double(data{1,2}); %variable qui dit si il faut fiter une gaussienne,2 ou 3 gaussienne
    fclose(tempfile);
    clear data
    
     
    clear Fit750 Fit696 Fit727 Fit826 Fit706 Fit738 Fit840 Fit794 Fit852 Fit751 Fit763 Fit800 Fit810 Fit866 Fit801 Fit842 Fit811
 
%% ###########################################################################################
%% ####################### DÉBUT DE L'ANALYSE EXPÉRIMENTALE DES SPECTRES #####################
%% ###########################################################################################

%% ========== Identification de tous les fichiers qui correspondent à Te dans le répertoire sélectionné ========== 
if filetype == 0                % Les spectres proviennent du Isoplane
    FilesTe=dir('*.csv');                     
elseif filetype == 0.1          % Les spectres proviennent du Isoplane
    FilesTe=dir('*.csv');                     
elseif filetype == 1            % Les spectres proviennent d'un Avantes
    FilesTe=dir('*.trt');                     
elseif filetype == 2            % Les spectres proviennent de Avantes Evolved
    FilesTe=dir('*.trtx');                     
elseif filetype == 3            % Les spectres proviennent de Reetesh
    FilesTe=dir('*.txt');                     
elseif filetype == 4            % Les spectres proviennent d'Avantes 8
    FilesTe=dir('*.TXT');                     
end
fileNames = {FilesTe.name};  % Extraction des noms de ces fichiers
clear FilesTe                % On n'a plus besoin de FilesTe

%% ========== Préallocation d'espace ==========
NumFichier      =zeros(length(fileNames),1);
mecanismeOptimal=zeros(length(fileNames),1);
NeOpt           =zeros(length(fileNames),1);
NeMin           =zeros(length(fileNames),1);
NeMax           =zeros(length(fileNames),1);
TeOptimal       =zeros(length(fileNames),1);
TeMin           =zeros(length(fileNames),1);
TeMax           =zeros(length(fileNames),1);
MoyenneOptimale =zeros(length(fileNames),1);
ErreurOptimal   =zeros(length(fileNames),1);

I_exp           =zeros(length(fileNames),length(E));
sig_I_exp       =zeros(length(fileNames),length(E));
FitCorrige      =zeros(length(fileNames),length(E));

%% ============= Extraction de l'intensité des raies de chaque spectre (i.e. fichier) =============
disp ('Extraction de l''intensité des raies des fichiers...')
h = waitbar(0,'Fit des spectres') ;
for k=1:length(fileNames)
    Fichier=fileNames{k}          %Affiche à l'écran quel fichier est actuellement analysé
    NumFichier(k)=k;              %Vecteur servant aux graphiques finaux de Te et Nm pour chaque fichier
   
    %% ============= Extraction des données brutes (pour fichier .csv) =============
    if filetype == 0
        tempdata=dlmread(fileNames{k});         % Extraction des data
        lambdaTe=tempdata(:,1);                 % Longueurs d'ondes des mesures
        intTe=tempdata(:,2);                    % Intensité des mesures          
        clear tempdata                          % On a plus besoin de tempdata  
        Int_time_correction=1;
    end   
    %% ============= Extraction des données brutes de Morgane (pour fichier .csv) =============
    if filetype == 0.1 %fuck it
        tempdata=dlmread(fileNames{k});         % Extraction des data
        lambdaTe=tempdata(:,1);                 % Longueurs d'ondes des mesures
        intTe=tempdata(:,3);                    % Intensité des mesures POUR MORGANE 
        clear tempdata                          % On a plus besoin de tempdata
        Int_time_correction=1;
    end    
    %% ============= Extraction des données brutes (pour fichier .trt) =============
    if filetype == 1
        file1=fopen(fileNames{k});
        tempdata=textscan(file1, '%f%f','delimiter',';', 'headerlines',8);  % Extraction des datas 
        lambdaTe=tempdata{1};                   % Longueurs d'ondes des mesures
        intTe=tempdata{2};                      % Intensité des mesures
        clear tempdata                          % On a plus besoin de tempdata   
        Int_time_correction=1;
    end    
    %% ============= Extraction des données brutes (pour fichier .trtx) =============
    if filetype == 2
        [Int_time_correction] = IntTime_trtx(fileNames{k});    %function qui prend les temp d'integration des deux spectros utilisés  
        file1=fopen(fileNames{k});
        tempdata=textscan(file1, '%f%f','delimiter',';', 'headerlines',10);  % Extraction des datas
        lambdaTe=tempdata{1};                   % Longueurs d'ondes des mesures
        intTe=tempdata{2};                      % Intensité des mesures
        clear tempdata                          % On a plus besoin de tempdata 
    end  
    %% ============= Extraction des données brutes (pour fichier .txt) =============
    if filetype == 3
        file1=fopen(fileNames{k});
        tempdata=textscan(file1, '%f%f','delimiter',' ');  % Extraction des datas
        lambdaTe=tempdata{1};                   % Longueurs d'ondes des mesures
        intTe=tempdata{2};                      % Intensité des mesures
        clear tempdata                          % On a plus besoin de tempdata 
        Int_time_correction=1;
    end     
    %% ============= Extraction des données brutes (pour fichier .TXT (Avantes8)) =============
    if filetype == 4
        file1=fopen(fileNames{k});
        tempdata=textscan(file1, '%s\t%s\t%s\t%s','delimiter',';', 'headerlines',9);  % Extraction des datas 
        lambdaTe=tempdata{1}  ;                 % Longueurs d'ondes des mesures
        lambdaTe=str2double(strrep(lambdaTe,',','.'));
        intTe=tempdata{2};                      % Intensité des mesures
        intTe=str2double(strrep(intTe,',','.'));
        clear tempdata                          % On a plus besoin de tempdata 
        Int_time_correction=0.1;
    end
    
    %% ============= Correction par la fonction de réponse si nécessaire =============
    if FctRep~= 0
       intTe=intTe./FctRep(:,2);    % Ireel=Imesure/Fct de reponse
    end

    %% ============= Correction pour un décallage éventuel en longueur d'onde fixé sur la 763nm =============
    %lambdaTe(769:end)=lambdaTe(769:end)+0.5;
%     w=1;
%     for j=1:length(lambdaTe)
%         if abs(lambdaTe(j)-763.51)<2
%             IntTemp(w)=intTe(j);
%             LambdaTemp(w)=lambdaTe(j);
%             w=w+1;
%         end
%     end
%     
%     Shift=763.51-LambdaTemp(IntTemp==max(IntTemp));
%     lambdaTe=lambdaTe+Shift;
%     clear Shift w j IntTemp LambdaTemp
    %% ============= Obtention de l'intensité des différentes raies (au max du pic avec fit Gaussien=============

    [IntPeak,sig_IntPeak,Fit] = IntensiteGaussMax_TRG(E,lambdaTe,intTe,GraphExp,Doublet,Overwrite); %La variable Fit dit quelles raies ont effectivement été fittées
    SaveIntPeak(k,:) =IntPeak;
    %% ============== Correction pour les raies qui ne doivent pas être considérée =============
    FitCorrige(k,:)=Fit.*Overwrite;
    I_exp(k,:)=IntPeak'.*FitCorrige(k,:);
    I_exp(k,1:6)=I_exp(k,1:6)*Int_time_correction; %correction pour le integration time du 500, les raies 1 à 6 sont sur le 500-700
    % meme chose pour l'incertitude
    sig_I_exp(k,:)=sig_IntPeak'.*FitCorrige(k,:);
    sig_I_exp(k,1:6)=sig_I_exp(k,1:6)*Int_time_correction; %correction pour le integration time du 500, les raies 1 à 6 sont sur le 500-700
    clear  lambdaTe intTe Fit IntPeak sig_IntPeak
    %check si l'intensité de la raie est plus grande que son incertitude
    for m=1:size(I_exp,2)
        if I_exp(k,m) < sig_I_exp(k,m)
           I_exp(k,m) = 0;
           sig_I_exp(k,m) = 0;
           FitCorrige(k,m) = 0;
           disp('Intensité fitté plus petit que l incertitude.')
        end
    end
       waitbar(k/length(fileNames),h) 
end %fin de la boucle sur les différents spectres
disp ('Intensités extraites des fichiers')
clear Doublet filetype
% save('SaveIntPeak.mat','SaveIntPeak')


    %% Calcul les pression partielle Par la diffusion différentes du mélange de gaz
    
    %calcul des facteurs de corrections
    [Pump,Pression_Calc] = TePumpingSpeed_TRG(P_tot,flow,Dimension)
   % Pump(2:5)=[1 1 1 1];
   % Pump(13)=1;
   % Pump(25)=0;
    P_tot=(133.3224e-3)*P_tot; %mtorr à pascal
    flow=flow.*Pump; %on les applique aux Flows
    flow_tot=sum(flow); %flow total
    P = (P_tot.*flow)/flow_tot; %Pressions partielles de gaz
    % Calcul de la densité de neutres à l'état fondamental avec p=n*k*Tg [m-3]
    %ng(gaz) car P est un vecteur avec les pressions partielles des gas utilisés
    ng=(P)./(1.38064852*10^(-23)*Tg);


%% ###########################################################################################
%% ########################### DÉBUT DE L'ANALYSE THÉORIQUE DU CODE ##########################
%% ###########################################################################################
disp ('Calcul de l''intensité théorique des raies...')


    %% Préallocation d'espace
    I_theo      =zeros(length(Ne),length(Te),length(E)); %les 41 raies dans le meme fichier
    sig_I_theo  =zeros(length(Ne),length(Te),length(E)); %les 41 raies dans le meme fichier
    Thetaij     =ones(length(Ne),length(Te),length(E)); 
    Doppler     =zeros(length(Ne),length(Te),length(E));
    VanDerWaals =zeros(length(Ne),length(Te),length(E));
    Resonant    =zeros(length(Ne),length(Te),length(E));

    Gains2p    =zeros(5,length(Ne),length(Te),10,5); %gaz #1à5 avec 10 2p chacuns
    Pertes2p   =zeros(5,length(Ne),length(Te),10,3);
    densite2p     =zeros(5,length(Ne),length(Te),10);
    ContributionFond    =zeros(5,length(Ne),length(Te),10);
    densite1s     =zeros(5,length(Ne),length(Te),5);
    sig_densite1s     =zeros(5,length(Ne),length(Te),5);

    Gains1s    =zeros(5,length(Ne),length(Te),5,5); %gaz #1à5 avec 10 2p chacuns
    Pertes1s   =zeros(5,length(Ne),length(Te),5,6);
    
    
    %% sélection des taux de réactions
    global K_1s_2p K_gs_1s K_quench_1s K_neutral K_gs_2p
    load Allrates1s_2p.mat 
    load AllratesGround_1s.mat 
    load AllratesQuenching_1s.mat
    load Allrates_neutral.mat
    if ChoixHautePression==0
      load AllratesGround_2p.mat    %AllratesGround_2p(#gas,2Px,Te)
      K_gs_2p = AllratesGround_2p;
    elseif ChoixHautePression==1 %si on est à Haute pression partielles d'Argon (>1mtorr)
      load AllratesGround_2p_HighP.mat    %AllratesGround_2p(#gas,2Px,Te)
      K_gs_2p = AllratesGround_2p_HighP;
    end
    
    K_1s_2p     =Allrates1s_2p;
    K_gs_1s     =AllratesGround_1s;
    K_quench_1s =AllratesQuenching_1s;
    K_neutral   =rates_neutral;
    
    
    %% Boucles sur les gaz, éléments théoriques du code.
    wait=0;
    global gaz_i gaz_f

for gaz=gaz_i:gaz_f %On fait le calcul théorique de l'intensité de raies pour 2=néon, 3=argon, 4=krypton et 5=xénon
    mecanismes=zeros(length(Ne),length(Te));
    %pour mettre les valeurs dans un vecteur avec tous les gaz inclus
       if gaz==2 %la raie de Neon
       index(1)=1; %position des Raies dans E
       index(2)=2; %position des Raies dans E
       end
       if gaz==3 %les 17 raies dargon
       index(1)=3;  %position des Raies dans E
       index(2)=21; %position des Raies dans E
       end
       if gaz==4 %les 11 raies de Krypton
       index(1)=22; %position des Raies dans E
       index(2)=33; %position des Raies dans E
       end
       if gaz==5 %les 7 raies de Krypton
       index(1)=34; %position des Raies dans E
       index(2)=41; %position des Raies dans E
       end
       waitbar(wait/(length(gaz_i:gaz_f)*length(Ne)),h,'Calcul théorique: sélection des Taux de réactions') 
    
    %% ============== Boucle d'initialisation sur Te pour la sélection des bons taux de réaction ==============
    %fonction qui prend les bon taux de réaction pour un gaz et le pas de
    %Te sélectionné
    [rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral] = TeTauxDeReaction_TRG(gaz,Te,ChoixHautePression,ChoixAutoabs,exposant);
    disp('Taux de réaction selectionnés')

    %%  ============== Première loop sur la densité électronique: COMPTEUR: j ==============
    for j=1:length(Ne)
        waitbar(wait/(length(gaz_i:gaz_f)*length(Ne)),h,'Calcul théorique') 
        wait=wait+1;
        %% Calcul du bilan de population: gains=pertes
       [densite_1s,sig_densite_1s,densite_2p,sig_density,Gains_2p,Pertes_2p,Emission,PopFond,Mecanisms,energie_1s,energie_2p,Gains_1s,Pertes_1s]=TeIntTheo_TRG(gaz,Ne(j),Te,rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral,Tg,longueur,AllIntegral,P,1,0,flow,ChoixAutoabs,sig_longueur,densite1s(3,j,:,:),sig_densite1s(3,j,:,:));
       %% Extraction et mise en mémoire des données pour chaque fichier, densité électronique et Te

       for i=1:length(Te)
           n=0;
            for m=index(1):index(2)                     % Pour chaque raie...
                n=n+1;% pour emission
                I_theo(j,i,m)=Emission(1,i,n);          % Intensité théorique
                sig_I_theo(j,i,m)=Emission(6,i,n);      % Intensité théorique
                Thetaij(j,i,m)=Emission(2,i,n);         % Escape Factor 
                Doppler(j,i,m)=Emission(3,i,n);         % Élargissement Doppler
                VanDerWaals(j,i,m)=Emission(4,i,n);     % Élargissement VanDerWaals
                Resonant(j,i,m)=Emission(5,i,n);        % Élargissement Resonant
            end  
            for l=1:10 
                %for k=1:size(Gains_2p,1)% Pour chaque niveau...         
                Gains2p(gaz,j,i,l,:)=Gains_2p(:,i,l);        %1=fond,2=Nm,3=coll2P,4=radtrap, 5 = Transfert Radiatif 
                Pertes2p(gaz,j,i,l,:)=Pertes_2p(:,i,l);      % 1=rad 2=coll2p 3=coll1s
                %end
                densite2p(gaz,j,i,l)=densite_2p(i,l);            % Densité des 2p
                ContributionFond(gaz,j,i,l)=PopFond(i,l);   % Pourcentage de la population en comparant seulement le fondamental et les métastables    
                
            end
             for l=1:5   
                 Gains1s(gaz,j,i,l,:)=Gains_1s(:,i,l);  
                 Pertes1s(gaz,j,i,l,:)=Pertes_1s(:,i,l); 
                 densite1s(gaz,j,i,l)=densite_1s(i,l);
                 sig_densite1s(gaz,j,i,l)=sig_densite_1s(i,l);
             end
             energie1s(gaz,:)=energie_1s(:);
             energie2p(gaz,:)=energie_2p(:);
           
            mecanismes(gaz,j,i)=Mecanisms(i);

       end   %Fin Boucle Te
    end    %Fin Boucle Ne

    
    %% mise en graphique des %population en fonction de Te ou Ne pour tout les niveaux

   
end %fin boucle sur les gaz

%TeGraphesPercentFond(gaz,Ne,Te,ContributionFond)
a(:,:)=I_theo(1,:,:);
    save('I_theo.mat','a')
%     stop

clear Allrates1s_2p AllratesGround_2p AllratesGround_1s AllratesQuenching_1s rates_neutral position
clear Emission Gains Pertes density PopFond
disp('Intensités théoriques calculées')
%% ###########################################################################################
%%  ####################### DÉBUT DE LA COMPARAISON THÉORIE-EXPÉRIENCE #######################
%%  ##########################################################################################
disp('Comparaison expérience-théorie...')

%% ============== Création des fichiers textes contenant les infos pertinentes pour chaque fichier analysé ==============
if exist('Te.out', 'file')==0 % n'écrit pas le Header si il est deja fait
    fileID = fopen('Te.out','a');               
    formatSpec = '%-22s\t %s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s\r\n';  
    fprintf(fileID,formatSpec,'------NomFichier------','TeOpt','TeMin','TeMax','NeOpt   ','% STD','Poids','exp','HighP','AutoAbs','Commentaire');
else
    formatSpec = '%-22s\t %s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s\r\n';  
    fileID = fopen('Te.out','a');   
end

fileID2 = fopen('Raies Utilisées.out','a');               
formatSpec2 = '%-22s \t%-8s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s\r\n';  
fprintf(fileID2,formatSpec2,'------NomFichier------','-#raies-','585','640','667','696','706','714','727','738','750','751','763','794','800','801','810','811','826','840','842','852','866','758','760','768','769','785','810','811','819','829','850','877','892','788','820','823','828','834','881','895','904');

if exist('Nm_Argon.out', 'file')==0 % n'écrit pas le Header si il est deja fait
    fileID3 = fopen( 'Nm_Argon.out','a');         
    formatSpec3 = '%.4e\t %e\t %.4e\t %e\r\n';
   % fprintf(fileID3,formatSpec3,'1s2','1s3','1s4','1s5');
else
    fileID3 = fopen('Nm_Argon.out','a');   
    formatSpec3 = '%.4e\t %e\t %.4e\t %e\r\n';
end
   

%% ========== Comparaison des spectres expérimentaux un à un avec la théorie ==========
for k=1:length(fileNames)
    Fichier=fileNames{k}          %Affiche à l'écran quel fichier est actuellement comparé 
    
    %% ===== Calcul du nombre total de raies à analyser =====
    nbrderaies(k)=sum(FitCorrige(k,:)); 
    %% Boucles d'évaluation de l'erreur de chaque raie aux températures et densite de 1s regardées
    if nbrderaies(k) > 1             %Vérification s'il y a au moins 2 raies à ajuster      
            %% Calcul de l'erreur
            [STDErr,NeOptimal,TempOptimal,mecanismeOptimal(k),MoyenneOptimale(k),ErreurOptimal(k),FitCorrige(k,:)] = TeCalculErreur_TRG(ChoixErreur,Ne(k),Te,mecanismes,I_exp(k,:),I_theo,FitCorrige(k,:),P,sig_I_theo,sig_I_exp(k,:));

            %% Extraction des résultats
            NeOpt(k)=NeOptimal(1); NeMin(k)=NeOptimal(2); NeMax(k)=NeOptimal(3); PosNeMin(1)=NeOptimal(4);
            TeOptimal(k)=TempOptimal(1); TeMin(k)=TempOptimal(2); TeMax(k)=TempOptimal(3); PosTeMin(1)=TempOptimal(4);
            
            %% 'ecriture des densité de metastable
%              x(:) = densite1s(3,PosNeMin(1),PosTeMin(1),:)
               fprintf(fileID3,formatSpec3,densite1s(3,PosNeMin(1),PosTeMin(1),2:5)); 
            %% Export de toute les fichiers...
%             file(:,:)=Gains2p(3,PosNeMin(1),PosTeMin(1),:,:);
%             file2(:,:)=Pertes2p(3,PosNeMin(1),PosTeMin(1),:,:);
%             file3(:,:)=Gains1s(3,PosNeMin(1),PosTeMin(1),:,:);
%             file4(:,:)=Pertes1s(3,PosNeMin(1),PosTeMin(1),:,:);
%             file5(:,1)=I_theo(PosNeMin(1),PosTeMin(1),:);
%             file5(:,2)=I_exp(1,:);
%             file5(:,3)=1-Thetaij(PosNeMin(1),PosTeMin(1),:);
%             save('Gains2P.mat','file')
%             save('Pertes2P.mat','file2')
%             save('Gains1s.mat','file3')
%             save('Pertes1s.mat','file4')
%             save('I_th-I_exp-AutoAbs.mat','file5')


              
              
             %% Obtention des graphiques à l'optimum si désiré dans les inputs
           if TeGraph2==1
                TeGraphes2_TRG(Ne,Te,STDErr,I_theo(PosNeMin(1),PosTeMin(1),:),sig_I_theo(PosNeMin(1),PosTeMin(1),:),Thetaij(PosNeMin(1),PosTeMin(1),:),I_exp(k,:),sig_I_exp(k,:),FitCorrige(k,:),ChoixErreur,MoyenneOptimale(k),P,E,ChoixAutoabs);
           end
                        %% Obtention des graphiques à l'optimum si désiré dans les inputs
            if TeGraph3==1
                TeGraphes3_TRG(densite2p(:,PosNeMin(1),PosTeMin(1),:),Gains2p(:,PosNeMin(1),PosTeMin(1),:,:),Pertes2p(:,PosNeMin(1),PosTeMin(1),:,:),ContributionFond(:,PosNeMin(1),PosTeMin(1),:),energie1s,energie2p,densite1s(:,PosNeMin(1),PosTeMin(1),:),Gains1s(:,PosNeMin(1),PosTeMin(1),:,:),Pertes1s(:,PosNeMin(1),PosTeMin(1),:,:),ng);
            end
    
            
    else %S'il n'y a pas 2 raies à comparer, on met tout à 0
        TeOptimal(k)=0;
        TeMin(k)=0;
        TeMax(k)=0;
        NeOpt(k)=0;
        NeMax(k)=0;
        NeMin(k)=0;
        mecanismeOptimal(k)=0;
      
        %% Écriture des paramètres optimaux dans un fichier
            formatSpec = '%-22s\t %-5s \t%-5s \t%-5s \t%-12s \t%-5s \t%-5s \t%-5s \t%-13s \t%-5s \t%-20s\r\n';  
            Table1{1}=char(fileNames(k));
            Table1{2}=num2str(TeOptimal(k),4);
            Table1{3}=num2str(TeMin(k),4);
            Table1{4}=num2str(TeMax(k),4);
            Table1{5}=num2str(NeOpt(k),5);
            Table1{6}=num2str(ErreurOptimal(k),4);
            Table1{7}=num2str(ChoixErreur);
            Table1{8}=num2str(exposant);
            Table1{9}=num2str(ChoixHautePression);
            Table1{10}=num2str(ChoixAutoabs);
            Table1{11}=Commentaire;
            fprintf(fileID,formatSpec,Table1{:});
      
      disp('Nombre de raies à comparer insuffisant. Tous les paramètres ont été mis à zéro. You Died.')
    end %fin boucle nbr de raies minimale
    
        %% Mise en graphique des rapports I_obs/_calc pour un Ne donné
     if TeGraph1==1 
         TeGraphes1_TRG(Ne,Te,I_theo(:,:,:),I_exp(1,:),FitCorrige(k,:),MoyenneOptimale(k),NeOptimal(4))
     end
      waitbar(k/length(fileNames),h,'Calcul d''erreur')   
    clear IntensiteSD metSD tempSD Emission Gains Pertes density PosTeMin PosN1s2Min SEGD
end %Fin boucle fichiers
            disp('Écriture des résultats')
for k=1:length(fileNames)
            %% Écriture des paramètres optimaux dans un fichier
            formatSpec = '%-22s\t %-5s \t%-5s \t%-5s \t%-10s \t%-6s \t%-5s \t%-5s \t%-5s \t%-6s \t%-20s\r\n';  
            Table1{1}=char(fileNames(k));
            Table1{2}=num2str(TeOptimal(k),4);
            Table1{3}=num2str(TeMin(k),4);
            Table1{4}=num2str(TeMax(k),4);
            Table1{5}=num2str(NeOpt(k),4);
            Table1{6}=num2str(ErreurOptimal(k),4);
            Table1{7}=num2str(ChoixErreur);
            Table1{8}=num2str(exposant);
            Table1{9}=num2str(ChoixHautePression);
            Table1{10}=num2str(ChoixAutoabs);
            Table1{11}=Commentaire;
            fprintf(fileID,formatSpec,Table1{:});

                %% Écriture des raies utilisées dans un fichier
    Table2{1}=char(fileNames(k));
    Table2{2}=num2str(nbrderaies(k),2);
    for i=1:41
    Table2{i+2}=num2str(FitCorrige(k,i),1);
    end
    fprintf(fileID2,formatSpec2,Table2{:});
            
end

disp('Comparaison effectuée et température trouvée')
fclose ('all');



%% Obtention des résultats pour tous les fichiers
if TeGrapheFinaux==1
    TeGraphesFinaux_TRG(NumFichier,NeOpt,NeMin,NeMax,TeOptimal,TeMin,TeMax,I_exp);
end
close(h) 
clear NumFichier i j k l m fileID formatSpec fileID2 formatSpec2 IntExp ...
     TeGraph2 TeGraph1 AllIntegral TeGrapheFinaux nbrderaie  Data Table2 
