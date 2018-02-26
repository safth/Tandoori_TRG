function Calcul_Te_TRG_ViaInterface(Overwrite,InfoGraphes,InfoTe,InfoNe,filetype,Tg,P_tot,longueur,FctRep,ChoixErreur,ChoixHautePression,ChoixAutoabs,Commentaire,flow,exposant,Dimension,sig_longueur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CECI EST LE CODE PRINCIPAL POUR CALCULER TE � PARTIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  D'UN SPECTRE CONTENANT UNE TRACE DE GAZ RARES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
warning('off') % ne pas affichier les avertissements lors du calcul des 1s (1s1 qui nexiste pas fait un NAN)
%% ---------- Notes d'utilisation ----------
%====== Les sous fonctions utilis�es sont: ======
% - Load l'indication de quels raies fitter "nearkrxe.txt"
% - Fit les raie des spectres avec la fonction IntensiteGaussMax_TRG.m
% - Calcul les flow effectifs avec TePumpingSpeed_TRG.m

% Loop Gaz
%   Selectionne les taux de r�action appropri�s (pr�calcul� avec CalculTauxDeReactions_TRG.m).
%   Loop Ne
%       Loop Te
%           Calcul des intensit�s th�oriques pour tout Ne Te (TeIntTheo_TRG.m)
%       End
%   End 
% End
% 
% TeCalculErreur_TRG calcul une %STD pond�r� par des poids trouv� par un calcul de propagation d'erreur et minimise I_theo/I_exp pour trouver Te et Ne
% Output les Valeurs � l'optimal
% Affiche les graphs si s�lectionn�s

% Tous ces fichiers ainsi que les sections efficaces sont dans le sous-dossier "Complements"
%------------------------------------------
%%  ###########################################################################################
%%  ######################################### INPUTS ##########################################
%%  ###########################################################################################

% Range des temp�ratures �lectroniques scann�es (eV)
    Te =InfoTe(1):InfoTe(3):InfoTe(2);   
    
% Range de densit� �lectronique scann�e (m^-3) en log space
   Ne =logspace(log10(InfoNe(1)/InfoNe(2)),log10(InfoNe(1)*InfoNe(2)),InfoNe(3));  %logspace(X,Y,Z)= 10^X � 10^Y, avec Z valeurs entre les deux 
    
% Choix d'affichage des diff�rents graphes (=0 n'affice pas, =1 affiche)
    GraphExp=InfoGraphes(1);         %Graphe de l'extraction des intensit� exp�rimentales des raies
    TeGraph1=InfoGraphes(2);         %Output toute en .mat (c'est pas un graph)
    TeGraph2=InfoGraphes(3);         %graphe des m�canismes pop depop � l'optimum
    TeGraph3=InfoGraphes(5);         %Graphes du fit des raies, 3D de l'erreur et STD
    TeGrapheFinaux=InfoGraphes(4);   %Graphes de l'�volution de NeOptimal et TeOptimal en fonction des fichiers
    
    %load le choix de transfert d'excitation 1=oui, 0=non
    global ChoixTransEx;

%% ###########################################################################################
%% ###################################### D�BUT DU CODE ######################################
%% ###########################################################################################
%% ============== Raies d'�mission utilis�es dans le calcul ==============
     if flow(3) > 10*flow(4) % si on a 10X plus dargon que de Krypton
         nearkrxe = 'nearkrxeAr.txt'; %met des doublet pis des triplet � des endroits diff�rents!!
     else
         nearkrxe = 'nearkrxe.txt';
     end
    tempfile=fopen(nearkrxe); %Longueurs d'onde et "Doublet"
    data = textscan(tempfile,'%s %s','delimiter',' ','headerlines',6);
    E=str2double(data{1,1});  % Variable avec les 41 raies � fiter en nm
    Doublet=str2double(data{1,2}); %variable qui dit si il faut fiter une gaussienne,2 ou 3 gaussienne
    fclose(tempfile);
    clear data
    
      
%% ###########################################################################################
%% ####################### D�BUT DE L'ANALYSE EXP�RIMENTALE DES SPECTRES #####################
%% ###########################################################################################

%% ========== Identification de tous les spectres � analyser dans le r�pertoire s�lectionn� (selon l'extension)========== 
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

%% ========== Pr�allocation d'espace ==========
NumFichier      =zeros(length(fileNames),1);
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


%% regarde si on a pas d�j� extrait les intensit�s th�oriques.
 disp ('Extraction de l''intensit� des raies des fichiers...')
    h = waitbar(0,'Fit des spectres') ; % affichage d'une WaitBar qui montre la progression du code.
    
% si non, on les extrait et save le resultats dans I_exp.mat
%JE VEUX QUE LE CODE FIT SEULEMENT LES RAIES QUI MANQUE ET QUIL RAJOUTE LE
%R�SULTAT DANS LA STRUCTURE. SINON, LE TE_ALL N'A PAS LE N�ON ET IL
%RECALCULE TOUT. DONC SEULEMENT RAJOUTER LE NEON QUAND JE FAIT TAIL.

if exist('I_exp.mat','file')==2 % si le fichier existe
    load 'I_exp.mat';
    % regarde si le overwrite match pour tous les fichiers.
    % si une raies n'a pas �t� fitt� quelque part, le overwrite va pas
    % match�
    Fit = Struc.Fit;
    for k=1:length(fileNames)
        FitCorrige(k,:)=Fit(k,:).*Overwrite;
        I_exp(k,:) = Struc.I_exp(k,:).*FitCorrige(k,:);
        sig_I_exp = Struc.sig_I_exp.*FitCorrige(k,:);
    end

else     
    temp_Overwrite = ones(1,41)
       
    %% ============= Extraction de l'intensit� des raies de chaque spectre (i.e. fichier) =============
   for k=1:length(fileNames)
        Fichier=fileNames{k}          %Affiche � l'�cran quel fichier est actuellement analys�
        NumFichier(k)=k;              %Vecteur servant aux graphiques finaux de Te et Nm pour chaque fichier

        %% ============= Extraction des donn�es brutes (pour fichier .csv) =============
        if filetype == 0
            tempdata=dlmread(fileNames{k});         % Extraction des data
            lambdaTe=tempdata(:,1);                 % Longueurs d'ondes des mesures
            intTe=tempdata(:,2);                    % Intensit� des mesures          
            clear tempdata                          % On a plus besoin de tempdata  
            Int_time_correction=1;
        end   
        %% ============= Extraction des donn�es brutes de Morgane (pour fichier .csv) =============
        if filetype == 0.1 %fuck it
            tempdata=dlmread(fileNames{k});         % Extraction des data
            lambdaTe=tempdata(:,1);                 % Longueurs d'ondes des mesures
            intTe=tempdata(:,3);                    % Intensit� des mesures POUR MORGANE 
            clear tempdata                          % On a plus besoin de tempdata
            Int_time_correction=1;
        end    
        %% ============= Extraction des donn�es brutes (pour fichier .trt) =============
        if filetype == 1
            file1=fopen(fileNames{k});
            tempdata=textscan(file1, '%f%f','delimiter',';', 'headerlines',8);  % Extraction des datas 
            lambdaTe=tempdata{1};                   % Longueurs d'ondes des mesures
            intTe=tempdata{2};                      % Intensit� des mesures
            clear tempdata                          % On a plus besoin de tempdata   
            Int_time_correction=1;
        end    
        %% ============= Extraction des donn�es brutes (pour fichier .trtx) =============
        if filetype == 2
            %function qui prend les temp d'integration des deux spectros utilis�s (pour le 500-900, merge, 2 temps d'integration diff�rents)
            [Int_time_correction] = IntTime_trtx(fileNames{k});  
            file1=fopen(fileNames{k});
            tempdata=textscan(file1, '%f%f','delimiter',';', 'headerlines',10);  % Extraction des datas
            lambdaTe=tempdata{1};                   % Longueurs d'ondes des mesures
            intTe=tempdata{2};                      % Intensit� des mesures
            clear tempdata                          % On a plus besoin de tempdata 
        end  
        %% ============= Extraction des donn�es brutes (pour fichier .txt) =============
        if filetype == 3
            file1=fopen(fileNames{k});
            tempdata=textscan(file1, '%f%f','delimiter',' ');  % Extraction des datas
            lambdaTe=tempdata{1};                   % Longueurs d'ondes des mesures
            intTe=tempdata{2};                      % Intensit� des mesures
            clear tempdata                          % On a plus besoin de tempdata 
            Int_time_correction=1;
        end     
        %% ============= Extraction des donn�es brutes (pour fichier .TXT (Avantes8)) =============
        if filetype == 4
            file1=fopen(fileNames{k});
            tempdata=textscan(file1, '%s\t%s\t%s\t%s','delimiter',';', 'headerlines',9);  % Extraction des datas 
            lambdaTe=tempdata{1}  ;                 % Longueurs d'ondes des mesures
            lambdaTe=str2double(strrep(lambdaTe,',','.'));
            intTe=tempdata{2};                      % Intensit� des mesures
            intTe=str2double(strrep(intTe,',','.'));
            clear tempdata                          % On a plus besoin de tempdata 
            Int_time_correction=0.1;
        end

        %% ============= Correction par la fonction de r�ponse si n�cessaire =============
        if FctRep~= 0
           intTe=intTe./FctRep(:,2);    % Ireel=Imesure/Fct de reponse
        end

        %% ============= Correction pour un d�callage �ventuel en longueur d'onde fix� sur la 763nm =============
        lambdaTe(769:end)=lambdaTe(769:end)+0.5;
        w=1;
        for j=1:length(lambdaTe)
            if abs(lambdaTe(j)-763.51)<2
                IntTemp(w)=intTe(j);
                LambdaTemp(w)=lambdaTe(j);
                w=w+1;
            end
        end

        Shift=763.51-LambdaTemp(IntTemp==max(IntTemp));
        lambdaTe=lambdaTe+Shift;
        clear Shift w j IntTemp LambdaTemp
        %% ============= Obtention de l'intensit� des diff�rentes raies (au max du pic avec fit Gaussien=============
        [IntPeak,sig_IntPeak,Fit] = IntensiteGaussMax_TRG(E,lambdaTe,intTe,GraphExp,Doublet,temp_Overwrite); 
        %La variable Fit dit quelles raies ont effectivement �t� fitt�es

        %SaveIntPeak(k,:) =IntPeak;
        %% ============== Correction pour les raies qui ne doivent pas �tre consid�r�e =============
        temp_FitCorrige(k,:)=Fit.*temp_Overwrite;
        FitCorrige(k,:)=Fit.*Overwrite;
        temp_I_exp(k,:)=IntPeak'.*temp_FitCorrige(k,:);
        I_exp(k,:)=IntPeak'.*FitCorrige(k,:);
        %correction pour le integration time du 500-700, les raies 1 � 6 sont sur le 500-700
        I_exp(k,1:6)=I_exp(k,1:6)*Int_time_correction; 
        temp_I_exp(k,1:6)=temp_I_exp(k,1:6)*Int_time_correction; 
        % meme chose pour l'incertitude
        sig_I_exp(k,:)=sig_IntPeak'.*FitCorrige(k,:);
        temp_sig_I_exp(k,:)=sig_IntPeak'.*temp_FitCorrige(k,:);
        sig_I_exp(k,1:6)=sig_I_exp(k,1:6)*Int_time_correction; %correction pour le integration time du 500, les raies 1 � 6 sont sur le 500-700
        temp_sig_I_exp(k,1:6)=temp_sig_I_exp(k,1:6)*Int_time_correction;
        clear  lambdaTe intTe Fit IntPeak sig_IntPeak
        %check si l'intensit� de la raie est plus grande que son incertitude.
        %Si oui, on l'enleve.
        for m=1:size(I_exp,2)
            if I_exp(k,m) < sig_I_exp(k,m)
                I_exp(k,m) = 0;
                sig_I_exp(k,m) = 0;
                FitCorrige(k,m) = 0;
                disp('Intensit� fitt� plus petit que l incertitude.')
            end
        end
           waitbar(k/length(fileNames),h) %la waitbar progresse apres chaque spectre fitt�.
    end %fin de la boucle sur les diff�rents spectres
    %% Sauvegarde des I_exp dans un structure. Si on a deja fit� une fois et qu'on reroule le code, pas besoins de refiter.
    %%
    
    Struc.I_exp = temp_I_exp;
    Struc.sig_I_exp = temp_sig_I_exp;
    Struc.Fit = temp_FitCorrige;
    save('I_exp.mat','Struc');
    
end



disp ('Intensit�s extraites des fichiers')
clear Doublet filetype
% save('SaveIntPeak.mat','SaveIntPeak')


    %% =============  Calcul les pression partielle Par la diffusion diff�rentes du m�lange de gaz ==========
    %calcul des facteurs de corrections et l'affiche, ainsi que la pression
    %calcul�. Si la pression calcul� est tr�s diff�rente, vos flux doivent
    %�tre corrig�s. par exemple si vous pass� par une needle valve, le flux est beaucoup plus petit que celui inject�.
    [Pump,Pression_Calc] = TePumpingSpeed_TRG(P_tot,flow,Dimension)
    P_tot=(133.3224e-3)*P_tot;  %mtorr � pascal
    flow=flow.*Pump;            %on applique les corrections aux Flows
    flow_tot=sum(flow);         %flow total
    P = (P_tot.*flow)/flow_tot; %Pressions partielles de gaz
    % Calcul de la densit� de neutres � l'�tat fondamental avec p=n*k*Tg [m-3]
    % ng(gaz) car P est un vecteur avec les pressions partielles des gaz utilis�s
    ng=(P)./(1.38064852*10^(-23)*Tg);


%% ###########################################################################################
%% ########################### D�BUT DE L'ANALYSE TH�ORIQUE DU CODE ##########################
%% ###########################################################################################
disp ('Calcul de l''intensit� th�orique des raies...')


    %% Pr�allocation d'espace
    I_theo      =zeros(length(Ne),length(Te),length(E)); %les 41 raies dans le meme fichier
    sig_I_theo  =zeros(length(Ne),length(Te),length(E)); %les 41 raies dans le meme fichier
    Thetaij     =ones(length(Ne),length(Te),length(E));  %autoabs des 41 raies
    Doppler     =zeros(length(Ne),length(Te),length(E)); %on s'en sert pas
    VanDerWaals =zeros(length(Ne),length(Te),length(E)); %on s'en sert pas
    Resonant    =zeros(length(Ne),length(Te),length(E)); %on s'en sert pas

    Gains2p    =zeros(5,length(Ne),length(Te),10,5); %(5Xgaz,Ne,Te,10x2p,5x1s)
    Pertes2p   =zeros(5,length(Ne),length(Te),10,3); %(5Xgaz,Ne,Te,10x2p,5x1s)
    densite2p     =zeros(5,length(Ne),length(Te),10); %(5Xgaz,Ne,Te,10x2p) 
    ContributionFond    =zeros(5,length(Ne),length(Te),10); 
    densite1s     =zeros(5,length(Ne),length(Te),5);  %(5Xgaz,Ne,Te,5x2p) 
    sig_densite1s     =zeros(5,length(Ne),length(Te),5);

    Gains1s    =zeros(5,length(Ne),length(Te),5,5);  %(5Xgaz,Ne,Te,5x2p,5Gains) 
    Pertes1s   =zeros(5,length(Ne),length(Te),5,6);  %(5Xgaz,Ne,Te,5x2p,6Pertes) 
    
    
    %% =============== s�lection des taux de r�actions ===============
    %load des taux pr�calcul�s et asignation d'une variable global � chacun
        %initialise des globals
    global K_1s_2p K_gs_1s K_quench_1s K_neutral K_gs_2p Choix_Taux 
    
    % prend des .mat moins lourds et donc, plus rapide!!
    if Choix_Taux ==1 % si on est maxwellien (exposant=1)
        load M_Allrates1s_2p.mat 
        load M_AllratesGround_1s.mat 
        load M_AllratesQuenching_1s.mat
        load M_Allrates_neutral.mat
        if ChoixHautePression==0
            load M_AllratesGround_2p.mat    %AllratesGround_2p(#gas,2Px,Te)
        elseif ChoixHautePression==1 %si on est � Haute pression partielles d'Argon (>1mtorr)
            load M_AllratesGround_2p_HighP.mat    %AllratesGround_2p(#gas,2Px,Te)
        end 
        % prend les fichiers lourds avec tous les exposants si on ne prend pas exposant=1
    else
        load Allrates1s_2p.mat 
        load AllratesGround_1s.mat 
        load AllratesQuenching_1s.mat
        load Allrates_neutral.mat
        if ChoixHautePression==0
            load AllratesGround_2p.mat    %AllratesGround_2p(#gas,2Px,Te)
        elseif ChoixHautePression==1 %si on est � Haute pression partielles d'Argon (>1mtorr)
            load AllratesGround_2p_HighP.mat    %AllratesGround_2p(#gas,2Px,Te)
        end
    end

    % assigne les Taux aux variables globals
    if ChoixHautePression==0
        K_gs_2p = AllratesGround_2p;
    elseif ChoixHautePression==1 %si on est � Haute pression partielles d'Argon (>1mtorr)
        K_gs_2p = AllratesGround_2p_HighP;
    end
    K_1s_2p     =Allrates1s_2p;
    K_gs_1s     =AllratesGround_1s;
    K_quench_1s =AllratesQuenching_1s;
    K_neutral   =rates_neutral;
    
    
    %% ====== Boucles sur les gaz, �l�ments th�oriques du code. =======
    wait=0; % wwait bar initilis� � 0 pour ce calcul (a priori le plus long)
    global gaz_i gaz_f
for gaz=gaz_i:gaz_f %On fait le calcul th�orique de l'intensit� de raies pour Gaz_i � Gaz_f
    
      % pour mettre les valeurs dans des vecteur avec tous les gaz inclus un � la suite de l'autre.
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
       %progression de la waitbar
       waitbar(wait/(length(gaz_i:gaz_f)*length(Ne)),h,'Calcul th�orique: s�lection des Taux de r�actions') 
    
    %% ============== Boucle d'initialisation sur Te pour la s�lection des bons taux de r�action ==============
    % fonction qui prend les taux de r�action correspondant aux Te et au step de Te s�lectionn�s.
    [rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral] = TeTauxDeReaction_TRG(gaz,Te,ChoixHautePression,ChoixAutoabs,exposant);
    disp('Taux de r�action selectionn�s')

    %%  ============== Premi�re loop sur la densit� �lectronique: COMPTEUR: j ==============
    for j=1:length(Ne)
        % progresse la waitbar.
        waitbar(wait/(length(gaz_i:gaz_f)*length(Ne)),h,'Calcul th�orique') 
        wait=wait+1;
        %% Calcul du bilan de population: gains=pertes
                     % !! tout se fait ici!!
        % ================================================
        % ================================================
        %       LOOP sur Te dans cette fonction!!
       [densite_1s,sig_densite_1s,densite_2p,sig_density,Gains_2p,Pertes_2p,Emission,PopFond,energie_1s,energie_2p,Gains_1s,Pertes_1s]=TeIntTheo_TRG(gaz,Ne(j),Te,rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral,Tg,longueur,P,ChoixAutoabs,sig_longueur,densite1s(3,j,:,:),sig_densite1s(3,j,:,:));
        % ================================================
        % ================================================
       %% Extraction et mise en m�moire des donn�es pour chaque fichier. Pour chaque Ne et Te
       for i=1:length(Te)
           n=0;
            for m=index(1):index(2)                     % Pour chaque raie...
                n=n+1;% pour emission
                I_theo(j,i,m)=Emission(1,i,n);          % Intensit� th�orique
                sig_I_theo(j,i,m)=Emission(6,i,n);      % Intensit� th�orique
                Thetaij(j,i,m)=Emission(2,i,n);         % Escape Factor 
                Doppler(j,i,m)=Emission(3,i,n);         % �largissement Doppler
                VanDerWaals(j,i,m)=Emission(4,i,n);     % �largissement VanDerWaals
                Resonant(j,i,m)=Emission(5,i,n);        % �largissement Resonant
            end  
            for l=1:10 
                Gains2p(gaz,j,i,l,:)=Gains_2p(:,i,l);        %1=fond,2=Nm,3=coll2P,4=radtrap,5 = Transfert Radiatif 
                Pertes2p(gaz,j,i,l,:)=Pertes_2p(:,i,l);      % 1=rad 2=coll2p 3=coll1s % les deux dernier sont toujours = 0
                densite2p(gaz,j,i,l)=densite_2p(i,l);        % Densit� des 2p
                ContributionFond(gaz,j,i,l)=PopFond(i,l);    % Pourcentage de la population en comparant seulement le fondamental et les m�tastables    
                
            end
             for l=1:5   
                 Gains1s(gaz,j,i,l,:)=Gains_1s(:,i,l);    %1=fond,2=autoabs,3=mixing,4=updown,5=Transfert d'excitation
                 Pertes1s(gaz,j,i,l,:)=Pertes_1s(:,i,l);  %1=rad,2=2p,3=mixing,4=super�lastique,5=ionisation,6=quenching
                 densite1s(gaz,j,i,l)=densite_1s(i,l);
                 sig_densite1s(gaz,j,i,l)=sig_densite_1s(i,l);
             end
             energie1s(gaz,:)=energie_1s(:); %l'�nergie des niveaux selon le gaz.
             energie2p(gaz,:)=energie_2p(:); %l'�nergie des niveaux selon le gaz.
           
       end   %Fin Boucle Te
    end    %Fin Boucle Ne
end %fin boucle sur les gaz

clear Allrates1s_2p AllratesGround_2p AllratesGround_1s AllratesQuenching_1s rates_neutral position
clear Emission Gains Pertes density PopFond
disp('Intensit�s th�oriques calcul�es')

%% ============== Cr�ation des fichiers textes contenant les infos pertinentes pour chaque fichier analys� ==============
if exist('Te.out', 'file')==0 % n'�crit pas le Header si il est deja fait
    fileID = fopen('Te.out','a');               
    formatSpec = '%-22s\t %s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s\r\n';  
    fprintf(fileID,formatSpec,'------NomFichier------','TeOpt','TeMin','TeMax','NeOpt   ','% STD','Poids','exp','HighP','AutoAbs','Transfert','Commentaire');
else
    formatSpec = '%-22s\t %s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s\r\n';  
    fileID = fopen('Te.out','a');   
end

fileID2 = fopen('Raies Utilis�es.out','a');               
formatSpec2 = '%-22s \t%-8s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s \t%-3s\r\n';  
fprintf(fileID2,formatSpec2,'------NomFichier------','-#raies-','585','640','667','696','706','714','727','738','750','751','763','794','800','801','810','811','826','840','842','852','866','758','760','768','769','785','810','811','819','829','850','877','892','788','820','823','828','834','881','895','904');

fileID3 = fopen('Nm_Argon.out','a');   
formatSpec3 = '%.4e\t %e\t %.4e\t %e\r\n';

%% ###########################################################################################
%%  ####################### D�BUT DE LA COMPARAISON TH�ORIE-EXP�RIENCE #######################
%%  ##########################################################################################
disp('Comparaison exp�rience-th�orie...')   

%% ========== Comparaison des spectres exp�rimentaux un � un avec la th�orie ==========
for k=1:length(fileNames) %loop sur tous les spectres
    Fichier=fileNames{k}  %Affiche � l'�cran quel fichier est actuellement compar� 
    nbrderaies(k)=sum(FitCorrige(k,:));   % Calcul du nombre total de raies � analyser
    
    %% Boucles d'�valuation de l'erreur des ratio I_exp/I_simul aux temp�ratures et densites �lectroniques
    if nbrderaies(k) > 1             %V�rification s'il y a au moins 2 raies � ajuster      
            %% Calcul de la d�viation standard
            [STDErr,NeOptimal,TempOptimal,MoyenneOptimale(k),ErreurOptimal(k),FitCorrige(k,:)] = TeCalculErreur_TRG(ChoixErreur,Ne,Te,I_exp(k,:),I_theo,FitCorrige(k,:),P,sig_I_theo,sig_I_exp(k,:));

            %% Extraction des r�sultats
            NeOpt(k)=NeOptimal(1); NeMin(k)=NeOptimal(2); NeMax(k)=NeOptimal(3); PosNeMin(1)=NeOptimal(4);
            TeOptimal(k)=TempOptimal(1); TeMin(k)=TempOptimal(2); TeMax(k)=TempOptimal(3); PosTeMin(1)=TempOptimal(4);
            
            %% ecriture des densit� de metastable d'Argon � l'optimal. on peut faire les autres si on veut...
               fprintf(fileID3,formatSpec3,densite1s(3,PosNeMin(1),PosTeMin(1),2:5)); 
            %% Export de toute les fichiers...
            if TeGraph1==1 
                file(:,:)=Gains2p(3,PosNeMin(1),PosTeMin(1),:,:);
                file2(:,:)=Pertes2p(3,PosNeMin(1),PosTeMin(1),:,:);
                file3(:,:)=Gains1s(3,PosNeMin(1),PosTeMin(1),:,:);
                file4(:,:)=Pertes1s(3,PosNeMin(1),PosTeMin(1),:,:);
                file5(:,1)=I_theo(PosNeMin(1),PosTeMin(1),:);
                file5(:,2)=I_exp(1,:);
                file5(:,3)=1-Thetaij(PosNeMin(1),PosTeMin(1),:);
                save('Gains2P.mat','file')
                save('Pertes2P.mat','file2')
                save('Gains1s.mat','file3')
                save('Pertes1s.mat','file4')
                save('I_th-I_exp-AutoAbs.mat','file5')
            end
         
             %% Obtention des graphiques � l'optimum si d�sir� dans les inputs
           if TeGraph2==1
                TeGraphes2_TRG(Ne,Te,STDErr,I_theo(PosNeMin(1),PosTeMin(1),:),sig_I_theo(PosNeMin(1),PosTeMin(1),:),Thetaij(PosNeMin(1),PosTeMin(1),:),I_exp(k,:),sig_I_exp(k,:),FitCorrige(k,:),ChoixErreur,MoyenneOptimale(k),P,E,ChoixAutoabs);
           end
             %% Obtention des graphiques � l'optimum si d�sir� dans les inputs
           if TeGraph3==1
               TeGraphes3_TRG(densite2p(:,PosNeMin(1),PosTeMin(1),:),Gains2p(:,PosNeMin(1),PosTeMin(1),:,:),Pertes2p(:,PosNeMin(1),PosTeMin(1),:,:),ContributionFond(:,PosNeMin(1),PosTeMin(1),:),energie1s,energie2p,densite1s(:,PosNeMin(1),PosTeMin(1),:),Gains1s(:,PosNeMin(1),PosTeMin(1),:,:),Pertes1s(:,PosNeMin(1),PosTeMin(1),:,:),ng);
           end
    
            
    else %S'il n'y a pas 2 raies � comparer, on met tout � 0
        TeOptimal(k)=0;
        TeMin(k)=0;
        TeMax(k)=0;
        NeOpt(k)=0;
        NeMax(k)=0;
        NeMin(k)=0;

      disp('Nombre de raies � comparer insuffisant. Tous les param�tres ont �t� mis � z�ro. You Died.')
    end %fin boucle nbr de raies minimale
    
      waitbar(k/length(fileNames),h,'Calcul d''erreur')   %on progresse la waitbar
    clear IntensiteSD metSD tempSD Emission Gains Pertes density PosTeMin PosN1s2Min SEGD
end %Fin boucle fichiers / fin du calcul d'erreur

            %% �criture des param�tres optimaux dans le fichier Te.out
            disp('�criture des r�sultats')
for k=1:length(fileNames)
            formatSpec = '%-22s\t %-5s \t%-5s \t%-5s \t%-10s \t%-6s \t%-5s \t%-5s \t%-5s \t%-10s \t%-10s \t%-20s\r\n';  
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
            Table1{11}=num2str(ChoixTransEx);
            Table1{12}=Commentaire;
            fprintf(fileID,formatSpec,Table1{:});

                %% �criture des raies utilis�es dans le fichier Raies Utilis�es.out
            Table2{1}=char(fileNames(k));
            Table2{2}=num2str(nbrderaies(k),2);
            for i=1:41
                Table2{i+2}=num2str(FitCorrige(k,i),1);
            end
            fprintf(fileID2,formatSpec2,Table2{:});
            
end

disp('Comparaison effectu�e et temp�rature trouv�e')
fclose ('all');


%% Obtention des r�sultats pour tous les fichiers
if TeGrapheFinaux==1
    TeGraphesFinaux_TRG(NumFichier,NeOpt,NeMin,NeMax,TeOptimal,TeMin,TeMax,I_exp);
end
close(h) 
clear NumFichier i j k l m fileID formatSpec fileID2 formatSpec2 IntExp ...
     TeGraph2 TeGraph1 TeGrapheFinaux nbrderaie  Data Table2 
 
 % bravo, vous etes maintenant un pro de la TRG-OES SIGNEZ
 % Simon Boivin
