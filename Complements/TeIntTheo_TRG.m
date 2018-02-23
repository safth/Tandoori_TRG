function [densite1s,sig_densite1s,densite2p,sig_densite2p,Gains2p,Pertes2p,Emission,ContributionFond2p,energie1s,energie2p,Gains1s,Pertes1s] = TeIntTheo_TRG(gaz,ne,Te,rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral,Tg,longueur,P,ChoixAutoabs,sig_longueur,nm_Ar,sig_nm_Ar)

%% =====================================================================================================
%% ================ Cette fonction calcule les intensit� th�orique de raies pour tous ==================
%% =============== les Te � pour une densit� �lectronique donn�e et pour un gaz donn�e  ================
%% =====================================================================================================

%% ====== INFOS SUR LA FONCTION ======
%La variable gaz indique pour quel gaz on cherche les intensit�s de raie
%La variable ne est la densit� �lectronique
%Le vecteur Te contient toutes les temp�ratures qui doivent �tre analys�es au cours de cette fonction
%La variable Tg contient la temp�rature du gaz    
%La variable longueur est la longueur du parcours optique dans le plasma pour le calcul d'auto-absorption
%Le vecteur P contient les pressions partielles
%------------------------------------------------------------------------------------
    %POUR INFO:
    %rateGround_2p = taux de r�action pour l'excitation des niveaux 2p par impact �lectronique sur le fondamental
    %rateGround_1s = Taux de r�action pour l'excitation du fondamental vers les 1s par impact �lectronique
    %rate1s_2p = taux de r�action pour l'excitation des niveaux 2p par impact �lectronique sur les niveaux 1s
    %rateQuenching = taux de r�action pour l'�change entre les 1s par impact �lectronique avec l'ionisation et le super�lastique en 6 et 7iem position
    %rateNeutral = Quenching des m�tastables vers le fondamental par les neutres ex : (3,4) = neutre Argon sur m�tastables de Krypton
    %En2px_1s = taux de transfert d'�nergie par collision avec les neutre des niveaux 2p (menant � 1s)
    %En2px_2py = taux de transfert d'�nergie par collision des niveaux 2p
    %poids2p = poids statistique des niveaux 2p
    %poids1s = poids statistique des niveaux 1s
    %LambdaTheo = longueurs d'onde des transition radiatives 2p-1s
    %Aij = les coefficients de desexcitation radiative spontan�e
    %energie1s = �nergie des niveaux 1s
%% ===== nomenclature =====
%les vecteur ont comme premi�re entr�e le num�ro du gaz utilis�
% ex rate1s_2p(4,3,10,i) krypton, 1s3->2P10 pour Te(i)
% 1 = He (inexistant)
% 2 = N�on
% 3 = Argon
% 4 = Krypton
% 5 = X�non
% 13 = O2 
%% ====== D�BUT DE LA FONCTION =======

% R�cup�ration des donn�es utiles / CONSTANTE
[En2px_1s,En2px_2py,poids2p,poids1s,LambdaTheo2p,LambdaTheo1s,Aij2p,Aij1s,energie1s,energie2p,fji2p,fji1s,Br2p]=TeInitialisation_TRG(gaz);
    
    %% Pr�allocation d'espace
    %selon le nombre de raie observ� par gaz, les vecteurs n'ont pas la meme dimension
    if gaz==2
    index=2;
    end
    if gaz==3
    index=19;
    end
    if gaz==4
    index=12;
    end
    if gaz==5
    index=8;
    end
    
    GainFond           =zeros(length(Te),10);
    GainNm             =zeros(length(Te),10);
    GainColl2p         =zeros(length(Te),10);
    GainRadTrap        =zeros(length(Te),10);
    GainAr             =zeros(length(Te),10);
    PerteRadiatif      =zeros(length(Te),10);
    PerteColl2p        =zeros(length(Te),10);
    PerteColl1s        =zeros(length(Te),10);
    ContributionFond2p =zeros(length(Te),10);   
    densite2p          =zeros(length(Te),10);
    densite1s          =zeros(length(Te),5);
    GainFond_1s        =zeros(length(Te),5);
    GainAutoAbs_1s     =zeros(length(Te),5);
    GainMixing_1s      =zeros(length(Te),5);
    GainUpdown_1s      =zeros(length(Te),5);
    GainPopAr_1s       =zeros(length(Te),5);
    PerteVers2p_1s     =zeros(length(Te),5);
    PerteRad_1s        =zeros(length(Te),5);
    PerteMixinxg_1s    =zeros(length(Te),5);
    PerteSuper_1s      =zeros(length(Te),5);
    PerteIonisation_1s =zeros(length(Te),5);
    PerteNeutre_1s     =zeros(length(Te),5);
    GainsRes_1s        =zeros(length(Te),5);
    GainsMetas_1s      =zeros(length(Te),5);
    sig_gain           =zeros(length(Te),5);
    sig_perte          =zeros(length(Te),5);
    gain               =zeros(length(Te),5);
    perte              =zeros(length(Te),5);
    sig_densite1s      =zeros(length(Te),5);
    sig_gain_2P        =zeros(length(Te),10);
    gain_2P            =zeros(length(Te),10);
    sig_densite2p      =zeros(length(Te),10);

    I_theo          =zeros(length(Te),index);
    sig_I_theo      =zeros(length(Te),index);
    Theta           =zeros(length(Te),index);
    Doppler         =zeros(length(Te),index);
    VanDerWaals     =zeros(length(Te),index);
    Resonant        =zeros(length(Te),index);
    
    Thetaij         =zeros(10,5);
    sig_Thetaij     =zeros(10,5);
    Thetaij_1s      =zeros(1,5);
    Dopp            =zeros(10,5);
    VDWaals         =zeros(10,5);
    Res             =zeros(10,5);
        
    n1sX            =zeros(1,5);
    
    %% Correction des taux collisionnels en temp�rature
    En2px_1s=En2px_1s*(sqrt(Tg/300)); %on pourrait faire ca pour le qucnhing aussi!! il depend de sqrt(Tg)
    En2px_2py=En2px_2py*(sqrt(Tg/300));    

%% Calcul de la densit� de neutres � l'�tat fondamental avec p=n*k*Tg [m-3]
   %ng(gaz) car P est un vecteur avec les pressions partielles des gas utilis�s
    ng=(P)./(1.38064852*10^(-23)*Tg);   
   
%% ================== Boucle sur la temp�rature �lectronique: COMPTEUR: t =================
for t=1:length(Te)

%% ##################################################################################################  
%  ##################################################################################################
%  ############# D'abord on calcule la densit� des niveaux 1s via un bilan gains-pertes #############
%  ##################################################################################################
%  ##################################################################################################
%% Calcul de l'�largissement et du pi�geage optique pour les transitions radiatives 1s-fondamental
    for i=1:5
        if Aij1s(i)~=0
            if ChoixAutoabs==1
               [Thetaij_1s(i)]=TeEscapeFactorDOPE_TRG(gaz,LambdaTheo1s(i),poids1s(i),1,Aij1s(i),ng(gaz),longueur,Tg,0,sig_longueur);        
            else
               Thetaij_1s(i) = 1; % on met l'autoabsorption � 0 si on a coch� de l'ignorer (autoabs = 1-thetaij)
            end
        end
    end 
%% M�canismes ind�pendants de la densit� des 1s
[PopFond]= TePopulation_Metastable_TRG(ng(gaz),ne,rateGround_1s(:,t));
    
%% M�canismes d�pendant de la densit� des 1s    
[Depop2p,DepopRadFond,PopAutoAbs,DepopMix,PopMix,DepopSuperelestique,DepopIonisation,DepopNeutre,PopUpDown,PopAr_1s] = TeDepopulation_Metastable_TRG(t,gaz,ng,ne,Aij1s,Thetaij_1s,rate1s_2p,rateQuenching,rateNeutral,Br2p,nm_Ar);

%% Obtention de la densit� des niveaux pour avoir l'etat stationnaire depopulation=population
% les gains dans "depop" sont n�gatifs. (pop=depop). Certains gains sont dans Depop puisqu'ils d�pendent de la denit� des 1s. on les mets alors n�gatifs
depopulation1s=Depop2p+DepopRadFond-PopAutoAbs+DepopMix-PopMix+DepopSuperelestique+DepopIonisation + DepopNeutre - PopUpDown;
population1s=PopFond+PopAr_1s;

n1sX=linsolve(depopulation1s,population1s);
%n1sX(2)=1e1; %negliger les r�sonnant
%n1sX(4)=1e1; %negliger les r�sonnant
densite1s(t,:)=n1sX;
densite1s(isnan(densite1s)) = 0 ; %enleve le NaN en 1s1


%% ============== Suivi de l'influence des processus de peuplement et de d�peuplement des niveaux 1s ==============   

 GainFond_1s(t,:)      =  100*PopFond                       ./(PopFond + PopAr_1s + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainAutoAbs_1s(t,:)   =  100*PopAutoAbs*densite1s(t,:)'    ./(PopFond + PopAr_1s + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainMixing_1s(t,:)    =  100*PopMix*densite1s(t,:)'        ./(PopFond + PopAr_1s + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainUpdown_1s(t,:)    =  100*PopUpDown*densite1s(t,:)'     ./(PopFond + PopAr_1s + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainPopAr_1s(t,:)     =  100*PopAr_1s                      ./(PopFond + PopAr_1s + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 PerteVers2p_1s(t,:)     =  100*Depop2p*densite1s(t,:)'             ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteRad_1s(t,:)        =  100*DepopRadFond*densite1s(t,:)'        ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteMixinxg_1s(t,:)    =  100*DepopMix*densite1s(t,:)'            ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteSuper_1s(t,:)      =  100*DepopSuperelestique*densite1s(t,:)' ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteIonisation_1s(t,:) =  100*DepopIonisation*densite1s(t,:)'     ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteNeutre_1s(t,:)     =  100*DepopNeutre*densite1s(t,:)'         ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');

 
 %% Calcul de l'importance des r�sonant sur les densit�s de m�tstables.
%      if gaz==3
%          GainsRes_1s(t,:)     =  100*(PopMix(:,2:2:4)+PopUpDown(:,2:2:4))*densite1s(t,2:2:4)'       ./(PopFond + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
%          GainsMetas_1s(t,:)   =  100*(PopMix(:,3:2:5)+PopUpDown(:,3:2:5))*densite1s(t,3:2:5)'       ./(PopFond + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
%          GainFond_1s(t,:);
%      end
%% =========================================================== 
%  =========================================================== 
%  ================ Calcul de l'incertitude ==================
%  =======================sur les 1s ========================= 
%  ===========================================================  
%  =========================================================== 
%calcul des incertide sur les taux de r�actions. On r�soud par m�thode
%it�rative puisque ce n'est plus lin�aire.
[sig_PopFond]= TePopulation_Metastable_TRG(ng(gaz),ne,sig_rateGround_1s(:,t));
[sig_Depop2p,sig_DepopRadFond,sig_PopAutoAbs,sig_DepopMix,sig_PopMix,sig_DepopSuperelestique,sig_DepopIonisation,sig_DepopNeutre,sig_PopUpDown,sig_PopAr_1s] = TeDepopulation_Metastable_TRG(t,gaz,ng,ne,Aij1s,Thetaij_1s,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral,Br2p,sig_nm_Ar);

%initialisation de variables
iter=0; % boucle qui compte le nombre d'it�ration
residu=10e21; %start le r�sidu � quelquechose de tr�s gros.
gain(t,:)=0;
sig_gain(t,:)=0;
perte(t,:)=0;
sig_perte(t,:)=0; 

% boucle tant que le r�sidu est n'est pas de 1/1000
while abs(residu) > sig_densite1s(t,5)/1000 %boucle car 1s3 d�pend de 1s5 non lin�airement et vice versa
for j=2:5 % boucle sur les 1s � r�soudre
    for i=2:5 % boucle sur les autres 1s
        if j~=i 
          gain(t,j)     = gain(t,j)     + densite1s(t,i)*PopMix(j,i) + densite1s(t,i)*PopUpDown(j,i);      
          sig_gain(t,j) = sig_gain(t,j) + densite1s(t,i).^2*(sum(sig_PopMix(j,i) + sig_PopUpDown(j,i))).^2 + sig_densite1s(i).^2*(sum(PopMix(j,i) + PopUpDown(j,i))).^2; % sqrt � la fin
        end
    end
    perte(t,j)        =    Depop2p(j,j)+DepopRadFond(j,j)+DepopMix(j,j)+DepopSuperelestique(j,j)+DepopIonisation(j,j)+DepopNeutre(j,j);
    sig_perte(t,j)    =    sqrt(sig_Depop2p(j,j)^2+sig_DepopRadFond(j,j)^2+sig_DepopMix(j,j)^2+sig_DepopSuperelestique(j,j)^2+sig_DepopIonisation(j,j)^2+sig_DepopNeutre(j,j)^2);
    gain(t,j)         =    gain(t,j)+PopFond(j) + PopAr_1s(j);
    sig_gain(t,j)     =    sqrt(sig_gain(t,j) + sig_PopFond(j).^2 + sig_PopAr_1s(j).^2);
end    
clear i j
iter=iter+1;
%calcul du r�sidu entre 2 it�rations succesives
residu = sig_densite1s(t,5) - densite1s(t,5).*sqrt((sig_gain(t,5)./gain(t,5)).^2 + (sig_perte(t,5)./perte(t,5)).^2);
%incertitude sur la Densit� de m�tsstables
sig_densite1s(t,:) =  densite1s(t,:).*sqrt((sig_gain(t,:)./gain(t,:)).^2 + (sig_perte(t,:)./perte(t,:)).^2);
end
sig_densite1s(isnan(sig_densite1s)) = 0 ; %enleve le NaN en 1s1
residu(isnan(residu)) = 0; %enleve le NaN en 1s1

%% met l'incertitude sur les r�sonannt  = � celle des m�tastables car trop grande!! � REVOIR

sig_densite1s(t,2)=sig_densite1s(t,3)*densite1s(t,2)/densite1s(t,3);
sig_densite1s(t,4)=sig_densite1s(t,5)*densite1s(t,4)/densite1s(t,5);
%100*sig_densite1s(t,:)./densite1s(t,:);

 %% affichage des incertitudes sur les 1s
% 
%     figure1=figure;
%     name = {'erreur He' 'erreur Ne' 'erreur Ar' 'erreur Kr' 'erreur Xe'};
%     set(figure1,'name',name{gaz},'numbertitle','off','Position',[200 600 1000 600])
%     subplot('position',[0.05 0.55 0.90 0.36]) 
%     sig_1s = 100*sig_densite1s(t,:)./densite1s(t,:);
%     bar2=bar(1:5,sig_1s,0.75,'stacked');
%     set(gca,'XTick',1:5,'Xdir','reverse')
%     ylim([0 50])
%     xlim([1.4 5.6]);
%     title('Incertitude sur les 1s')
%     xlabel('Niveau 1s_x','FontSize',11,'fontweight','bold');
%      
%% ###########################################################  
%  ###########################################################
%  ############# Ensuite on peut r�soudre les 2p #############
%  ###########################################################
%  ########################################################### 
%% Calcul de l'�largissement et du pi�geage optique pour les transitions radiatives 2p-1s
    for i=1:10  %Boucle sur les 2p
        for j=1:5   %Boucle sur les 1s
            if Aij2p(i,j)~=0
                if ChoixAutoabs==1
                   [Thetaij(i,j),Dopp(i,j),sig_Thetaij(i,j)]=TeEscapeFactorDOPE_TRG(gaz,LambdaTheo2p(i,j),poids2p(i),poids1s(j),Aij2p(i,j),n1sX(j),longueur,Tg,sig_densite1s(t,j),sig_longueur);  
                else
                   Thetaij(i,j) = 1; %si on n�glige l'autoabsorption
                end
            end
        end      
    end

    
    

%% Calcul des processus de population/depopulation d�pendant de la densit� des niveaux 2p (inclut: �mission et absorption radiative, processus collisionnel entre les 2p et vers les 1s)
    [DepopRad2p,PopRad2p,DepopColl2p,DepopColl1s,PopColl2p]=TeDepopulation_TRG(Aij2p,Thetaij,En2px_1s,En2px_2py,ng(gaz)); 
  
%% Calcul des processus de population des niveaux 2p via impact �lectronique sur le fondamental et les niveaux 1s
    [PopFond2p,PopNm2p,Pop1s,PopAr]=TePopulation_TRG(t,ng,ne,n1sX,rate1s_2p,rateGround_2p,gaz,nm_Ar);

%% R�organisation des matrices: Les fonction flipud et rot90 sont utilis�es car � priori les matrices contiennent en premi�re position 2p10, 
%en 2e position 2p9, etc. Elles permettent donc de mettre ce qui correspond � 2p1 en premi�re position, 2p2 en seconde position, etc
    DepopRad2p=rot90(DepopRad2p,2); %Le ,2 est parce qu'on fait une rotation de 2x 90deg
    PopRad2p=rot90(PopRad2p,2);
    DepopColl2p=rot90(DepopColl2p,2);
    DepopColl1s=rot90(DepopColl1s,2);
    PopColl2p=rot90(PopColl2p,2);
    PopFond2p=flipud(PopFond2p);
    PopNm2p=flipud(PopNm2p);
    PopAr=flipud(PopAr);
    
    %% Obtention de la densit� des niveaux pour avoir l'etat stationnaire depopulation=population
    depopulation=DepopRad2p-PopRad2p+DepopColl2p+DepopColl1s-PopColl2p;
    %depopulation=DepopRad2p-PopRad2p; % on prend juste les trans Radiative
    population=PopNm2p+PopFond2p+PopAr;
    
    densite2p(t,:)=linsolve(depopulation,population);              %depopulation*density=population
    
    
    %% ============== Suivi de l'influence des processus de peuplement et de d�peuplement des niveaux 2p ==============   
    GainFond(t,:)=   100*PopFond2p                 ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p+PopAr);
    GainNm(t,:)=     100*PopNm2p                   ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p+PopAr);
    GainColl2p(t,:)= 100*PopColl2p*densite2p(t,:)' ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p+PopAr);    
    GainRadTrap(t,:)=100*PopRad2p*densite2p(t,:)'  ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p+PopAr);
    GainAr(t,:)  =   100*PopAr                     ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p+PopAr);
    
    PerteRadiatif(t,:)= (100*DepopRad2p*densite2p(t,:)')./((DepopRad2p+DepopColl2p+DepopColl1s)*densite2p(t,:)');
    PerteColl2p(t,:)=(100*DepopColl2p*densite2p(t,:)')./((DepopRad2p+DepopColl2p+DepopColl1s)*densite2p(t,:)');
    PerteColl1s(t,:)=(100*DepopColl1s*densite2p(t,:)')./((DepopRad2p+DepopColl2p+DepopColl1s)*densite2p(t,:)');
    
    %% En particulier, on regarde pour quels param�tres l'excitation du fondamental (ou � l'inverse via les 1s) domine
    ContributionFond2p(t,:)=100*PopFond2p./(PopFond2p+PopNm2p);
    
    %% =========================================================== 
    %  =========================================================== 
    %  ================ Calcul de l'incertitude ==================
    %  =======================sur les 2P ========================= 
    %  ===========================================================  
    %  =========================================================== 
    [sig_PopFond2p,sig_PopNm2p,sig_Pop1,sig_PopAr]=TePopulation_TRG(t,ng,ne,n1sX,sig_rate1s_2p,sig_rateGround_2p,gaz,sig_nm_Ar);
    [sig_DepopRad2p,sig_PopRad2p,sig_DepopColl2p,sig_DepopColl1s,sig_PopColl2p]=TeDepopulation_TRG(Aij2p,sig_Thetaij,En2px_1s,En2px_2py,ng(gaz));   %Consid�re desexc. rad. + transfert coll vers le haut et le bas

    %Aucune incertitude sur le Depop. Enfait, une petite sur les Aij du Xe
    %qui pourrait etre ajout�. Je ne crois pas que ca change vram
    %quelquechose surtout quelle est consid�r� sur l'intensit� des raies
for iter=1:3
    sig_gain_2P_a = zeros(10,1);
    sig_gain_2P(t,:) = 0;
    gain_2P(t,:) = 0;
        for j=1:10 % boucle sur les 2P � r�soudre
            for i=2:5 % Boucle sur les 1s
              sig_gain_2P_a(j) = sig_gain_2P_a(j) + (sig_densite1s(t,i)*Pop1s(j,i))^2;
            end
            sig_gain_2P(t,11-j) = sqrt(sig_gain_2P_a(j) + sig_PopFond2p(j)^2 + sig_PopNm2p(j)^2 + sig_PopAr(j)^2 + (densite2p(t,j)*sig_PopRad2p(j))^2 + (sig_densite2p(t,j)*PopRad2p(j))^2);%11-j car gain(1)=2p10 et gain(10)=2p1
            gain_2P(t,j) = PopFond2p(j)+ PopNm2p(j) + densite2p(t,j)*PopRad2p(j)+PopAr(j);

        end    
    clear i j
    %Valeur final de lincertitude sur les densit� de 2P
    %reste =  sig_densite2p(t,:) - densite2p(t,:).*(sig_gain_2P(t,:)./gain_2P(t,:))
    sig_densite2p(t,:) = densite2p(t,:).*(sig_gain_2P(t,:)./gain_2P(t,:));
    sig_densite2p(isnan(sig_densite2p)) = 0; %enleve les NaN
end

%% affichage des diff�rentes contributions � l'erreur
%  if IndexOffset>1
%     subplot('position',[0.05 0.07 0.90 0.36])
% 
%      erreur(:,1) = 100*sqrt(sig_gain_2P_a)./(sqrt(sig_gain_2P_a) + sig_PopFond2p + sig_PopNm2p)
%      erreur(:,2) = 100*sig_PopFond2p./(sqrt(sig_gain_2P_a) + sig_PopFond2p + sig_PopNm2p)
%      erreur(:,3) = 100*sig_PopNm2p./(sqrt(sig_gain_2P_a) + sig_PopFond2p + sig_PopNm2p)
%      sig_2P = 100*sig_densite2p./ densite2p;
%      
%  
%      
%     bar2=bar(1:10,sig_2P,0.75,'stacked');
%     set(gca,'XTick',1:10,'Xdir','reverse')
%     ylim([0 120])
%     xlim([0.4 10.6]);
%     title('Incertitude sur les 2P')
%     xlabel('Niveau 2p_x','FontSize',11,'fontweight','bold');
% 
%  end

    %% ==================================================
    %% ============== Intensit� des raies  ==============
    %% ==================================================
        if gaz==2 % Raie de n�on
            I585=densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2);             %2p1-1s2
            I640=densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5);             %2p9-1s5

            I_theo(t,:)=[I585 I640];
            
            %calcul d'erreur
            I585 =sqrt((sig_densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2))^2 + (densite2p(t,1)*Aij2p(1,2)*sig_Thetaij(1,2))^2);             %2p1-1s2
            I640 =sqrt((sig_densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5))^2 + (densite2p(t,9)*Aij2p(9,5)*sig_Thetaij(9,5))^2);             %2p9-1s5
            
            sig_I_theo(t,:) = [I585 I640];
                  
            %% Ainsi qu'� quel point elles sont effectivement sorties du plasma (sans �tre auto-absorb�es)
            Theta(t,:)=[Thetaij(1,2) Thetaij(9,5)] ;

            clear  I585 I640
        elseif gaz ==3 % Raies d'argon
             I667=densite2p(t,1)*Aij2p(1,4)*Thetaij(1,4);             %2p1-1s4
             I696=densite2p(t,2)*Aij2p(2,5)*Thetaij(2,5);             %2p2-1s5
             I706=densite2p(t,3)*Aij2p(3,5)*Thetaij(3,5);             %2p3-1s5
             I714=densite2p(t,4)*Aij2p(4,5)*Thetaij(4,5);             %2p3-1s5
             I727=densite2p(t,2)*Aij2p(2,4)*Thetaij(2,4);             %2p2-1s4
             I738=densite2p(t,3)*Aij2p(3,4)*Thetaij(3,4);             %2p3-1s4
             I750=densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2);             %2p1-1s2
             I751=densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4);             %2p5-1s4
             I763=densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5);             %2p6-1s5
             I794=densite2p(t,4)*Aij2p(4,3)*Thetaij(4,3);             %2p4-1s3
             I800=densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4);             %2p6-1s4
             I801=densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5);             %2p8-1s5 
             I810=densite2p(t,7)*Aij2p(7,4)*Thetaij(7,4);             %2p7-1s4
             I811=densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5);             %2p9-1s5
             I826=densite2p(t,2)*Aij2p(2,2)*Thetaij(2,2);             %2p2-1s2
             I840=densite2p(t,3)*Aij2p(3,2)*Thetaij(3,2);             %2p3-1s2
             I842=densite2p(t,8)*Aij2p(8,4)*Thetaij(8,4);             %2p8-1s4
             I852=densite2p(t,4)*Aij2p(4,2)*Thetaij(4,2);             %2p4-1s2
             I866=densite2p(t,7)*Aij2p(7,3)*Thetaij(7,3);             %2p7-1s3

             I_theo(t,:)=[I667 I696 I706 I714 I727 I738 I750 I751 I763 I794 I800 I801 I810 I811 I826 I840 I842 I852 I866]; %pump de Donnelly1.1;
             
             %calcul d'erreur
             I667=sqrt( (sig_densite2p(t,1)*Aij2p(1,4)*Thetaij(1,4))^2 + (densite2p(t,1)*Aij2p(1,4)*sig_Thetaij(1,4))^2);             %2p1-1s4
             I696=sqrt( (sig_densite2p(t,2)*Aij2p(2,5)*Thetaij(2,5))^2 + (densite2p(t,2)*Aij2p(2,5)*sig_Thetaij(2,5))^2);             %2p2-1s5
             I706=sqrt( (sig_densite2p(t,3)*Aij2p(3,5)*Thetaij(3,5))^2 + (densite2p(t,3)*Aij2p(3,5)*sig_Thetaij(3,5))^2);             %2p3-1s5
             I714=sqrt( (sig_densite2p(t,4)*Aij2p(4,5)*Thetaij(4,5))^2 + (densite2p(t,4)*Aij2p(4,5)*sig_Thetaij(4,5))^2);             %2p3-1s5
             I727=sqrt( (sig_densite2p(t,2)*Aij2p(2,4)*Thetaij(2,4))^2 + (densite2p(t,2)*Aij2p(2,4)*sig_Thetaij(2,4))^2);             %2p2-1s4
             I738=sqrt( (sig_densite2p(t,3)*Aij2p(3,4)*Thetaij(3,4))^2 + (densite2p(t,3)*Aij2p(3,4)*sig_Thetaij(3,4))^2);             %2p3-1s4
             I750=sqrt( (sig_densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2))^2 + (densite2p(t,1)*Aij2p(1,2)*sig_Thetaij(1,2))^2);             %2p1-1s2
             I751=sqrt( (sig_densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4))^2 + (densite2p(t,5)*Aij2p(5,4)*sig_Thetaij(5,4))^2);             %2p5-1s4
             I763=sqrt( (sig_densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5))^2 + (densite2p(t,6)*Aij2p(6,5)*sig_Thetaij(6,5))^2);             %2p6-1s5
             I794=sqrt( (sig_densite2p(t,4)*Aij2p(4,3)*Thetaij(4,3))^2 + (densite2p(t,4)*Aij2p(4,3)*sig_Thetaij(4,3))^2);             %2p4-1s3
             I800=sqrt( (sig_densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4))^2 + (densite2p(t,6)*Aij2p(6,4)*sig_Thetaij(6,4))^2);             %2p6-1s4
             I801=sqrt( (sig_densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5))^2 + (densite2p(t,8)*Aij2p(8,5)*sig_Thetaij(8,5))^2);             %2p8-1s5 
             I810=sqrt( (sig_densite2p(t,7)*Aij2p(7,4)*Thetaij(7,4))^2 + (densite2p(t,7)*Aij2p(7,4)*sig_Thetaij(7,4))^2);             %2p7-1s4
             I811=sqrt( (sig_densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5))^2 + (densite2p(t,9)*Aij2p(9,5)*sig_Thetaij(9,5))^2);             %2p9-1s5
             I826=sqrt( (sig_densite2p(t,2)*Aij2p(2,2)*Thetaij(2,2))^2 + (densite2p(t,2)*Aij2p(2,2)*sig_Thetaij(2,2))^2);             %2p2-1s2
             I840=sqrt( (sig_densite2p(t,3)*Aij2p(3,2)*Thetaij(3,2))^2 + (densite2p(t,3)*Aij2p(3,2)*sig_Thetaij(3,2))^2);             %2p3-1s2
             I842=sqrt( (sig_densite2p(t,8)*Aij2p(8,4)*Thetaij(8,4))^2 + (densite2p(t,8)*Aij2p(8,4)*sig_Thetaij(8,4))^2);             %2p8-1s4
             I852=sqrt( (sig_densite2p(t,4)*Aij2p(4,2)*Thetaij(4,2))^2 + (densite2p(t,4)*Aij2p(4,2)*sig_Thetaij(4,2))^2);             %2p4-1s2
             I866=sqrt( (sig_densite2p(t,7)*Aij2p(7,3)*Thetaij(7,3))^2 + (densite2p(t,7)*Aij2p(7,3)*sig_Thetaij(7,3))^2);             %2p7-1s3

             sig_I_theo(t,:)=[I667 I696 I706 I714 I727 I738 I750 I751 I763 I794 I800 I801 I810 I811 I826 I840 I842 I852 I866]; 

            %% Ainsi qu'� quel point elles sont effectivement sorties du plasma (sans �tre auto-absorb�es)
            Theta(t,:)=[Thetaij(1,4) Thetaij(2,5) Thetaij(3,5) Thetaij(4,5) Thetaij(2,4) Thetaij(3,4) Thetaij(1,2)...
                        Thetaij(5,4) Thetaij(6,5) Thetaij(4,3) Thetaij(6,4) Thetaij(8,5) Thetaij(7,4) Thetaij(9,5)...
                        Thetaij(2,2) Thetaij(3,2) Thetaij(8,4) Thetaij(4,2) Thetaij(7,3)];

            clear I667 I696 I706 I714 I727 I738 I750 I751 I763 I794 I800 I801 I810 I811 I826 I840 I842 I852 I866

        elseif gaz==4
             I758=densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4);             %2p5-1s4
             I760=densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5);             %2p6-1s5
             I768=densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2);             %2p1-1s2
             I769=densite2p(t,7)*Aij2p(7,5)*Thetaij(7,5);             %2p7-1s5
             I785=densite2p(t,3)*Aij2p(3,3)*Thetaij(3,3);             %2p3-1s3
             I810=densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5);             %2p8-1s5
             I811=densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5);             %2p9-1s5
             I819=densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4);             %2p6-1s4
             I829=densite2p(t,7)*Aij2p(7,4)*Thetaij(7,4);             %2p7-1s4 
             I850=densite2p(t,4)*Aij2p(4,2)*Thetaij(4,2);             %2p4-1s2
             I877=densite2p(t,8)*Aij2p(8,4)*Thetaij(8,4);             %2p8-1s4
             I892=densite2p(t,10)*Aij2p(10,5)*Thetaij(10,5);          %2p8-1s4

             I_theo(t,:)=[I758 I760 I768 I769 I785 I810 I811 I819 I829 I850 I877 I892];
             
             %calcul erreur   
             I758=sqrt( (sig_densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4))^2 + (densite2p(t,5)*Aij2p(5,4)*sig_Thetaij(5,4))^2);             %2p5-1s4
             I760=sqrt( (sig_densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5))^2 + (densite2p(t,6)*Aij2p(6,5)*sig_Thetaij(6,5))^2);             %2p6-1s5
             I768=sqrt( (sig_densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2))^2 + (densite2p(t,1)*Aij2p(1,2)*sig_Thetaij(1,2))^2);             %2p1-1s2
             I769=sqrt( (sig_densite2p(t,7)*Aij2p(7,5)*Thetaij(7,5))^2 + (densite2p(t,7)*Aij2p(7,5)*sig_Thetaij(7,5))^2);             %2p7-1s5
             I785=sqrt( (sig_densite2p(t,3)*Aij2p(3,3)*Thetaij(3,3))^2 + (densite2p(t,3)*Aij2p(3,3)*sig_Thetaij(3,3))^2);             %2p3-1s3
             I810=sqrt( (sig_densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5))^2 + (densite2p(t,8)*Aij2p(8,5)*sig_Thetaij(8,5))^2);             %2p8-1s5
             I811=sqrt( (sig_densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5))^2 + (densite2p(t,9)*Aij2p(9,5)*sig_Thetaij(9,5))^2);             %2p9-1s5
             I819=sqrt( (sig_densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4))^2 + (densite2p(t,6)*Aij2p(6,4)*sig_Thetaij(6,4))^2);             %2p6-1s4
             I829=sqrt( (sig_densite2p(t,7)*Aij2p(7,4)*Thetaij(7,4))^2 + (densite2p(t,7)*Aij2p(7,4)*sig_Thetaij(7,4))^2);             %2p7-1s4 
             I850=sqrt( (sig_densite2p(t,4)*Aij2p(4,2)*Thetaij(4,2))^2 + (densite2p(t,4)*Aij2p(4,2)*sig_Thetaij(4,2))^2);             %2p4-1s2
             I877=sqrt( (sig_densite2p(t,8)*Aij2p(8,4)*Thetaij(8,4))^2 + (densite2p(t,8)*Aij2p(8,4)*sig_Thetaij(8,4))^2);             %2p8-1s4
             I892=sqrt( (sig_densite2p(t,10)*Aij2p(10,5)*Thetaij(10,5))^2 + (densite2p(t,10)*Aij2p(10,5)*sig_Thetaij(10,5))^2);          %2p8-1s4

             sig_I_theo(t,:)=[I758 I760 I768 I769 I785 I810 I811 I819 I829 I850 I877 I892];%pump de Donnelly 1.48

            %% Ainsi qu'� quel point elles sont effectivement sorties du plasma (sans �tre auto-absorb�es)
            Theta(t,:)=[Thetaij(5,4) Thetaij(6,5) Thetaij(1,2) Thetaij(7,5) Thetaij(3,3) Thetaij(8,5)...
                        Thetaij(9,5) Thetaij(6,4) Thetaij(7,4) Thetaij(4,2) Thetaij(8,4) Thetaij(10,5)];

            clear I758 I760 I768 I769 I785 I810 I811 I819 I829 I850 I877 I892
        elseif gaz==5
              
                 I788=densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2);             %2p1-1s2 
                 I820=densite2p(t,4)*Aij2p(4,3)*Thetaij(4,3);             %2p4-1s3 
                 I823=densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5);             %2p6-1s5
                 I828=densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4);             %2p5-1s4 
                 I834=densite2p(t,3)*Aij2p(3,2)*Thetaij(3,2);             %2p3-1s2 
                 I881=densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5);             %2p5-1s4 
                 I895=densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4);             %2p6-1s5 
                 I904=densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5);             %2p6-1s5 
                
                  I_theo(t,:)=[I788 I820 I823 I828 I834 I881 I895 I904]; 
                  
                  %erreur avec une erreur sur les Aij puisque la seul
                  %source avec des Aij pour le Xenon en donne 2 diff�rents
                 I788=sqrt((sig_densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2))^2 + (densite2p(t,1)*Aij2p(1,2)*Thetaij(1,2)*0.18)^2 + (densite2p(t,1)*Aij2p(1,2)*sig_Thetaij(1,2))^2); % 18% sur Aij     
                 I820=sqrt((sig_densite2p(t,4)*Aij2p(4,3)*Thetaij(4,3))^2 + (densite2p(t,4)*Aij2p(4,3)*Thetaij(4,3)*0.32)^2 + (densite2p(t,4)*Aij2p(4,3)*sig_Thetaij(4,3))^2); % 32% sur Aij                
                 I823=sqrt((sig_densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5))^2 + (densite2p(t,6)*Aij2p(6,5)*Thetaij(6,5)*0.10)^2 + (densite2p(t,6)*Aij2p(6,5)*sig_Thetaij(6,5))^2); % 10% sur Aij               
                 I828=sqrt((sig_densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4))^2 + (densite2p(t,5)*Aij2p(5,4)*Thetaij(5,4)*0.32)^2 + (densite2p(t,5)*Aij2p(5,4)*sig_Thetaij(5,4))^2); % 32% sur Aij              
                 I834=sqrt((sig_densite2p(t,3)*Aij2p(3,2)*Thetaij(3,2))^2 + (densite2p(t,3)*Aij2p(3,2)*Thetaij(3,2)*0.27)^2 + (densite2p(t,3)*Aij2p(3,2)*sig_Thetaij(3,2))^2); % 27% sur Aij                 
                 I881=sqrt((sig_densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5))^2 + (densite2p(t,8)*Aij2p(8,5)*Thetaij(8,5)*0.24)^2 + (densite2p(t,8)*Aij2p(8,5)*sig_Thetaij(8,5))^2); % 24% sur Aij              
                 I895=sqrt((sig_densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4))^2 + (densite2p(t,6)*Aij2p(6,4)*Thetaij(6,4)*0.11)^2 + (densite2p(t,6)*Aij2p(6,4)*sig_Thetaij(6,4))^2); % 11% sur Aij               
                 I904=sqrt((sig_densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5))^2 + (densite2p(t,9)*Aij2p(9,5)*Thetaij(9,5)*0.09)^2 + (densite2p(t,9)*Aij2p(9,5)*sig_Thetaij(9,5))^2); %  9% sur Aij                 
              

            sig_I_theo(t,:)=[I788 I820 I823 I828 I834 I881 I895 I904];

            %% Ainsi qu'� quel point elles sont effectivement sorties du plasma (sans �tre auto-absorb�es)
            Theta(t,:)=[Thetaij(1,2) Thetaij(4,3) Thetaij(6,5) Thetaij(5,4) Thetaij(3,2) Thetaij(9,5) Thetaij(6,4) Thetaij(8,5)];

            clear I788 I820 I823 I828 I834 I881 I904  
        end

end %Fin Boucle Te
clear Thetaij depopulation Dopp VDWaals Res

%% Regroupement de l'information pour faciliter l'exportation et la lisibilit�
    Gains2p=zeros(5,length(Te),10);
    Pertes2p=zeros(3,length(Te),10);
    Gains1s=zeros(5,length(Te),5);
    Pertes1s=zeros(6,length(Te),5);
    Emission=zeros(6,length(Te),index);
    



    for i=1:length(Te)
        for j=1:10
            Gains2p(1,i,j)=GainFond(i,j);
            Gains2p(2,i,j)=GainNm(i,j);
            Gains2p(3,i,j)=GainColl2p(i,j);
            Gains2p(4,i,j)=GainRadTrap(i,j);
            Gains2p(5,i,j)=GainAr(i,j);
            Pertes2p(1,i,j)=PerteRadiatif(i,j);
            Pertes2p(2,i,j)=PerteColl2p(i,j);
            Pertes2p(3,i,j)=PerteColl1s(i,j);
        end
        for j=1:5
            Gains1s(1,i,j)=GainFond_1s(i,j);
            Gains1s(2,i,j)=GainAutoAbs_1s(i,j);
            Gains1s(3,i,j)=GainMixing_1s(i,j);
            Gains1s(4,i,j)=GainUpdown_1s(i,j);
            Gains1s(5,i,j)=GainPopAr_1s(i,j);
            Pertes1s(1,i,j)=PerteRad_1s(i,j);
            Pertes1s(2,i,j)=PerteVers2p_1s(i,j);
            Pertes1s(3,i,j)=PerteMixinxg_1s(i,j);
            Pertes1s(4,i,j)=PerteSuper_1s(i,j);
            Pertes1s(5,i,j)=PerteIonisation_1s(i,j);
            Pertes1s(6,i,j)=PerteNeutre_1s(i,j);      
        end
        %[lenght1 length2]=size(I_theo);
        for k=1:index 
            Emission(1,i,k)=I_theo(i,k);
            Emission(2,i,k)=Theta(i,k);
            Emission(3,i,k)=Doppler(i,k);
            Emission(4,i,k)=VanDerWaals(i,k);
            Emission(5,i,k)=Resonant(i,k);
            Emission(6,i,k)=sig_I_theo(i,k); %Nouveau! l'incertitude
        end
    end

end
