function [densite1s,densite2p,Gains2p,Pertes2p,Emission,ContributionFond2p,Mecanismes,energie1s,energie2p,Gains1s,Pertes1s] = TeIntTheo_TRG(gaz,ne,Te,rateGround_1s,rateGround_2p,rate1s_2p,rateQuenching,rateNeutral,sig_rateGround_1s,sig_rateGround_2p,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral,Tg,longueur,AllIntegral,P,IndexOffset,max,flow)

%% =====================================================================================================
%% ================ Cette fonction calcule les intensité théorique de raies pour tous ==================
%% =============== les Te à une densité des niveaux 1s et une température de gaz donnés ================
%% =====================================================================================================

%% ====== INFOS SUR LA FONCTION ======
%La variable gaz indique pour quel gaz on cherche les intensités de raie
%La variable ne est la densité électronique
%La variable Te contient toutes les températures qui doivent être analysées au cours de cette fonction
%Les variables rateGround_2p et rate1s_2p sont les taux de réaction des
    %collisions e-fondamental et e-1s
%La variable Tg contient la température du gaz    
%La variable longueur est la longueur du parcours optique dans le plasma
    %pour le calcul d'auto-absorption
%La variable AllIntegral contient l'information sur le résultat de
    %l'intégrale du profil Voigt pour différents beta possibles
%La variable P contient la pression d'opération
%La variable IndexOffset sert à sélectionner le bon taux de réaction
    %lorsque le Te optimal a été trouvé
%------------------------------------------------------------------------------------
%% ===== nomenclature =====
%les vecteur ont comme première entrée le numéro du gaz utilisé
% ex rate1s_2p(4,3,10.i) krypton, 1s3->2P10 pour Te(i)
% 1 = He (inexistant)
% 2 = Néon
% 3 = Argon
% 4 = Krypton
% 5 = Xénon
% 13 = O2
%% ====== DÉBUT DE LA FONCTION =======

%% Récupération des données utiles
[En2px_1s,En2px_2py,poids2p,poids1s,LambdaTheo2p,LambdaTheo1s,Aij2p,Aij1s,energie1s,energie2p,fji2p,fji1s]=TeInitialisation_TRG(gaz);
    %POUR INFO:
    %rateGround_2p = taux de réaction pour l'excitation des niveaux 2p par impact électronique sur le fondamental
    %rateGround_1s = Taux de réaction pour l'excitation du fondamental vers les 1s par impact électronique
    %rate1s_2p = taux de réaction pour l'excitation des niveaux 2p par impact électronique sur les niveaux 1s
    %rateQuenching = taux de réaction pour l'échange entre les 1s par impact électronique avec l'ionisation et le superélastique en 6 et 7iem position
    %rateNeutral = Quenching des métastables vers le fondamental par les neutres ex : (3,4) = neutre Argon sur métastables de Krypton
    %En2px_1s = taux de transfert d'énergie par collision avec les neutre des niveaux 2p (menant à 1s)
    %En2px_2py = taux de transfert d'énergie par collision des niveaux 2p
    %poids2p = poids statistique des niveaux 2p
    %poids1s = poids statistique des niveaux 1s
    %LambdaTheo = longueurs d'onde des transition radiatives 2p-1s
    %Aij = Aij, les coefficients de desexcitation radiative spontanée
    %energie1s = énergie des niveaux 1s
    
    %% Préallocation d'espace
    %selon le nombre de raie observé par gaz
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
    PerteVers2p_1s     =zeros(length(Te),5);
    PerteRad_1s        =zeros(length(Te),5);
    PerteMixinxg_1s    =zeros(length(Te),5);
    PerteSuper_1s      =zeros(length(Te),5);
    PerteIonisation_1s =zeros(length(Te),5);
    PerteNeutre_1s     =zeros(length(Te),5);
    %tout ce qui est incertitudes
    sig_gain           =zeros(length(Te),5);
    sig_perte          =zeros(length(Te),5);
    gain               =zeros(length(Te),5);
    perte              =zeros(length(Te),5);
    sig_densite1s      =zeros(length(Te),5);
    sig_gain_2P        =zeros(length(Te),10);
    sig_gain_2P_a      =zeros(10,1);
    gain_2P            =zeros(length(Te),10);
    sig_densite2p      =zeros(length(Te),10);

    
    
    temp_I_theo     =zeros(length(Te),index);
    I_theo          =zeros(length(Te),index);
    sig_I_theo      =zeros(length(Te),index);
    Theta           =zeros(length(Te),index);
    Doppler         =zeros(length(Te),index);
    VanDerWaals     =zeros(length(Te),index);
    Resonant        =zeros(length(Te),index);
    
    Thetaij         =zeros(10,5);
    Thetaij_1s      =zeros(1,5);
    Dopp            =zeros(10,5);
    VDWaals         =zeros(10,5);
    Res             =zeros(10,5);
    
    Mecanismes   =2*ones(length(Te),1);
    
    
    %% Correction des taux collisionnels en température
    En2px_1s=En2px_1s*(sqrt(Tg/300));
    En2px_2py=En2px_2py*(sqrt(Tg/300));    

%% Calcul de la densité de neutres à l'état fondamental avec p=n*k*Tg [m-3]
   %ng(gaz) car P est un vecteur avec les pressions partielles des gas utilisés
    ng=(P)./(1.38064852*10^(-23)*Tg);
    

   
%% ================== Boucle sur la température électronique: COMPTEUR: t =================
for t=1:length(Te)

%% ##################################################################################################  
%  ##################################################################################################
%  ############# D'abord on calcule la densité des niveaux 1s via un bilan gains-pertes #############
%  ##################################################################################################
%  ##################################################################################################
%% Calcul de l'élargissement et du piégeage optique pour les transitions radiatives 1s-fondamental
    for i=1:5
        if Aij1s(i)~=0
            Thetaij_1s(i) = 1; % on met l'autoabsorption à 0 pour tester Donnelly (autoabs = 1-thetaij)
            %[Thetaij_1s(i)]=TeEscapeFactorDOPE_TRG(gaz,LambdaTheo1s(i),poids1s(i),0,energie1s(i),1,0,0,Aij1s(i),ng(gaz),longueur,Tg,ng(gaz),fji1s(i),AllIntegral);        
        end
    end 
%% Mécanismes indépendants de la densité des 1s
[PopFond]= TePopulation_Metastable_TRG(ng(gaz),ne,rateGround_1s(:,t+IndexOffset-1));
population1s=PopFond;
    
%% Mécanismes dépendant de la densité des 1s    
[Depop2p,DepopRadFond,PopAutoAbs,DepopMix,PopMix,DepopSuperelestique,DepopIonisation,DepopNeutre,PopUpDown] = TeDepopulation_Metastable_TRG(t+IndexOffset-1,ng,ne,Aij1s,Thetaij_1s,rate1s_2p,rateQuenching,rateNeutral,Aij2p);

depopulation1s=Depop2p+DepopRadFond-PopAutoAbs+DepopMix-PopMix+DepopSuperelestique+DepopIonisation+DepopNeutre-PopUpDown;



%% Obtention de la densité des niveaux pour avoir l'etat stationnaire depopulation=population
% n1sX=zeros(1,5);
% for iter=1:10
%     gains=zeros(1,5);
%     pertes=zeros(1,5);
%     for j=2:1:5
%         for i=2:5 % comme sur les 1s pour les n1sX(j) à résoudre
%             if j~=i 
%               gains(j) = gains(j) + n1sX(i)*PopMix(j,i) + n1sX(i)*PopUpDown(j,i);
%             end
%         end
%         pertes(j) = DepopSuperelestique(j,j) + DepopIonisation(j,j) + DepopNeutre(j,j) + Depop2p(j,j) + DepopMix(j,j);
%         gains(j)=gains(j)+PopFond(j);     
%         n1sX(j)=gains(j)/pertes(j);
%     end
% end

n1sX=linsolve(depopulation1s,population1s);

 n1sX(2)=0; %negliger les résonnant
 n1sX(4)=0; %negliger les résonnant

densite1s(t,:)=n1sX;

densite1s(isnan(densite1s)) = 0 ; %enleve le NaN en 1s1
%% ============== Suivi de l'influence des processus de peuplement et de dépeuplement des niveaux 1s ==============   

 GainFond_1s(t,:)      =  100*PopFond                       ./(PopFond + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainAutoAbs_1s(t,:)   =  100*PopAutoAbs*densite1s(t,:)'    ./(PopFond + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainMixing_1s(t,:)    =  100*PopMix*densite1s(t,:)'        ./(PopFond + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 GainUpdown_1s(t,:)    =  100*PopUpDown*densite1s(t,:)'     ./(PopFond + PopAutoAbs*densite1s(t,:)' + PopMix*densite1s(t,:)' + PopUpDown*densite1s(t,:)');
 PerteVers2p_1s(t,:)     =  100*Depop2p*densite1s(t,:)'             ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteRad_1s(t,:)        =  100*DepopRadFond*densite1s(t,:)'        ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteMixinxg_1s(t,:)    =  100*DepopMix*densite1s(t,:)'            ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteSuper_1s(t,:)      =  100*DepopSuperelestique*densite1s(t,:)' ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteIonisation_1s(t,:) =  100*DepopIonisation*densite1s(t,:)'     ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');
 PerteNeutre_1s(t,:)     =  100*DepopNeutre*densite1s(t,:)'         ./(Depop2p*densite1s(t,:)' + DepopRadFond*densite1s(t,:)' + DepopMix*densite1s(t,:)' + DepopSuperelestique*densite1s(t,:)' + DepopIonisation*densite1s(t,:)' + DepopNeutre*densite1s(t,:)');

%% =========================================================== 
%  =========================================================== 
%  ================ Calcul de l'incertitude ==================
%  =======================sur les 1s ========================= 
%  ===========================================================  
%  =========================================================== 
%calcul des incertide sur les taux de réactions
[sig_PopFond]= TePopulation_Metastable_TRG(ng(gaz),ne,sig_rateGround_1s(:,t+IndexOffset-1));
[sig_Depop2p,sig_DepopRadFond,sig_PopAutoAbs,sig_DepopMix,sig_PopMix,sig_DepopSuperelestique,sig_DepopIonisation,sig_DepopNeutre,sig_PopUpDown] = TeDepopulation_Metastable_TRG(t+IndexOffset-1,ng,ne,Aij1s,Thetaij_1s,sig_rate1s_2p,sig_rateQuenching,sig_rateNeutral,Aij2p);

iter=0; % boucle qui compte le nombre d'itération
residu=10e21; %start le résidu à quelquechose de très gros.
% boucle tant que la différence entre 2 résultats consécutif est >1000 fois
% celui-ci
%while iter<3
while abs(residu) > sig_densite1s(t,5)/1000 %boucle car 1s3 dépend de 1s5 et vice versa
for j=2:5 % boucle sur les 1s à résoudre
    for i=2:5 % comme sur les 1s pour les n1sX(j) à résoudre
        if j~=i 
          gain(t,j)     = gain(t,j)     + densite1s(t,i)*PopMix(j,i) + densite1s(t,i)*PopUpDown(i,j);      
          sig_gain(t,j) = sig_gain(t,j) + densite1s(t,i).^2*(sum(sig_PopMix(j,i) + sig_PopUpDown(j,i))).^2 + sig_densite1s(i).^2*(sum(PopMix(j,i) + PopUpDown(j,i))).^2; % sqrt à la fin
        end
    end
    perte(t,j)        =    Depop2p(j,j)+DepopRadFond(j,j)+DepopMix(j,j)+DepopSuperelestique(j,j)+DepopIonisation(j,j)+DepopNeutre(j,j);
    sig_perte(t,j)    =    sqrt(sig_Depop2p(j,j)^2+sig_DepopRadFond(j,j)^2+sig_DepopMix(j,j)^2+sig_DepopSuperelestique(j,j)^2+sig_DepopIonisation(j,j)^2+sig_DepopNeutre(j,j)^2);
    gain(t,j)         =    gain(t,j)+PopFond(j);
    sig_gain(t,j)     =    sqrt(sig_gain(t,j)+sig_PopFond(j).^2);
end    
clear i j
iter=iter+1;
%calcul du résidu entre 2 itérations succesives
residu = sig_densite1s(t,5) - densite1s(t,5).*sqrt((sig_gain(t,5)./gain(t,5)).^2 + (sig_perte(t,5)./perte(t,5)).^2);
%incertitude sur la Densité de métsstables
sig_densite1s(t,:) =  densite1s(t,:).*sqrt((sig_gain(t,:)./gain(t,:)).^2 + (sig_perte(t,:)./perte(t,:)).^2);
end
sig_densite1s(isnan(sig_densite1s)) = 0 ; %enleve le NaN en 1s1
residu(isnan(residu)) = 0; %enleve le NaN en 1s1

 %% affichage des incertitudes sur les 1s
%  if IndexOffset>2
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
%  end


%iter %montre le nombre d'itération faites
%% ###########################################################  
%  ###########################################################
%  ############# Ensuite on peut résoudre les 2p #############
%  ###########################################################
%  ########################################################### 
%% Calcul de l'élargissement et du piégeage optique pour les transitions radiatives 2p-1s
    for i=1:10  %Boucle sur les 2p
        for j=1:5   %Boucle sur les 1s
            if Aij2p(i,j)~=0
                Thetaij(i,j) = 1; %on néglige l'autoabsorption
               % [Thetaij(i,j),Dopp(i,j)]=TeEscapeFactorDOPE_TRG(gaz,LambdaTheo2p(i,j),poids2p(i),1,energie2p(i),poids1s(j),0,energie1s(j),Aij2p(i,j),n1sX(j),longueur,Tg,ng(gaz),fji2p(i,j),AllIntegral); 
            end
        end      
    end

%% Calcul des processus de population/depopulation dépendant de la densité des niveaux 2p (inclut: émission et absorption radiative, processus collisionnel entre les 2p et vers les 1s)
    [DepopRad2p,PopRad2p,DepopColl2p,DepopColl1s,PopColl2p]=TeDepopulation_TRG(Aij2p,Thetaij,En2px_1s,En2px_2py,ng(gaz));   %Considère desexc. rad. + transfert coll vers le haut et le bas
  
%% Calcul des processus de population des niveaux 2p via impact électronique sur le fondamental et les niveaux 1s
    [PopFond2p,PopNm2p,Pop1s]=TePopulation_TRG(t+IndexOffset-1,ng(gaz),ne,n1sX,rate1s_2p,rateGround_2p,gaz);

%% Réorganisation des matrices: Les fonction flipud et rot90 sont utilisées car à priori les matrices contiennent en première position 2p10, 
%en 2e position 2p9, etc. Elles permettent donc de mettre ce qui correspond à 2p1 en première position, 2p2 en seconde position, etc
    DepopRad2p=rot90(DepopRad2p,2); %Le ,2 est parce qu'on fait une rotation de 2x 90deg
    PopRad2p=rot90(PopRad2p,2);
    DepopColl2p=rot90(DepopColl2p,2);
    DepopColl1s=rot90(DepopColl1s,2);
    PopColl2p=rot90(PopColl2p,2);
    PopFond2p=flipud(PopFond2p);
    PopNm2p=flipud(PopNm2p);
    
    %% Obtention de la densité des niveaux pour avoir l'etat stationnaire depopulation=population
    %depopulation=DepopRad2p-PopRad2p+DepopColl2p+DepopColl1s-PopColl2p;
    depopulation=DepopRad2p; % on prend juste les trans Radiative
    population=PopNm2p+PopFond2p;
    
    densite2p(t,:)=linsolve(depopulation,population);              %depopulation*density=population
    
    
    %% ============== Suivi de l'influence des processus de peuplement et de dépeuplement des niveaux 2p ==============   
    GainFond(t,:)=   100*PopFond2p                 ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p);
    GainNm(t,:)=     100*PopNm2p                   ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p);
    GainColl2p(t,:)= 100*PopColl2p*densite2p(t,:)'  ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p);    
    GainRadTrap(t,:)=100*PopRad2p*densite2p(t,:)'     ./(PopColl2p*densite2p(t,:)'+PopRad2p*densite2p(t,:)'+PopFond2p+PopNm2p);
    
    PerteRadiatif(t,:)= (100*DepopRad2p*densite2p(t,:)')./((DepopRad2p+DepopColl2p+DepopColl1s)*densite2p(t,:)');
    PerteColl2p(t,:)=(100*DepopColl2p*densite2p(t,:)')./((DepopRad2p+DepopColl2p+DepopColl1s)*densite2p(t,:)');
    PerteColl1s(t,:)=(100*DepopColl1s*densite2p(t,:)')./((DepopRad2p+DepopColl2p+DepopColl1s)*densite2p(t,:)');
    
    %% En particulier, on regarde pour quels paramètres l'excitation du fondamental (ou à l'inverse via les 1s) domine
    ContributionFond2p(t,:)=100*PopFond2p./(PopFond2p+PopNm2p);
    
    ground=0; met=0;    %Compteurs
    for m=1:10 %Pour tous les niveaux
       if ContributionFond2p(t,m)>90  %On regarde s'il y en a peuplés à plus de 90% par le fondamental
           ground=ground+1;
       end
       if ContributionFond2p(t,m)<10  %Et s'il y en a pleuplés à moins de 10% par le fondamental
           met=met+1;
       end
    end
    
    if ground==10   %Si tous les niveaux sont peuplés à plus de 90% du fondamental
        Mecanismes(t)=1;  %Alors le fondamental domine!
    end        
    %Si au contraire tous les niveaux sont peuplés à moins de 10% du le fondamental
    if met==10
        Mecanismes(t)=3;  %Alors les métastables dominent!
    end
    

    %% =========================================================== 
%  =========================================================== 
%  ================ Calcul de l'incertitude ==================
%  =======================sur les 2P ========================= 
%  ===========================================================  
%  =========================================================== 
    [sig_PopFond2p,sig_PopNm2p,sig_Pop1s]=TePopulation_TRG(t+IndexOffset-1,ng(gaz),ne,n1sX,sig_rate1s_2p,sig_rateGround_2p,gaz);
    %Aucune incertitude sur le Depop 
   
for j=1:10 % boucle sur les 2P à résoudre
    for i=2:5 % Boucle sur les 1s
      sig_gain_2P_a(j) = sig_gain_2P_a(j) + (sig_densite1s(t,i)*Pop1s(j,i))^2;
    end
    sig_gain_2P(t,11-j) = sqrt(sig_gain_2P_a(j) + sig_PopFond2p(j)^2 + sig_PopNm2p(j)^2);%11-j car gain(1)=2p10 et gain(10)=2p1
    gain_2P(t,j) = PopFond2p(j)+ PopNm2p(j);
end    
clear i j
%Valeur final de lincertitude sur les densité de 2P
sig_densite2p(t,:) = densite2p(t,:).*(sig_gain_2P(t,:)./gain_2P(t,:));
sig_densite2p(isnan(sig_densite2p)) = 0; %enleve les NaN


%% affichage des différentes contributions à l'erreur
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

    %% ============== Intensité des raies tenant compte de l'auto-absorption ==============
    for iter=1:1 %calcul I_theo et sig_I_the0
        if iter==2 %si on calcul sig_I_theo
            Xe_densite2p = densite2p; %besoin pour le Xe, terme croisé avec les Aij et la densité de 2P
            densite2p = sig_densite2p; %on met la densité de 2P comme Sigma
            
        end  
        if gaz==2 % Raie de néon
            I585=PopFond2p(1)*Aij2p(1,2)/sum(Aij2p(1,:));             %2p1-1s2
            I640=PopFond2p(9)*Aij2p(9,5)/sum(Aij2p(9,:));             %2p9-1s5

            temp_I_theo(t,:)=[I585 I640];

            %% Ainsi qu'à quel point elles sont effectivement sorties du plasma (sans être auto-absorbées)
            Theta(t,:)=[Thetaij(1,2) Thetaij(9,5)] ;

            %% Et leur diverses composantes à l'élargissement
            %à refaire!!!! pas fini avec les nouvelles raiees de TRG
            %Doppler(t)=Dopp(1,2);
            %VanDerWaals(t)=VDWaals(1,2);
            %Resonant(t)=Res(1,2);    

            clear  I585 I640
        elseif gaz ==3 % Raies d'argon
             I667=PopFond2p(1)*1      ; %2p1-1s4
             I696=PopFond2p(2)*1     ; %2p2-1s5
             I706=PopFond2p(3)*1    ;         %2p3-1s5
             I714=PopFond2p(4)*1   ;         %2p4-1s5
             I727=PopFond2p(2)*1  ;          %2p2-1s4
             I738=PopFond2p(3)*1 ;            %2p3-1s4
             I750=PopFond2p(1)*0.99      ;       %2p1-1s2
             I751=PopFond2p(5)*0.95     ;       %2p5-1s4
             I763=PopFond2p(6)*0.7     ;        %2p6-1s5
             I794=PopFond2p(4)*0.55   ;          %2p4-1s3
             I800=PopFond2p(6)*0.13 ;        %2p6-1s4 
             I801=PopFond2p(8)*0.28;%*1.1 ;         %2p8-1s5   %FUDGE 1.1
             I810=PopFond2p(7)*0.72;            %2p7-1s4
             I811=PopFond2p(9)*1         ;   %2p9-1s5
             I826=PopFond2p(2)*0.45     ;      %2p2-1s2
             I840=PopFond2p(3)*0.67    ;        %2p3-1s2
             I842=PopFond2p(8)*0.67   ;       %2p8-1s4
             I852=PopFond2p(4)*0.43  ;           %2p4-1s2
             I866=PopFond2p(7)*0.083;            %2p7-1s3

             temp_I_theo(t,:)=[I667 I696 I706 I714 I727 I738 I750 I751 I763 I794 I800 I801 I810 I811 I826 I840 I842 I852 I866]; %pump de Donnelly1.1;

            %% Ainsi qu'à quel point elles sont effectivement sorties du plasma (sans être auto-absorbées)
            Theta(t,:)=[Thetaij(1,4) Thetaij(2,5) Thetaij(3,5) Thetaij(4,5) Thetaij(2,4) Thetaij(3,4) Thetaij(1,2)...
                        Thetaij(5,4) Thetaij(6,5) Thetaij(4,3) Thetaij(6,4) Thetaij(8,5) Thetaij(7,4) Thetaij(9,5)...
                        Thetaij(2,2) Thetaij(3,2) Thetaij(8,4) Thetaij(4,2) Thetaij(7,3)];

            %% Et leur diverses composantes à l'élargissement
            %à refaire!!!! pas fini avec les nouvelles raiees de TRG
            %Doppler(t,:)=[Dopp(1,4) Dopp(2,5) Dopp(3,5) Dopp(4,5) Dopp(2,4) Dopp(3,4) Dopp(1,2) Dopp(5,4) Dopp(6,5) Dopp(4,3) Dopp(6,4) Dopp(8,5) Dopp(7,4) Dopp(9,5) Dopp(2,2) Dopp(3,2) Dopp(8,4) Dopp(4,2) Dopp(7,3)];
            %VanDerWaals(t,:)=[VDWaals(1,4) VDWaals(2,5) VDWaals(3,5) VDWaals(4,5) VDWaals(2,4) VDWaals(3,4) VDWaals(1,2) VDWaals(5,4) VDWaals(6,5) VDWaals(4,3) VDWaals(6,4) VDWaals(8,5) VDWaals(7,4) VDWaals(9,5) VDWaals(2,2) VDWaals(3,2) VDWaals(8,4) VDWaals(4,2) VDWaals(7,3)];
            %Resonant(t,:)=[Res(1,2) Res(2,5) Res(2,4) Res(2,2) Res(3,5) Res(3,4) Res(3,2) Res(4,3) Res(4,2) Res(5,4) Res(6,5) Res(6,4) Res(7,4) Res(7,3) Res(8,5) Res(8,4) Res(9,5)];    

            clear I667 I696 I706 I714 I727 I738 I750 I751 I763 I794 I800 I801 I810 I811 I826 I840 I842 I852 I866

        elseif gaz==4
             I758=PopFond2p(5)*1         ;    %2p5-1s4
             I760=PopFond2p(6)*0.76     ;       %2p6-1s5
             I768=PopFond2p(1)*1       ;   %2p1-1s2
             I769=PopFond2p(7)*0.14   ;          %2p7-1s5
             I785=PopFond2p(3)*0.55  ;         %2p3-1s3
             I810=PopFond2p(8)*0.34     ;        %2p8-1s5
             I811=PopFond2p(9)*1       ;     %2p9-1s5
             I819=PopFond2p(6)*0.25   ;         %2p6-1s4
             I829=PopFond2p(7)*0.86;%*(0.7/0.6);           %2p7-1s4 FUDGE Donnelly
             I850=PopFond2p(4)*0.54 ;            %2p4-1s2
             I877=PopFond2p(8)*0.66;            %2p8-1s4
             I892=PopFond2p(10)*Aij2p(10,5)/sum(Aij2p(10,:));          %2p8-1s4

             temp_I_theo(t,:)=[I758 I760 I768 I769 I785 I810 I811 I819 I829 I850 I877 I892];%pump de Donnelly 1.48

            %% Ainsi qu'à quel point elles sont effectivement sorties du plasma (sans être auto-absorbées)
            Theta(t,:)=[Thetaij(5,4) Thetaij(6,5) Thetaij(1,2) Thetaij(7,5) Thetaij(3,3) Thetaij(8,5)...
                        Thetaij(9,5) Thetaij(6,4) Thetaij(7,4) Thetaij(4,2) Thetaij(8,4) Thetaij(10,5)];
            %% Et leur diverses composantes à l'élargissement
            %à refaire!!!! pas fini avec les nouvelles raiees de TRG
            %Doppler(t,:)=[Dopp(5,4) Dopp(6,5) Dopp(1,2) Dopp(7,5) Dopp(3,3) Dopp(8,5) Dopp(9,5) Dopp(6,4) Dopp(7,4) Dopp(4,2) Dopp(8,4)];
            %VanDerWaals(t,:)=[VDWaals(5,4) VDWaals(6,5) VDWaals(1,2) VDWaals(7,5) VDWaals(3,3) VDWaals(8,5) VDWaals(9,5) VDWaals(6,4) VDWaals(7,4) VDWaals(4,2) VDWaals(8,4)];
            %Resonant(t,:)=[Res(1,4) Res(2,5) Res(3,5) Res(4,5) Res(2,4) Res(3,4) Res(1,2) Res(5,4) Res(6,5) Res(4,3) Res(6,4) Res(8,5) Res(7,4) Res(9,5) Res(2,2) Res(3,2) Res(8,4) Res(4,2) Res(7,3)];  

            clear I758 I760 I768 I769 I785 I810 I811 I819 I829 I850 I877 I892
        elseif gaz==5
                if iter==1 % si on fait I_theo, normal
                 I788=PopFond2p(1)*0.72       ;      %2p1-1s2 %corrigé par simon
                 I820=PopFond2p(4)*0.41      ;       %2p4-1s3 %corrigé par simon
                 I823=PopFond2p(6)*0.75     ;      %2p6-1s5 %corrigé par simon
                 I828=PopFond2p(5)*1       ;      %2p5-1s4 %corrigé par simon
                 I834=PopFond2p(3)*0.86   ;          %2p3-1s2 %corrigé par simon
                 I881=PopFond2p(8)*1     ;        %2p5-1s4 %corrigé par simon
                 I895=PopFond2p(6)*1    ;         %2p6-1s5 %corrigé par simon
                 I904=PopFond2p(9)*0.38;             %2p6-1s5 %corrigé par simon
                else %calcul d'erreur avec les Aij incertains du Xe
                 I788=sqrt((PopFond2p(1)*Aij2p(1,2)*Thetaij(1,2))^2 + (Xe_PopFond2p(1)*Aij2p(1,2)*Thetaij(1,2)*0.18)^2); % 18% sur Aij     
                 I820=sqrt((PopFond2p(4)*Aij2p(4,3)*Thetaij(4,3))^2 + (Xe_PopFond2p(4)*Aij2p(4,3)*Thetaij(4,3)*0.32)^2); % 32% sur Aij                
                 I823=sqrt((PopFond2p(6)*Aij2p(6,5)*Thetaij(6,5))^2 + (Xe_PopFond2p(6)*Aij2p(6,5)*Thetaij(6,5)*0.10)^2); % 10% sur Aij               
                 I828=sqrt((PopFond2p(5)*Aij2p(5,4)*Thetaij(5,4))^2 + (Xe_PopFond2p(5)*Aij2p(5,4)*Thetaij(5,4)*0.32)^2); % 32% sur Aij              
                 I834=sqrt((PopFond2p(3)*Aij2p(3,2)*Thetaij(3,2))^2 + (Xe_PopFond2p(3)*Aij2p(3,2)*Thetaij(3,2)*0.27)^2); % 27% sur Aij                 
                 I881=sqrt((PopFond2p(8)*Aij2p(8,5)*Thetaij(8,5))^2 + (Xe_PopFond2p(8)*Aij2p(8,5)*Thetaij(8,5)*0.24)^2); % 24% sur Aij              
                 I895=sqrt((PopFond2p(6)*Aij2p(6,4)*Thetaij(6,4))^2 + (Xe_PopFond2p(6)*Aij2p(6,4)*Thetaij(6,4)*0.11)^2); % 11% sur Aij               
                 I904=sqrt((PopFond2p(9)*Aij2p(9,5)*Thetaij(9,5))^2 + (Xe_PopFond2p(9)*Aij2p(9,5)*Thetaij(9,5)*0.09)^2); %  9% sur Aij                 
                end

             temp_I_theo(t,:)=[I788 I820 I823 I828 I834 I881 I895 I904]; 

            %% Ainsi qu'à quel point elles sont effectivement sorties du plasma (sans être auto-absorbées)
            Theta(t,:)=[Thetaij(1,2) Thetaij(4,3) Thetaij(6,5) Thetaij(5,4) Thetaij(3,2) Thetaij(9,5) Thetaij(6,4) Thetaij(8,5)];

            %% Et leur diverses composantes à l'élargissement
            %à refaire!!!! pas fini avec les nouvelles raiees de TRG
            %Doppler(t,:)=[Dopp(2,5) Dopp(3,5) Dopp(2,4) Dopp(3,4) Dopp(1,2) Dopp(5,4) Dopp(6,5)];
            %VanDerWaals(t,:)=[VDWaals(2,5) VDWaals(3,5) VDWaals(2,4) VDWaals(3,4) VDWaals(1,2) VDWaals(5,4) VDWaals(6,5)];
            %Resonant(t,:)=[Res(2,5) Res(3,5) Res(2,4) Res(3,4) Res(1,2) Res(5,4) Res(6,5)];    
            clear I788 I820 I823 I828 I834 I881 I904  
        end
            if iter==1 %premiere itération calcul I_theo.
             I_theo(t,:) = temp_I_theo(t,:);
            
             sig_I_theo(t,:) = 0.15*temp_I_theo(t,:);     
            end
        
    end
    %%écriture des intensité des UV pour Pierre
    if max==1
      if gaz==3
          x(:) = densite1s(t,:)
          fopen( 'Intensité_UV.txt','a')
              fileID = fopen( 'Intensité_UV.txt','at'); 
              formatSpec = '%.4e\t %d\t %.4e\t %d\r\n';
              fprintf(fileID,formatSpec,x(2:5));
              fclose(fileID);
          clear Ratio
      end
    end
    
end %Fin Boucle Te
clear Thetaij depopulation Dopp VDWaals Res

%% Regroupement de l'information pour faciliter l'exportation et la lisibilité
    Gains2p=zeros(4,length(Te),10);
    Pertes2p=zeros(3,length(Te),10);
    Gains1s=zeros(4,length(Te),5);
    Pertes1s=zeros(6,length(Te),5);
    Emission=zeros(6,length(Te),index);
    



    for i=1:length(Te)
        for j=1:10
            Gains2p(1,i,j)=GainFond(i,j);
            Gains2p(2,i,j)=GainNm(i,j);
            Gains2p(3,i,j)=GainColl2p(i,j);
            Gains2p(4,i,j)=GainRadTrap(i,j);
            Pertes2p(1,i,j)=PerteRadiatif(i,j);
            Pertes2p(2,i,j)=PerteColl2p(i,j);
            Pertes2p(3,i,j)=PerteColl1s(i,j);
        end
        for j=1:5
            Gains1s(1,i,j)=GainFond_1s(i,j);
            Gains1s(2,i,j)=GainAutoAbs_1s(i,j);
            Gains1s(3,i,j)=GainMixing_1s(i,j);
            Gains1s(4,i,j)=GainUpdown_1s(i,j);
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
