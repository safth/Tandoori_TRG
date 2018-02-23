%function [rate1s_2p_Ar, rateGround_2p_Ar]=TeTauxDeReactions(Te)
%% ============================================================================
%% ======= Cette fonction permet le calcul initial des taux de réaction =======
%% ======================= pour tous les Te considérés ========================
%% ============================================================================
 %% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%% TRES IMPORTANT DE CHANGER DANS L'ENDROIT QUI ASSOCIE LES TAUX DE
%%% RÉACTION AU PAS SÉLECTIONNÉ SI ONT CHANGE DE PAS DE O.01 À PLUS PETIT.
%%% DANS TeTauxDeReaction_TRG.m LIGNE 31
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Tous les taux de réactions sont ceux de TEOES 2004 pour l'instant. 
% plus ceux de gs->résonant, on a eux de Ajello,
% pour les résonant vers 2P, on a YORK
% pour le mixing des résonant, on a LXCAT NGFSRDW
% et du time reversal pour le reste qui est pas sur LXCAT.

 Te = 0.05:0.01:15;
 exposant=0.1:0.01:2;


%%calcul des constante de normalisation de K
e=1.602176462*10^(-19);
me=9.10938188*10^(-31);             %Mass of electron
kb=1.3806503*10^(-23);               %Boltzmann Constant

tic
%% ################################################################################################################
%##################################################################################################################
%##################################################################################################################
%############################################# Calcul pour l'Argon ################################################
%##################################################################################################################
%##################################################################################################################
for e=1:length(exposant)
       disp('Début du code')

        energie1s=[0 11.828 11.723 11.624 11.549]; %1s2 à 1s5
        poids=[0 3 1 3 5]; %0 1s2 ..1s5

    for i=1:length(Te)
        % =============== Rate excitation processes ============================
        %.txt file : First column : Energy // Second column : Cross section
        % here I use ratetrapz_TRG function because the corresponding cross section was already in m2 unit unlike others where we need to
        % multiply by 1e-23 in the function

        %=====Donnelly data for metastable excitation=====

        rate1s_2p(e,3,5,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,5,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-2p1.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,3,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-2p1.txt',   exposant(e), 0, 0);


        %=====Donnelly data for ground state excitation======
        rateGround_2p(e,3,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p10.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p9.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p8.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p7.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p6.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p5.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p4.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p3.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p2.txt',   exposant(e), 0, 0);
        rateGround_2p(e,3,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p1.txt',   exposant(e), 0, 0);

        rateGround_2p_HighP(e,3,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p10.txt',  exposant(e), 1, 10);
        rateGround_2p_HighP(e,3,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p9.txt',  exposant(e), 1, 9);
        rateGround_2p_HighP(e,3,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p8.txt',  exposant(e), 1, 8);
        rateGround_2p_HighP(e,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p7.txt',  exposant(e), 1, 7);
        rateGround_2p_HighP(e,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p6.txt',  exposant(e), 1, 6);
        rateGround_2p_HighP(e,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p5.txt',  exposant(e), 1, 5);
        rateGround_2p_HighP(e,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p4.txt',  exposant(e), 1, 4);
        rateGround_2p_HighP(e,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p3.txt',  exposant(e), 1, 3);
        rateGround_2p_HighP(e,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p2.txt',  exposant(e), 1, 2);
        rateGround_2p_HighP(e,3,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p1.txt',  exposant(e), 1, 1);
        
        rateGround_1s(e,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-1s3.txt',   exposant(e), 0, 0);
        rateGround_1s(e,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-1s5.txt',   exposant(e), 0, 0);


        %=============================================================
        %====== 6 pour le ground et 7 pour l'ionisation ======
        %=============================================================
         rateQuenching_1s(e,3,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-1s2.txt',   exposant(e), 0, 0);
         rateQuenching_1s(e,3,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-1s4.txt',   exposant(e), 0, 0);
         rateQuenching_1s(e,3,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-1s5.txt',   exposant(e), 0, 0);
         rateQuenching_1s(e,3,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-gs.txt',   exposant(e), 0, 0);
         rateQuenching_1s(e,3,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s3-ionization.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,3,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,3,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,3,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-1s3.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,3,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,3,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/1s5-ionization.txt',   exposant(e), 0, 0);

        %=============================================================
        % ===================Ajout des résonant=======================
        %=============================================================
            % plus ceux de Antoine de york
        rate1s_2p(e,3,2,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p1.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,2,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars2-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p1.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,3,4,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/ars4-2p10.txt',   exposant(e), 0, 0);

        rateGround_1s(e,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/Ajellogs-1s2.txt',   exposant(e), 0, 0);
        rateGround_1s(e,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/Ajellogs-1s4.txt',   exposant(e), 0, 0);

         %resonant  quelques sections de LXCAT  NGFSRDW et du time reversal

         rateQuenching_1s(e,3,4,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/NGFSRDW1s4-1s2.txt',   exposant(e), 0, 0);
         rateQuenching_1s(e,3,4,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/NGFSRDW1s4-1s3.txt',   exposant(e), 0, 0);
         rateQuenching_1s(e,3,4,5,i)=rateQuenching_1s(e,3,5,4,i)*exp(abs(energie1s(5)-energie1s(4))/Te(i))*(poids(5)/poids(4));
         rateQuenching_1s(e,3,4,6,i)=rateGround_1s(e,3,4,i)*exp(abs(energie1s(4))/Te(i))*(1/poids(4));
         rateQuenching_1s(e,3,4,7,i)=(rateQuenching_1s(e,3,3,7,i)+rateQuenching_1s(e,3,5,7,i))/2 ; %moyenne des deux   

         rateQuenching_1s(e,3,2,3,i)=rateQuenching_1s(e,3,3,2,i)*exp(abs(energie1s(3)-energie1s(2))/Te(i))*(poids(3)/poids(2));
         rateQuenching_1s(e,3,2,4,i)=rateQuenching_1s(e,3,4,2,i)*exp(abs(energie1s(4)-energie1s(2))/Te(i))*(poids(4)/poids(2));
         rateQuenching_1s(e,3,2,5,i)=rateQuenching_1s(e,3,5,2,i)*exp(abs(energie1s(5)-energie1s(2))/Te(i))*(poids(5)/poids(2));
         rateQuenching_1s(e,3,2,6,i)=rateGround_1s(e,3,2,i)*exp(abs(energie1s(2))/Te(i))*(1/poids(2));
         rateQuenching_1s(e,3,2,7,i)=(rateQuenching_1s(e,3,3,7,i)+rateQuenching_1s(e,3,5,7,i))/2;

        %y3=y1.*exp(abs(energie1s(5)-energie1s(3))./Te)*5 ; %dégé du 1s5=5 dégé du 1s3=1 pour 1s3 à 1s5


    end

    'fin calcul Argon'

    %% ################################################################################################################
    %##################################################################################################################
    %##################################################################################################################
    %############################################# Calcul pour le Neon ################################################
    %##################################################################################################################
    %##################################################################################################################


    for i=1:length(Te)
        % =============== Rate excitation processes ============================
        %.txt file : First column : Energy // Second column : Cross section
        % here I use ratetrapz_TRG function because the corresponding cross section was already in m2 unit unlike others where we need to
        % multiply by 1e-23 in the function

        % ON UTILISE MAINTENANT LES SECTIONS EFFICACES RÉVISÉ PNe ED


        %=====Donnelly data for metastable excitation=====

        rate1s_2p(e,2,5,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,5,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-2p1.txt',   exposant(e), 0, 0);
        rate1s_2p(e,2,3,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-2p1.txt',   exposant(e), 0, 0);

        %=====Lxcat BSR (Quantum-mechanical calculations by O. Zatsarinny and K. Bartschat 2012======
        rateGround_2p(e,2,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p10.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p9.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p8.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p7.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p6.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p5.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p4.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p3.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p2.txt',   exposant(e), 0, 0);
        rateGround_2p(e,2,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-2p1.txt',   exposant(e), 0, 0);
        rateGround_1s(e,2,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-1s3.txt',   exposant(e), 0, 0);
        rateGround_1s(e,2,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/gs-1s5.txt',   exposant(e), 0, 0);
        
        rateGround_2p_HighP(e,2,10,i) = rateGround_2p(e,2,10,i);
        rateGround_2p_HighP(e,2,9,i) = rateGround_2p(e,2,9,i);
        rateGround_2p_HighP(e,2,8,i) = rateGround_2p(e,2,8,i);
        rateGround_2p_HighP(e,2,7,i) = rateGround_2p(e,2,7,i);
        rateGround_2p_HighP(e,2,6,i) = rateGround_2p(e,2,6,i);
        rateGround_2p_HighP(e,2,5,i) = rateGround_2p(e,2,5,i);
        rateGround_2p_HighP(e,2,4,i) = rateGround_2p(e,2,4,i);
        rateGround_2p_HighP(e,2,3,i) = rateGround_2p(e,2,3,i);
        rateGround_2p_HighP(e,2,2,i) = rateGround_2p(e,2,2,i);
        rateGround_2p_HighP(e,2,1,i) = rateGround_2p(e,2,1,i);
        %Donnelly
        %=============================================================
        %======6 pour le ground et 7 pour l'ionisation======
        %=============================================================
        rateQuenching_1s(e,2,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-1s5.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s3-ionization.txt',   exposant(e), 0, 0);

        rateQuenching_1s(e,2,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-1s3.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,2,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Neon/1s5-ionization.txt',   exposant(e), 0, 0);


    end


    'fin calcul Néon'

    %% ################################################################################################################
    %##################################################################################################################
    %##################################################################################################################
    %############################################# Calcul pour le Krypton #############################################
    %##################################################################################################################
    %##################################################################################################################


    for i=1:length(Te)
        % =============== Rate excitation processes ============================
        %.txt file : First column : EKrrgy // Second column : Cross section
        % here I use ratetrapz_TRG function because the corresponding cross section was already in m2 unit unlike others where we Kred to
        % multiply by 1e-23 in the function

        % ON UTILISE MAINTENANT LES SECTIONS EFFICACES RÉVISÉ PKr ED


        %=====Donnelly data for metastable excitation=====

        rate1s_2p(e,4,5,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,5,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-2p1.txt',   exposant(e), 0, 0);
        rate1s_2p(e,4,3,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-2p1.txt',   exposant(e), 0, 0);

        %=====Donnelly data for ground state excitation======
        rateGround_2p(e,4,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p10.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p9.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p8.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p7.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p6.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p5.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p4.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p3.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p2.txt',   exposant(e), 0, 0);
        rateGround_2p(e,4,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-2p1.txt',   exposant(e), 0, 0);
        
        rateGround_2p_HighP(e,4,10,i) = rateGround_2p(e,4,10,i);
        rateGround_2p_HighP(e,4,9,i) = rateGround_2p(e,4,9,i);
        rateGround_2p_HighP(e,4,8,i) = rateGround_2p(e,4,8,i);
        rateGround_2p_HighP(e,4,7,i) = rateGround_2p(e,4,7,i);
        rateGround_2p_HighP(e,4,6,i) = rateGround_2p(e,4,6,i);
        rateGround_2p_HighP(e,4,5,i) = rateGround_2p(e,4,5,i);
        rateGround_2p_HighP(e,4,4,i) = rateGround_2p(e,4,4,i);
        rateGround_2p_HighP(e,4,3,i) = rateGround_2p(e,4,3,i);
        rateGround_2p_HighP(e,4,2,i) = rateGround_2p(e,4,2,i);
        rateGround_2p_HighP(e,4,1,i) = rateGround_2p(e,4,1,i);
        
        rateGround_1s(e,4,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-1s3.txt',   exposant(e), 0, 0);
        rateGround_1s(e,4,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/gs-1s5.txt',   exposant(e), 0, 0);

        %=============================================================
        %======6 pour le ground et 7 pour l'ionisation======
        %=============================================================
        rateQuenching_1s(e,4,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-1s5.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s3-ionization.txt',   exposant(e), 0, 0);

        rateQuenching_1s(e,4,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-1s3.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,4,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Krypton/1s5-ionization.txt',   exposant(e), 0, 0);


    end


    'fin calcul Krypton'

    %##################################################################################################################
    %##################################################################################################################
    %##################################################################################################################
    %############################################# Calcul pour le Xenon ################################################
    %##################################################################################################################
    %##################################################################################################################


    for i=1:length(Te)

        %=====DonXelly data for metastable excitation=====

        rate1s_2p(e,5,5,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p10.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p9.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p8.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p7.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p6.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p5.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p4.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p3.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p2.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,5,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-2p1.txt',   exposant(e), 0, 0);
        rate1s_2p(e,5,3,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-2p1.txt',   exposant(e), 0, 0);

        %=====DonXelly data for ground state excitation======
        rateGround_2p(e,5,10,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p10.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,9,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p9.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,8,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p8.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p7.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p6.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p5.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p4.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p3.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p2.txt',   exposant(e), 0, 0);
        rateGround_2p(e,5,1,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-2p1.txt',   exposant(e), 0, 0);
            
        rateGround_2p_HighP(e,5,10,i) = rateGround_2p(e,5,10,i);
        rateGround_2p_HighP(e,5,9,i) = rateGround_2p(e,5,9,i);
        rateGround_2p_HighP(e,5,8,i) = rateGround_2p(e,5,8,i);
        rateGround_2p_HighP(e,5,7,i) = rateGround_2p(e,5,7,i);
        rateGround_2p_HighP(e,5,6,i) = rateGround_2p(e,5,6,i);
        rateGround_2p_HighP(e,5,5,i) = rateGround_2p(e,5,5,i);
        rateGround_2p_HighP(e,5,4,i) = rateGround_2p(e,5,4,i);
        rateGround_2p_HighP(e,5,3,i) = rateGround_2p(e,5,3,i);
        rateGround_2p_HighP(e,5,2,i) = rateGround_2p(e,5,2,i);
        rateGround_2p_HighP(e,5,1,i) = rateGround_2p(e,5,1,i);
        
        rateGround_1s(e,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-1s3.txt',   exposant(e), 0, 0);
        rateGround_1s(e,5,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/gs-1s5.txt',   exposant(e), 0, 0);

        %=============================================================
        %======6 pour le ground et 7 pour l'ionisation======
        %=============================================================
        rateQuenching_1s(e,5,3,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,3,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,3,5,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-1s5.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,3,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,3,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s3-ionization.txt',   exposant(e), 0, 0);

        rateQuenching_1s(e,5,5,2,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-1s2.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,5,4,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-1s4.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,5,3,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-1s3.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,5,6,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-gs.txt',   exposant(e), 0, 0);
        rateQuenching_1s(e,5,5,7,i)=ratetrapz_TRG(Te(i),'sections_efficaces/Xenon/1s5-ionization.txt',   exposant(e), 0, 0); 

    end

    'fin calcul Xenon'
    toc
disp(strcat(num2str(100*(e/length(exposant))),' % de Fait'))
    
end
%% Save des Datas
Allrates1s_2p = rate1s_2p;
AllratesGround_2p = rateGround_2p;
AllratesGround_2p_HighP = rateGround_2p_HighP;
AllratesGround_1s = rateGround_1s;
AllratesQuenching_1s = rateQuenching_1s;

save('Allrates1s_2p.mat','Allrates1s_2p')
save('AllratesGround_2p.mat','AllratesGround_2p')
save('AllratesGround_2p_HighP.mat','AllratesGround_2p_HighP')
save('AllratesGround_1s.mat','AllratesGround_1s')
save('AllratesQuenching_1s.mat','AllratesQuenching_1s')

 
 
 %du code de Donnely
%% ####################################################################################################
%######################################################################################################
%######################################################################################################
%############################################# Quenching des neutres ###################################
%######################################################################################################
%######################################################################################################
%rate  de quenching des métastables par les neutres (m^3 s^-1) on considère
%que c'est le MEME taux pour tous les niveaux 1s puisque donnelly le fait
%pour les 2 métastables.
%gas to He,Ne,Ar,Kr,Xe
rates_neutral=zeros(25,5);
rates_neutral(1,:)=[0.0, 5.0e-19, 1.0e-19, 1.0e-21, 4.1e-23];%He
rates_neutral(2,:)=[1.0e-17, 0.0, 1.0e-21, 5.0e-21, 1.85e-22];%Ne
rates_neutral(3,:)=[1.3e-16, 1.1e-16, 1.5e-21, 6.0e-22, 6.5e-22];%Ar
rates_neutral(4,:)=[1.0e-16, 1.0e-16, 1.0e-16, 1.0e-21, 1.0e-22];%Kr
rates_neutral(5,:)=[1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-21];%Xe
rates_neutral(13,:)=[0.0, 0.0, 2.1e-16, 1.6e-16, 2.2e-16];%O2
rates_neutral(24,:)=[0.0, 0.0, 2.36e-16, 0.0, 0.0];%HMDSO 
save('Allrates_neutral.mat','rates_neutral')

 'fin du code'
toc
