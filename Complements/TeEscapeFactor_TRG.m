%% Cette fonction n'est pas utilisé, elle vient du code d'argon pression atmosphérique

function [Yij,Doppler,VdW,Resonant] = TeEscapeFactor_TRG(gaz,lambda, gi,li,EnergieI, gj,lj,EnergieJ, Aij, nj,rho,Tg,ng,fji,AllIntegral)
%% =============================================================================================
%% ========= Cette fonction calcule le facteur d'échappement à appliquer pour corriger =========
%% ======== l'intensité des raies auto-absorbées dans le cas général d'un élargissement ========
%% ======== de raie de profil Voigt (résonant + Van der Waals + Doppler (Stark négligé))========
%% =============================================================================================

%% ====== INFOS SUR LA FONCTION ======
%La variable gaz est le numero associé au gaz dont la transition est étudiée
%La variable lambda est la longueur d'onde de la raie potentiellement auto-absorbée
%Les variables gi et gj sont les poids statistiques des niveaux impliqués dans la transition
%Les variables li et lj sont les nombres quantiques l des niveaux impliqués dans la transition
%La variable Aij est le coefficient de transition radiative spontanée
%La variable nj est la densité du niveau inférieur pouvant ré-absorber la raie
%La variable rho est la longueur sur laquelle peut se produire l'auto-absorption
%Les variables Tg et ng sont la température et la densité des neutres à l'état fondamental
%La variable fij est la force d'oscillateur d'un niveau 1s avec le fondamental pour le calcul d'élargissement résonant
%La variable AllIntegral contient l'information sur le résultat de l'intégrale du profil Voigt pour différents beta possibles
%------------------------------------------------------------------------------------

%% Références, formules tirées de:
%Élargissement Doppler: "MSchulze, A Yanguas-Gil1, A von Keudell and P Awakowicz Center - J. Phys. D: Appl. Phys. 41 (2008) 065206"
%Élargissement Van der Waals: "C. Yubero, M.S. Dimitrijevic, M.C. García, M.D. Calzada - Spectrochimica Acta Part B 62 (2007) 169–176" 
%Convolution Voigt: "Eduardo Castaños-Martínez, Michel Moisan - Spectrochimica Acta Part B 65 (2010) 199–209"
%Élargissement résonant: formule tirée du livre Braodening mecanisms chp. 9.2 p. 157 eqn. 9.21

%% ====== DÉBUT DE LA FONCTION ======= 
%% Définitions des constantes
c = 299792458; %vitesse de la lumière
Eh=13.56;   %Energie d'ionisation de l'hydrogène

if  gaz==2      %Néon
    M=20.1797;  %Masse en unité de masse atomique
    Eip=21.56;  %Énergie d'ionisation en eV
    alpha=0.392*10^(-24); %Polarisabilité 
    lambdaRes=74;   %Longueur d'onde de la transition resonante en nm
    
elseif gaz==3  %Argon
    M=39.948;      
    Eip=15.76;     
    alpha=1.654*10^(-24);   
    lambdaRes=105;
    
elseif gaz== 4 %Krypton
    M=83.798;
    Eip=14.00;
    alpha=2.465*10^(-24);
    lambdaRes=119; 
    
elseif gaz==5 %Xenon
    M=131.293;
    Eip=12.13;
    alpha=4.01*10^(-24);
    lambdaRes=137;
end

%% Élargissement Doppler
    Doppler=7.16*10^(-7)*lambda*sqrt(Tg/M); %Doppler en nm

%% Élargissement van der Waals
    %Niveau inférieur
    l=lj; EL=EnergieJ; 
    netoileL=Eh/(Eip-EL);
    RL=netoileL/2+(5*netoileL+1-3*l*(l+1));
    
    %Niveau supérieur
    l=li; EU=EnergieI;
    netoileU=Eh/(Eip-EU);
    RU=netoileU/2+(5*netoileU+1-3*l*(l+1));

    R=RU-RL;
    mu=(M/2);
    VdW=8.18*10^(-26)*(lambda^2)*((alpha*R)^(2/5))*((Tg/mu)^(3/10))*(ng/10^6);  %DVdW en cm
    VdW=10^(7)*VdW; %Conversion de cm a nm 

%% Élargissement Résonant
    %Poids statistiques du niveau fondamental
    gground=1;
    %Poids statistique des niveaux absorbeurs pour chaque raie
    g1s=gj;
    %Passage en fréquence angulaire dela transition résonant-fondamental
    omegaTheo=2*pi*c./(lambdaRes*10^(-9));
    %Calcul de l'élargissement en fréquence angulaire
    Resonant=fji*4*pi*sqrt(gground/g1s)*2.81794*10^(-15)*(c^2/omegaTheo)*ng ;
    %Retour à un élargissement en longueur d'onde
    Resonant = 10^9*((lambda*10^-9)^2)*(Resonant)/(2*pi*c);     %Resonant en nm

%% Convolution en profil Voigt   
    Lorentz = VdW + Resonant;
    beta=(Lorentz/Doppler)*sqrt(log(2));
    %Résultat de la convolution via la variable AllIntegral déjà préchargée dans le code principal 
    position = round(beta./0.001);       %Pour AllIntegral_new le step est de 0.001
   % [valeur position]=(min(abs(beta./AllIntegral(:,1)-1)));
    if position<1
        position=1;
    end
    %clear valeur
    integrale = AllIntegral(position,2); %Trouve la valeur correspondante à la position de la valeur de beta dans la colonne des valeurs d'integrale

    %Conversion d'unités
    lambda = lambda*10^-9; % en metre
    Dopplerfreq = (c/lambda^2)*(Doppler*10^-9); % conversion du Doppler en fréquence (s-1)

    %Obtention du coefficient d'absorption global au centre de la raie
    k=(nj/beta)*(sqrt(log(2))/4)*(gi/gj)*(lambda^2*Aij/(Dopplerfreq*(integrale)));

clear  integrale position Dopplerfreq     

%% Calcul du facteur d'échappement
Yij = (2-exp(-(k*rho)*.001))/(1+k*rho);
end