function [Yij,Doppler,sig_Yij] = TeEscapeFactor_TRG(gaz,lambda,gi,gj, Aij, nj,rho,Tg,sig_nj,sig_rho)
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
global c

if  gaz==2      %Néon
    M=20.1797;  %Masse en unité de masse atomique    
elseif gaz==3  %Argon
    M=39.948;       
elseif gaz== 4 %Krypton
    M=83.798;
elseif gaz==5 %Xenon
    M=131.293;
end

%% Élargissement Doppler
    Doppler=7.16*10^(-7)*lambda*sqrt(Tg/M); %Doppler en nm
    Doppler = Doppler*1e-9; %Doppler en metre
    sig_Tg=0;
    sig_Doppler = Doppler*(sig_Tg/Tg);
    lambda = lambda*10^-9; % en metre

    %Obtention du coefficient d'absorption global au centre de la raie

    k=(2*lambda^4*sqrt(log(2))*gi*nj*Aij)/(c*Doppler*sqrt(pi)*gj*8*pi);
    
    sig_k = k*sqrt( (sig_Doppler/Doppler)^2 + (sig_nj/nj)^2 );

%% Calcul du facteur d'échappement
Yij = (2-exp(-(k*rho)*.001))/(1+k*rho);

%changement de variable pour simplifier.
A=k*rho;
sig_A = A*sqrt( (sig_k/k)^2 + (sig_rho/rho)^2);
a=1000;
sig_Yij =abs(exp(-A/a)*((-2*a*exp(A/a)+a+A+1)/(a*(A+1)^2)))*sig_A;

%% plot de Yij en fonction de rho pour Ar résonant.
% if lambda == 104.82*10^-9
% L=linspace(0.00001,1,100);
% test = (2-exp(-(k.*L)*.001))./(1+k.*L);
% figure
% plot(L,test)
% end
end