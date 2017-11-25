function [Yij,Doppler,sig_Yij] = TeEscapeFactor_TRG(gaz,lambda,gi,gj, Aij, nj,rho,Tg,sig_nj,sig_rho)
%% =============================================================================================
%% ========= Cette fonction calcule le facteur d'�chappement � appliquer pour corriger =========
%% ======== l'intensit� des raies auto-absorb�es dans le cas g�n�ral d'un �largissement ========
%% ======== de raie de profil Voigt (r�sonant + Van der Waals + Doppler (Stark n�glig�))========
%% =============================================================================================
%% ====== INFOS SUR LA FONCTION ======
%La variable gaz est le numero associ� au gaz dont la transition est �tudi�e
%La variable lambda est la longueur d'onde de la raie potentiellement auto-absorb�e
%Les variables gi et gj sont les poids statistiques des niveaux impliqu�s dans la transition
%Les variables li et lj sont les nombres quantiques l des niveaux impliqu�s dans la transition
%La variable Aij est le coefficient de transition radiative spontan�e
%La variable nj est la densit� du niveau inf�rieur pouvant r�-absorber la raie
%La variable rho est la longueur sur laquelle peut se produire l'auto-absorption
%Les variables Tg et ng sont la temp�rature et la densit� des neutres � l'�tat fondamental
%La variable fij est la force d'oscillateur d'un niveau 1s avec le fondamental pour le calcul d'�largissement r�sonant
%La variable AllIntegral contient l'information sur le r�sultat de l'int�grale du profil Voigt pour diff�rents beta possibles
%------------------------------------------------------------------------------------

%% R�f�rences, formules tir�es de:
%�largissement Doppler: "MSchulze, A Yanguas-Gil1, A von Keudell and P Awakowicz Center - J. Phys. D: Appl. Phys. 41 (2008) 065206"
%�largissement Van der Waals: "C. Yubero, M.S. Dimitrijevic, M.C. Garc�a, M.D. Calzada - Spectrochimica Acta Part B 62 (2007) 169�176" 
%Convolution Voigt: "Eduardo Casta�os-Mart�nez, Michel Moisan - Spectrochimica Acta Part B 65 (2010) 199�209"
%�largissement r�sonant: formule tir�e du livre Braodening mecanisms chp. 9.2 p. 157 eqn. 9.21

%% ====== D�BUT DE LA FONCTION ======= 
%% D�finitions des constantes
global c

if  gaz==2      %N�on
    M=20.1797;  %Masse en unit� de masse atomique    
elseif gaz==3  %Argon
    M=39.948;       
elseif gaz== 4 %Krypton
    M=83.798;
elseif gaz==5 %Xenon
    M=131.293;
end

%% �largissement Doppler
    Doppler=7.16*10^(-7)*lambda*sqrt(Tg/M); %Doppler en nm
    Doppler = Doppler*1e-9; %Doppler en metre
    sig_Tg=0;
    sig_Doppler = Doppler*(sig_Tg/Tg);
    lambda = lambda*10^-9; % en metre

    %Obtention du coefficient d'absorption global au centre de la raie

    k=(2*lambda^4*sqrt(log(2))*gi*nj*Aij)/(c*Doppler*sqrt(pi)*gj*8*pi);
    
    sig_k = k*sqrt( (sig_Doppler/Doppler)^2 + (sig_nj/nj)^2 );

%% Calcul du facteur d'�chappement
Yij = (2-exp(-(k*rho)*.001))/(1+k*rho);

%changement de variable pour simplifier.
A=k*rho;
sig_A = A*sqrt( (sig_k/k)^2 + (sig_rho/rho)^2);
a=1000;
sig_Yij =abs(exp(-A/a)*((-2*a*exp(A/a)+a+A+1)/(a*(A+1)^2)))*sig_A;

%% plot de Yij en fonction de rho pour Ar r�sonant.
% if lambda == 104.82*10^-9
% L=linspace(0.00001,1,100);
% test = (2-exp(-(k.*L)*.001))./(1+k.*L);
% figure
% plot(L,test)
% end
end