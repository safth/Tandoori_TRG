%% Cette fonction n'est pas utilis�, elle vient du code d'argon pression atmosph�rique

function [Yij,Doppler,VdW,Resonant] = TeEscapeFactor_TRG(gaz,lambda, gi,li,EnergieI, gj,lj,EnergieJ, Aij, nj,rho,Tg,ng,fji,AllIntegral)
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
c = 299792458; %vitesse de la lumi�re
Eh=13.56;   %Energie d'ionisation de l'hydrog�ne

if  gaz==2      %N�on
    M=20.1797;  %Masse en unit� de masse atomique
    Eip=21.56;  %�nergie d'ionisation en eV
    alpha=0.392*10^(-24); %Polarisabilit� 
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

%% �largissement Doppler
    Doppler=7.16*10^(-7)*lambda*sqrt(Tg/M); %Doppler en nm

%% �largissement van der Waals
    %Niveau inf�rieur
    l=lj; EL=EnergieJ; 
    netoileL=Eh/(Eip-EL);
    RL=netoileL/2+(5*netoileL+1-3*l*(l+1));
    
    %Niveau sup�rieur
    l=li; EU=EnergieI;
    netoileU=Eh/(Eip-EU);
    RU=netoileU/2+(5*netoileU+1-3*l*(l+1));

    R=RU-RL;
    mu=(M/2);
    VdW=8.18*10^(-26)*(lambda^2)*((alpha*R)^(2/5))*((Tg/mu)^(3/10))*(ng/10^6);  %DVdW en cm
    VdW=10^(7)*VdW; %Conversion de cm a nm 

%% �largissement R�sonant
    %Poids statistiques du niveau fondamental
    gground=1;
    %Poids statistique des niveaux absorbeurs pour chaque raie
    g1s=gj;
    %Passage en fr�quence angulaire dela transition r�sonant-fondamental
    omegaTheo=2*pi*c./(lambdaRes*10^(-9));
    %Calcul de l'�largissement en fr�quence angulaire
    Resonant=fji*4*pi*sqrt(gground/g1s)*2.81794*10^(-15)*(c^2/omegaTheo)*ng ;
    %Retour � un �largissement en longueur d'onde
    Resonant = 10^9*((lambda*10^-9)^2)*(Resonant)/(2*pi*c);     %Resonant en nm

%% Convolution en profil Voigt   
    Lorentz = VdW + Resonant;
    beta=(Lorentz/Doppler)*sqrt(log(2));
    %R�sultat de la convolution via la variable AllIntegral d�j� pr�charg�e dans le code principal 
    position = round(beta./0.001);       %Pour AllIntegral_new le step est de 0.001
   % [valeur position]=(min(abs(beta./AllIntegral(:,1)-1)));
    if position<1
        position=1;
    end
    %clear valeur
    integrale = AllIntegral(position,2); %Trouve la valeur correspondante � la position de la valeur de beta dans la colonne des valeurs d'integrale

    %Conversion d'unit�s
    lambda = lambda*10^-9; % en metre
    Dopplerfreq = (c/lambda^2)*(Doppler*10^-9); % conversion du Doppler en fr�quence (s-1)

    %Obtention du coefficient d'absorption global au centre de la raie
    k=(nj/beta)*(sqrt(log(2))/4)*(gi/gj)*(lambda^2*Aij/(Dopplerfreq*(integrale)));

clear  integrale position Dopplerfreq     

%% Calcul du facteur d'�chappement
Yij = (2-exp(-(k*rho)*.001))/(1+k*rho);
end