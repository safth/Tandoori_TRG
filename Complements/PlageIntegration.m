%fonction appel� par IntensiteGaussMax_TRG qui prend les point autour d'une
%raie � fitter et qui lisse.
%Prend E(q) en entr�e
function [IntTemp,LambdaTemp,moyenne] = PlageIntegration(E,Lambda,Int,plage,bkgd)
%% D�claration des vvecteurs de sorties
w=1;
j=1; %compte les points � rejeter
k=0; %compte les points � mettre dans la moyenne
%% On scan le spectre � la recherche de la raie
for o=1:length(Lambda)
    if  abs(Lambda(o)-E)<plage      %On regarde � +/-1nm pour �tre certain de bien voir les doublets
                IntTemp(w)=Int(o);          %On enregistre les valeurs d'intensit� de la raie
                LambdaTemp(w)=Lambda(o) ;   %on enregistre aussi les longueurs d'ondes dans la plage permise par facteur*FWHM
                w=w+1;
    end
end
clear o w
 %% Soustraction pour rapporter la baseline � 0 en moyennant sur les 5 premiers et derniers pts
   % qui sont plus petites que le back pour ne pas pogner de raies parasites. 
somme=0;
 for i=1:15
     %check les point  droite si il ne sont pas trop gros
    if abs(IntTemp(i)) < 2*abs(bkgd)   
        somme = somme+IntTemp(i);
        k=k+1;
    else %sinon, ont marque le point dans un vecteur
         a(j) = i;
         j=j+1;
         %check les point � gauche
    end
    if IntTemp(length(IntTemp)-i+1) < 2*abs(bkgd)        
        somme = somme+IntTemp(length(IntTemp)-i+1);
        k=k+1;
        
    else 
        a(j) = length(IntTemp)-i+1;
        j=j+1;
    end
 end

%% enleve la moyenne calcul� et met les points trop haut � la moyenne
    %si aucun point est pogn� (k=o), il faut mettre ca � zero sinon c'est NAN pis ca bug...
if k>0
    moyenne=somme/(k); %%% i = 10 for sure????????????++++++++++++++++++++++++++++++++++++
    IntTemp=IntTemp-moyenne;
    
else
    moyenne = 0.01
    IntTemp=IntTemp-moyenne;
end
 for n=1:j-1
     
     IntTemp(a(n))=moyenne;  
 end
             clear somme i k j n
 
