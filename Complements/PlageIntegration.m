%fonction appelé par IntensiteGaussMax_TRG qui prend les point autour d'une
  %raie à fitter et qui lisse.
%Prend E(q) en entrée
function [IntTemp,LambdaTemp,moyenne] = PlageIntegration(E,Lambda,Int,plage,bkgd)
%% Déclaration des vvecteurs de sorties
w=1;
j=1; %compte les points à rejeter
k=0; %compte les points à mettre dans la moyenne
%% On scan le spectre à la recherche de la raie
for o=1:length(Lambda)
    if  abs(Lambda(o)-E)<plage      %On regarde à +/-1nm pour être certain de bien voir les doublets
                IntTemp(w)=Int(o);          %On enregistre les valeurs d'intensité de la raie
                LambdaTemp(w)=Lambda(o) ;   %on enregistre aussi les longueurs d'ondes dans la plage permise par facteur*FWHM
                w=w+1;
    end
end
clear o w
 %% Soustraction pour rapporter la baseline à 0 en moyennant sur les 5 premiers et derniers pts
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
         %check les point è gauche
    end
    if IntTemp(length(IntTemp)-i+1) < 2*abs(bkgd)        
        somme = somme+IntTemp(length(IntTemp)-i+1);
        k=k+1;
        
    else 
        a(j) = length(IntTemp)-i+1;
        j=j+1;
    end
 end

%% enleve la moyenne calculé et met les points trop haut à la moyenne
    %si aucun point est pogné (k=o), il faut mettre ca à zero sinon c'est NAN pis ca bug...
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
 
