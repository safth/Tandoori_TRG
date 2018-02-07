function [Depop2p,DepopRadFond,PopAutoAbs,DepopMix,PopMix,DepopSuperelestique,DepopIonisation,DepopNeutre,PopUpDown,PopAr] = ...
    TeDepopulation_Metastable_TRG(i,gaz,ng,ne,Aij1s,Thetaij,rate1s_2p,rateQuenching,rateNeutral,Br2p,nm_Ar)

%% 1) Depop impact électronique vers les 2p
Depop2p = ne*[0 0 0 0 0; %RIEN
              0 rate1s_2p(2,1,i)*(1-Br2p(2,1))+rate1s_2p(2,2,i)*(1-Br2p(2,2))+rate1s_2p(2,3,i)*(1-Br2p(2,3))+rate1s_2p(2,4,i)*(1-Br2p(2,4))+rate1s_2p(2,5,i)*(1-Br2p(2,5))+rate1s_2p(2,6,i)*(1-Br2p(2,6))+rate1s_2p(2,7,i)*(1-Br2p(2,7))+rate1s_2p(2,8,i)*(1-Br2p(2,8))+rate1s_2p(2,9,i)*(1-Br2p(2,9))+rate1s_2p(2,10,i)*(1-Br2p(2,10)) 0 0 0; %1s2->2p1+2p2+...
              0 0 rate1s_2p(3,1,i)*(1-Br2p(3,1))+rate1s_2p(3,2,i)*(1-Br2p(3,2))+rate1s_2p(3,3,i)*(1-Br2p(3,3))+rate1s_2p(3,4,i)*(1-Br2p(3,4))+rate1s_2p(3,5,i)*(1-Br2p(3,5))+rate1s_2p(3,6,i)*(1-Br2p(3,6))+rate1s_2p(3,7,i)*(1-Br2p(3,7))+rate1s_2p(3,8,i)*(1-Br2p(3,8))+rate1s_2p(3,9,i)*(1-Br2p(3,9))+rate1s_2p(3,10,i)*(1-Br2p(3,10)) 0 0; %1s3
              0 0 0 rate1s_2p(4,1,i)*(1-Br2p(4,1))+rate1s_2p(4,2,i)*(1-Br2p(4,2))+rate1s_2p(4,3,i)*(1-Br2p(4,3))+rate1s_2p(4,4,i)*(1-Br2p(4,4))+rate1s_2p(4,5,i)*(1-Br2p(4,5))+rate1s_2p(4,6,i)*(1-Br2p(4,6))+rate1s_2p(4,7,i)*(1-Br2p(4,7))+rate1s_2p(4,8,i)*(1-Br2p(4,8))+rate1s_2p(4,9,i)*(1-Br2p(4,9))+rate1s_2p(4,10,i)*(1-Br2p(4,10)) 0; %1s4
              0 0 0 0 rate1s_2p(5,1,i)*(1-Br2p(5,1))+rate1s_2p(5,2,i)*(1-Br2p(5,2))+rate1s_2p(5,3,i)*(1-Br2p(5,3))+rate1s_2p(5,4,i)*(1-Br2p(5,4))+rate1s_2p(5,5,i)*(1-Br2p(5,5))+rate1s_2p(5,6,i)*(1-Br2p(5,6))+rate1s_2p(5,7,i)*(1-Br2p(5,7))+rate1s_2p(5,8,i)*(1-Br2p(5,8))+rate1s_2p(5,9,i)*(1-Br2p(5,9))+rate1s_2p(5,10,i)*(1-Br2p(5,10))];%1s5

%% 2) Depop Radiative vers le fondamental 
DepopRadFond = [0 0      0  0      0; %1s1 existe pas :(
                0 Aij1s(2) 0  0      0; %1s2 résonnant
                0 0      0  0      0; %1s3 métastable
                0 0      0  Aij1s(4) 0; %1s4 résonnant
                0 0      0  0      0];%1s5 métastable
                 
%% 2.1) Pop par auto-absorption d'une Depop Radiative (1-Thetaij = Pourcentage des photons réabsorbés)
Theta=1-Thetaij;
PopAutoAbs = [0 0               0  0                0; %1s1 existe pas :(
              0 Aij1s(2)*Theta(2) 0  0                0; %1s2 résonnant
              0 0               0  0                0; %1s3 métastable
              0 0               0  Aij1s(4)*Theta(4)  0; %1s4 résonnant
              0 0               0  0                0];%1s5 métastable

%% 3) Depop mixing entre les 1s par impact électronique
DepopMix=ne*[0 0 0 0 0;
             0 rateQuenching(2,3,i)+rateQuenching(2,4,i)+rateQuenching(2,5,i) 0 0 0;  %1s2                                                              %1s2
             0 0 rateQuenching(3,2,i)+rateQuenching(3,4,i)+rateQuenching(3,5,i) 0 0;  %1s3    
             0 0 0 rateQuenching(4,2,i)+rateQuenching(4,3,i)+rateQuenching(4,5,i) 0;  %1s4
             0 0 0 0 rateQuenching(5,2,i)+rateQuenching(5,3,i)+rateQuenching(5,4,i)]; %1s5
                  
%% 3.1) pop mixing entre les 1s par impact électronique                 
PopMix=ne*[0        0                        0                    0                    0          ;
           0        0               rateQuenching(3,2,i) rateQuenching(4,2,i) rateQuenching(5,2,i);
           0 rateQuenching(2,3,i)            0           rateQuenching(4,3,i) rateQuenching(5,3,i);
           0 rateQuenching(2,4,i)   rateQuenching(3,4,i)          0           rateQuenching(5,4,i);
           0 rateQuenching(2,5,i)   rateQuenching(3,5,i) rateQuenching(4,5,i)          0         ];                  

%% 4) Depop collision superélastique des électrons
DepopSuperelestique=ne*[0           0               0               0               0        ;  
                        0    rateQuenching(2,6,i)   0               0               0        ;
                        0           0      rateQuenching(3,6,i)     0               0        ;
                        0           0               0   rateQuenching(4,6,i)        0        ;
                        0           0               0               0   rateQuenching(5,6,i)];
                          
%% 5) Depop ionisation des 1s par impact électronique
DepopIonisation = ne*[0           0               0               0               0        ;  
                      0    rateQuenching(2,7,i)   0               0               0        ;
                      0           0      rateQuenching(3,7,i)     0               0        ;
                      0           0               0   rateQuenching(4,7,i)        0        ;
                      0           0               0               0   rateQuenching(5,7,i)];

%% 6) Depop impact avec un neutre
% Les taux de réaction sont considèrés les memes pour les 2 métastables par
% donnelly. Nous considèrons que c'est les mêmes pour les résonnants aussi.
global nu_hmdso ;

  DepopNeutre = [0 0 0 0 0;
                 0 sum(ng'.*rateNeutral)+nu_hmdso 0 0 0;
                 0 0 sum(ng'.*rateNeutral)+nu_hmdso 0 0;
                 0 0 0 sum(ng'.*rateNeutral)+nu_hmdso 0;
                 0 0 0 0 sum(ng'.*rateNeutral)+nu_hmdso];



%% 7) Pop venant des 2P se désexcitant apres impact avec un autre 1s
  %calcul fait par donnelly lorsqu'il résoud les 1s!!!
  Cascade = zeros(5,5);
  %on fait une somme sur les 10 2p du niveau métastable m vers n
  for n=2:5             %boucle sur le niveux qui peuple le "n"
      for m=2:5         %boucle sur les 1s peuplés
          for k=1:10    %comme sur les 2p                           %Branching Ratio
              Cascade(n,m) = Cascade(n,m) + rate1s_2p(m,k,i) * Br2p(n,k); %rate du niv m et Br vers le n
          end
      end
  end


PopUpDown = ne*[0         0           0               0              0     ; %RIEN
                 0        0        Cascade(2,3)    Cascade(2,4)   Cascade(2,5); %ce qui peuple le 1s2
                 0 Cascade(3,2)        0           Cascade(3,4)   Cascade(3,5); %1s3
                 0 Cascade(4,2)    Cascade(4,3)        0          Cascade(4,5); %1s4
                 0 Cascade(5,2)    Cascade(5,3)    Cascade(5,4)       0      ]; %1s5
             
             
%% 8) Transfert d'excitation pour le Kr de l' Ar(1s5)     
% if gaz == 4
%     PopAr = (nm_Ar(1,1,i,5)+ng(gaz)*1e-16)*ng(gaz)*[0 ; % rajouter n_g et sum(n_metast_Argon) J
%                                     0;
%                                     0;
%                                     5.6e-18*Br2p(4,6)+0.64e-18*Br2p(4,7) ; % 0.64e-18 ; %2p7
%                                     5.6e-18*Br2p(5,6)+0.64e-18*Br2p(5,7)]; % 5.6e-18  ; %2p6   
% 
% elseif gaz == 5
%     %Electronic Energy Transfer from Metastable Argon (4s3P2,0) to Xenon,
%     %Oxygen and Chlorine Atoms. David L king 1973
%    Gain_Ar1s3 = (nm_Ar(1,1,i,3)+ng(gaz)*1e-16)*ng(gaz)*1e-17*[ 10.53*0.0775 ;    %7d->6p
%                                                 4.12*0.03*0.775; %5f->7d->6p
%                                                 0.69*0.505];     %9s->6p
%                     
%     Gain_Ar1s5 = (nm_Ar(1,1,i,5)+ng(gaz)*1e-16)*ng(gaz)*1e-17*[16.78*0.68             ; %8d->6p
%                                                16.78*0.001*0.02*0.78  ; %8d-9p-7d->6p
%                                                16.78*0.001*0.2*0.51   ; %8d-9p-9s->6p
%                                                8.82*0.022             ; %6f-8d->6p
%                                                8.82*0.007*0.775       ; %6f-7d->6p
%                                                2.96*0.021*0.775       ; %9p-7d->6p
%                                                2.96*0.188*0.505       ; %9p-9s->6p
%                                                1.42*0.47              ; %10s->6p
%                                                1.42*0.021*0.775       ; %10s-9p-7d->6p
%                                                1.42*0.188*0.505]      ; %10s-9p-9s->6p 
%                                            %tot 3e-16
%    Gain_tot = (sum(Gain_Ar1s3) + sum(Gain_Ar1s5))/6;        %/6 car c'est des Br moyen vers tous les 2p10à5             
%     
%    PopAr = Gain_tot*[0;
%                      sum(Br2p(2,5:10));
%                      sum(Br2p(3,5:10));
%                      sum(Br2p(4,5:10));
%                      sum(Br2p(5,5:10))];
%    
% else
    PopAr=[0;
           0;
           0;
           0;
           0];
    
%end

end                                                                 
                                        