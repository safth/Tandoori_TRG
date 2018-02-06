function [PopFond,DepopMixAutoAbs,PopMixAutoAbs]= TePopulation_Metastable_TRG(ng,ne,rateGround_1s,Thetaij_2p,Aij_2p,Br2p)

%On néglige le peuplement par désexcitation radiatives des 2Px puisqu'elle sont comprise en cascading.
%there is no FUDGE for those rates in donnelly's files
%pop par le fondamental par impact électronique
PopFond = ne*ng*[        0       ;  %RIEN
                 rateGround_1s(2);  %1s2 
                 rateGround_1s(3);  %1s3 
                 rateGround_1s(4);  %1s4
                 rateGround_1s(5)]; %1s5

             
             %% 8) depop autoabsorption
%quand un 1si peuple un sp par autoabs et qu'il se désintégre vers un autre niveau!

AutoAbs=1-Thetaij_2p;
AutoAbs=AutoAbs'; %flip
Aij_2p=Aij_2p'; %flip
DepopMixAutoAbs = [0 ;
                   AutoAbs(2,1)*(1-Br2p(2,1))*Aij_2p(2,1) + AutoAbs(2,2)*(1-Br2p(2,2))*Aij_2p(2,2) + AutoAbs(2,3)*(1-Br2p(2,3))*Aij_2p(2,3) + AutoAbs(2,4)*(1-Br2p(2,4))*Aij_2p(2,4) + AutoAbs(2,5)*(1-Br2p(2,5))*Aij_2p(2,5) + AutoAbs(2,6)*(1-Br2p(2,6))*Aij_2p(2,6) + AutoAbs(2,7)*(1-Br2p(2,7))*Aij_2p(2,7) + AutoAbs(2,8)*(1-Br2p(2,8))*Aij_2p(2,8) + AutoAbs(2,9)*(1-Br2p(2,9))*Aij_2p(2,9) + AutoAbs(2,10)*(1-Br2p(2,10))*Aij_2p(2,10) ;
                   AutoAbs(3,1)*(1-Br2p(3,1))*Aij_2p(3,1) + AutoAbs(3,2)*(1-Br2p(3,2))*Aij_2p(3,2) + AutoAbs(3,3)*(1-Br2p(3,3))*Aij_2p(3,3) + AutoAbs(3,4)*(1-Br2p(3,4))*Aij_2p(3,4) + AutoAbs(3,5)*(1-Br2p(3,5))*Aij_2p(3,5) + AutoAbs(3,6)*(1-Br2p(3,6))*Aij_2p(3,6) + AutoAbs(3,7)*(1-Br2p(3,7))*Aij_2p(3,7) + AutoAbs(3,8)*(1-Br2p(3,8))*Aij_2p(3,8) + AutoAbs(3,9)*(1-Br2p(3,9))*Aij_2p(3,9) + AutoAbs(3,10)*(1-Br2p(3,10))*Aij_2p(3,10) ;
                   AutoAbs(4,1)*(1-Br2p(4,1))*Aij_2p(4,1) + AutoAbs(4,2)*(1-Br2p(4,2))*Aij_2p(4,2) + AutoAbs(4,3)*(1-Br2p(4,3))*Aij_2p(4,3) + AutoAbs(4,4)*(1-Br2p(4,4))*Aij_2p(4,4) + AutoAbs(4,5)*(1-Br2p(4,5))*Aij_2p(4,5) + AutoAbs(4,6)*(1-Br2p(4,6))*Aij_2p(4,6) + AutoAbs(4,7)*(1-Br2p(4,7))*Aij_2p(4,7) + AutoAbs(4,8)*(1-Br2p(4,8))*Aij_2p(4,8) + AutoAbs(4,9)*(1-Br2p(4,9))*Aij_2p(4,9) + AutoAbs(4,10)*(1-Br2p(4,10))*Aij_2p(4,10) ;
                   AutoAbs(5,1)*(1-Br2p(5,1))*Aij_2p(5,1) + AutoAbs(5,2)*(1-Br2p(5,2))*Aij_2p(5,2) + AutoAbs(5,3)*(1-Br2p(5,3))*Aij_2p(5,3) + AutoAbs(5,4)*(1-Br2p(5,4))*Aij_2p(5,4) + AutoAbs(5,5)*(1-Br2p(5,5))*Aij_2p(5,5) + AutoAbs(5,6)*(1-Br2p(5,6))*Aij_2p(5,6) + AutoAbs(5,7)*(1-Br2p(5,7))*Aij_2p(5,7) + AutoAbs(5,8)*(1-Br2p(5,8))*Aij_2p(5,8) + AutoAbs(5,9)*(1-Br2p(5,9))*Aij_2p(5,9) + AutoAbs(5,10)*(1-Br2p(5,10))*Aij_2p(5,10)];
popAuto = zeros(5);
  for n=2:5             %boucle sur le niveux qui peuple le "n"
      for m=2:5         %boucle sur les 1s peuplés
              for k=1:10    %comme sur les 2p                           %Branching Ratio
                popAuto(n)  = popAuto(n) + AutoAbs(m,k)*Aij_2p(m,k)*(Br2p(n,k)); %rate du niv m et Br vers le n
          end
      end
  end  
  
  
    
  
PopMixAutoAbs = [0;
                 popAuto(2) ;
                 popAuto(3) ;
                 popAuto(4) ;
                 popAuto(5)];
