function [En2px_1s,En2px_2py,poids2p,poids1s,LambdaTheo2p,LambdaTheo1s,Aij2p,Aij1s,energie1s,energie2p,fji2p,fji1s,Br2p]=TeInitialisation_TRG(gaz)

%Nomenclature des gaz à la Donnelly:
% 2 = Ne ; 3 = Ar ;  4 = Kr ; 5 = Xe

%% Définition du poids statistique des niveaux
poids2p=[1 3 5 3 1 5 3 5 7 3];  %2p1 ... 2p10
poids1s=[0 3 1 3 5];            %0 1s2 ..1s5

%% Définition des taux de transfert d'énergie par collision (Energy transfer quenching) des niveaux 2p (en m3/sec)
En2px_1s=zeros(1,10);                
En2px_2py=zeros(10,10);

%% ###################################################################
%############################### Néon ################################
%#####################################################################
if gaz == 2       
    
%% Definition de l'énergie des niveaux 
energie2p=[18.966 18.726 18.711 18.704 18.693 18.637 18.613 18.576 18.555 18.382];

energie1s=[0 16.848 16.715 16.671 16.619]; 
   
%% Définition des transitions théorique 2p->1s. Source: Nist
LambdaTheo2p =  [0   585.2488         0  540.0562         0; %2P1 ->1sy
                 0   659.8953  616.3594  602.9997  588.1895; %2P2 ->1sy
                 0   665.2092         0  607.4338         0; %2P5 ->1sy
                 0   667.8277         0  609.6163  594.4834; %2P3 ->1sy           
                 0   671.7043  626.6495  612.8450  597.5534; %2P4 ->1sy       
                 0   692.9467         0  630.4789  614.3063; %2P6 ->1sy
                 0   702.4050  653.2882  638.2991  621.7281; %2P7 ->1sy
                 0   717.3938         0  650.6528  633.4428; %2P8 ->1sy
                 0   0                0         0  640.2248; %2P9 ->1sy
                 0   808.2458  743.8898  724.5167  703.2413];%2P10->1sy
             
 LambdaTheo1s = [0   73.58960         0   74.3719         0];            
              
%% Source: NIST
 Aij2p =  [0 68200000           0      900000           0;%1
           0 23200000    14600000     5610000    11500000;%2
           0 290000             0    60300000           0;%5
           0 23300000           0    18100000    11300000;%3
           0 21700000    24900000      670000     3510000;%4
           0 17400000           0     4160000    28200000;%6
           0 1890000     10800000    32100000     6370000;%7
           0 2870000            0    30000000    16100000;%8
           0 0                  0           0    51400000;%9
           0 120000       2310000     9350000    25300000];%10
             

Aij1s  =  [0 611000000          0    47600000           0]; 

% Définition des Forces d'oscillateurs pour les raies théoriques équivalentes       
%Force d'oscillateur correspondants entre le ground state et le niveau le plus bas de la transition associé à cette raie;  
%soit 1s2 ou 1s4. 0 pour les niveaux non résonants ce qui donnera un élargissement résonant nul pour ces raies.        %%Data NIST
fji2p    =  [0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0.149    0      0.0118           0;
             0  0        0           0           0;
             0  0.149    0      0.0118           0];
         
fji1s   =   [0  0.149    0      0.0118           0];

for j=1:5      
    for i=1:10
        Br2p(j,i) = Aij2p(i,j)/sum(Aij2p(i,:)); %Br2p(vers1s,du 2p)
    end
end

end   

%% ####################################################################
%############################### Argon ################################
%######################################################################
if gaz == 3

% Definition de l'énergie des niveaux
energie2p=[13.48 13.33 13.30 13.28 13.27 13.17 13.15 13.09 13.08 12.91]; %2p1 à 2p10
energie1s=[0 11.828 11.723 11.624 11.549]; %1s2 à 1s5

% Définition des transitions théorique 2p->1s et 1s-ground
LambdaTheo2p  = [0    750.39  0       667.73       0;          % 0..0 2p1->1s2 2p1->1s3 2p1->1s4 2p1->1s5    
                 0    826.45  772.42  727.29  696.54;
                 0    840.82  0       738.4   706.72;
                 0    852.14  794.82  0       714.70;
                 0    0       0       751.47       0;
                 0    922.45  0       800.62  763.51;
                 0    935.42  866.79  810.37  772.38;
                 0    978.45  0       842.47  801.48;
                 0    0       0       0       811.53;
                 0    1148.8  1047    965.78  912.30];         % 0..0 2p10->1s2 2p10->1s3 2p10->1s4 2p10->1s5
             
LambdaTheo1s=   [0    104.82  0       106.67  0     ];         % RIEN, 1s2->fondamental ... 1s5->fondamental
        



% Définition du Aij des transitions sans tenir compte du piégeage (en s-1)
Aij2p =[0   4.45e7  0        2.36e5        0;      % 0..0 2p1->1s2 2p1->1s3 2p1->1s4 2p1->1s5    
        0   1.53e7  1.17e7   1.83e6   6.39e6;
        0   2.23e7  0        8.47e6   3.8e6;
        0   1.39e7  1.86e7   0        6.25e5;
        0   0       0        4.02e7        0;
        0   5.03e6  0        4.9e6    2.45e7;
        0   1e6     2.43e6   2.5e7    5.18e6;
        0   1.47e6  0        2.15e7   9.28e6;
        0   0       0        0        3.31e7;
        0   1.9e5   9.8e5    5.43e6   1.89e7];     % 0..0 2p10->1s2 2p10->1s3 2p10->1s4 2p10->1s5
    
%Ajustement des Aij pour correspondre aux branching Ratio de Donnelly Teoes
%2007
%826 0.45 et non 0.43
fac=0.94;
Aij2p(2,3:5)=Aij2p(2,3:5)*fac; 
%840 0.67 aulieu de 0.64
fac=0.895;
Aij2p(3,3:5)=Aij2p(3,3:5)*fac; 
%852 0.43 aulieu de 0.42
%794 0.55 aulieu de 0.56
fac=0.95;
Aij2p(4,3:5)=Aij2p(4,3:5)*fac; 
%751 0.95 aulieu de 100%
Aij2p(5,4)=Aij2p(5,4)*(0.95); 
%763 0.7 aulieu de 0.71
%800.6 0.13 aulieu de 0.15

%Aij2p(6,2)=Aij2p(6,2)*1.2; 
%Aij2p(6,5)=Aij2p(6,5)*1.07;

%810.4 0.72 aulieu de 0.74
%866 0.083 aulieu de 0.072
Aij2p(7,2)=Aij2p(7,2)*1.1; 
Aij2p(7,5)=Aij2p(7,5)*1.1;
Aij2p(7,3)=Aij2p(7,3)*1.2;
%965 0.18 aulieu de 0.21
fac=1.23;
Aij2p(10,2)=Aij2p(10,2)*fac; 
Aij2p(10,3)=Aij2p(10,3)*fac; 
Aij2p(10,5)=Aij2p(10,5)*fac; 
%les autres raies sont correct, bon Branching Ratio.
    
    

        
Aij1s= [0   5.10e8  0        1.19e8        0];   
 

% Définition des Forces d'oscillateurs pour les raies théoriques équivalentes       
%Force d'oscillateur correspondants entre le ground state et le niveau le plus bas de la transition associé à cette raie;  
%soit 1s2 ou 1s4. 0 pour les niveaux non résonants ce qui donnera un élargissement résonant nul pour ces raies.
fji2p  =[0   0.25     0   0.0609   0;
         0   0.25     0   0.0609   0;
         0   0.25     0   0.0609   0;
         0   0.25     0   0        0;
         0   0        0   0.0609   0;
         0   0.25     0   0.0609   0;
         0   0.25     0   0.0609   0;
         0   0.25     0   0.0609   0;
         0   0        0   0        0;
         0   0.25     0   0.0609   0];
        
fji1s= [0   0.25      0   0.0609   0]; % RIEN, 1s2->fondamental ... 1s5->fondamental

for j=1:5      
    for i=1:10
        Br2p(j,i) = Aij2p(i,j)/sum(Aij2p(i,:)); %Br2p(vers1s,du 2p)
    end
end


%% Seulement pour l'argon car la littérature le permet
% Définition des taux de transfert d'énergie par collision (Energy transfer quenching) des niveaux 2p (en m3/sec)
% Collision 2p-1s: A simple collisional–radiative model for low-temperature argon discharges with pressure ranging from 1?Pa to atmospheric pressure:
% kinetics of Paschen 1s and 2p levels,Journal of Physics D: Applied Physics,Zhu, Xi-Ming Pu, Yi-Kang,2009
% Mis à part 2p5 et 2p6 qui eux viennent de pas J. Chem. Phys. 69(9), (1978) 3885. The only modification in this is I used quenching for 2p6 level 3.0 instead of zero
% reported in the paper.For  2p5 levels quenching has been taken from  J. Chem. Phys. 76, 977 (1982) (sadeghi p980)
En2px_1s=1e-17*[3.0 5.3 4.7 3.9 2.0 2.2 6.1 3.0 3.5 2.0];                

% 13009 Nguyen, T. D., and Sadeghi, N., Rate Coefficients for Collisional Population Transfer between 3p54p Argon Levels at 300°K, Phys. Rev. A 18, 1388-1395 (1978)
load En2px_2py.mat
end



%% ####################################################################
%############################### Krypton ################################
%######################################################################
if gaz == 4
    
% Definition de l'énergie des niveaux
energie2p=[12.2662 12.1533 12.1500 12.1099 11.6753 11.5550 11.5352 11.4537 11.4521 11.3124];
energie1s=[0 10.6521 10.5708 10.0403 9.9231];

% Définition des transitions théorique 2p->1s et 1s-ground (Nist)
LambdaTheo2p  = [0  768.5246         0   557.3130         0; %  2p1->1s2 2p1->1s3 2p1->1s4 2p1->1s5
                 0  826.3243         0   587.0916  556.2225;
                 0  828.1052  785.4823   587.9900  557.0289;
                 0  850.8873  805.9505   599.3850  567.2451;
                 0         0         0   758.7414         0;
                 0         0         0   819.0057  760.1546;
                 0         0         0   829.8110  769.4540;
                 0         0         0   877.6751  810.4366;
                 0         0         0          0  811.2901;
                 0         0         0   975.1761  892.8693];
             
LambdaTheo1s =  [0  116.4867         0   123.5838         0];

            
    %Prise du Nist et du atomic table...(comme Xe) pour celles qu'il manquait.
 Aij2p  =   1e6*[0   4      0        0.44         0;
                 0   35     0        1.8       0.28;
                 0   19     23       0.2        2.1;
                 0   24     19       0.08     0.005;
                 0   0      0        51           0;
                 0   0      0        11          31;
                 0   0      0        32         5.6;
                 0   0      0        27          13;
                 0   0      0        0           36;
                 0   0      0        3.4         37];

%Ajustement des Aij pour correspondre aux branching Ratio de Donnelly Teoes
%2007
%768 100% 
Aij2p(1,4)= 0;  
%785 0.55 aulieu de 0.53
fac=0.9;
Aij2p(3,2)=Aij2p(3,2)*fac;
Aij2p(3,4)=Aij2p(3,4)*fac;
Aij2p(3,5)=Aij2p(3,5)*fac;
 %806 0.42 aulieu de 0.44
fac=18.00;
Aij2p(4,2)=Aij2p(4,2);%850 0.54 aulieu de 0.56 2p4 aussi
Aij2p(4,4)=Aij2p(4,4)*fac;
Aij2p(4,5)=Aij2p(4,5)*fac;
%760 0.76 aulieu de 0.74
fac=0.89;
Aij2p(6,4)=Aij2p(6,4)*fac;
%769 0.14 aulieu de 0.15
fac=1.08;
Aij2p(7,4)=Aij2p(7,4)*fac;%829 0.86 aulieu de 0.85
%810 0.34 aulieu de 0.33
fac=0.965;
Aij2p(8,4)=Aij2p(8,4)*fac;%877 è 0.66 aulieu de 0.68
%892 0.9 aulieu de 0.92 
fac=1.2;
Aij2p(10,4)=Aij2p(10,4)*fac;
             
             
             
             
Aij1s    =  1e6*[0   290    0        257          0];   

           % Valeur du Atomic data table en "lenght" correct selon toutes
           % les valeurs de "oscillator strengths of Kr I and Xe I
           % resonance lines (1995) J.C. Molino Garcia
 fji2p    =     [0   0.18   0        0.18     0;
                 0   0.18   0        0.18     0;
                 0   0.18   0        0.18     0;
                 0   0.18   0        0.18     0;
                 0   0      0        0.18     0;
                 0   0      0        0.18     0;
                 0   0      0        0.18     0;
                 0   0      0        0.18     0;
                 0   0      0        0        0;
                 0   0      0        0.18     0];            
   
fji1s     =     [0   0.18   0        0.18     0];         

for j=1:5      
    for i=1:10
        Br2p(j,i) = Aij2p(i,j)/sum(Aij2p(i,:)); %Br2p(vers1s,du 2p)
    end
end

end 

%% #################################################################################
%############################### Xenon #############################################
%####################################################################################
if gaz == 5

% Definition de l'énergie des niveaux
energie2p=[11.1500 11.0779 11.0635 10.9663 9.9414 9.8289 9.7971 9.7284 9.6933 9.5877];
energie1s=[0 9.5773 9.4547 8.4432 8.3219];

% Définition des transitions théorique 2p->1s et 1s->gs source: NIST
LambdaTheo2p  =  [0   788.7393          0  458.2747         0;
                  0   826.6520  764.2024   470.8210  450.0978;
                  0   834.6822          0  473.4152  452.4681;
                  0   893.0830  820.6336   491.6507  469.0970;
                  0           0         0  828.0116         0;
                  0           0         0  895.2251  823.1634;
                  0           0         0  916.2652  840.9189;
                  0           0         0         0  881.9411;
                  0           0         0  992.3198  904.5447;
                  0           0         0         0  979.9697];
    
LambdaTheo1s = [0   129.5588          0  146.9610         0];


%Valeurs de Donnely, dans le theorical transition probabilities and lifetimes in Kr I and Xe I spectra de M. Aymar et
%M.couloube (1978) Atomic Data and nuclear Data table 21, 537-566
Aij2p   =  1e6*[0   0.22e2    0          0.30e1           0;
                0   0.13e2    0.15e2     0.45e-1       0.16;
                0   0.26e2    0          0.11e1        0.87;
                0   0.10e2    0.93e1     0.12e1        0.39;
                0   0         0          0.45e2           0;   
                0   0         0          0.11e2      0.25e2;
                0   0         0          0.32e2      0.18e1;
                0   0         0          0           0.39e2;
                0   0         0          0.20e2      0.11e2;
                0   0         0          0           0.31e2];
            


%Ajustement des Aij pour correspondre aux branching Ratio de Donnelly Teoes
%2007
% pour matcher les Br de Donnelly
%2P1 788  Br à 0.72
fac=2.85;
Aij2p(1,4)=Aij2p(1,4)*fac; 
%764 0.5 aulieu de 0.53
fac=1.13;
Aij2p(2,2)=Aij2p(2,2)*fac; 
Aij2p(2,4)=Aij2p(2,4)*fac; 
Aij2p(2,5)=Aij2p(2,5)*fac; 
%2P3 834 Br à 0.86 et non 0.93
fac=2;
Aij2p(3,4)=Aij2p(3,4)*fac;
Aij2p(3,5)=Aij2p(3,5)*fac;
%820 0.41 et pas 0.45
fac=1.155;
Aij2p(4,2)=Aij2p(4,2)*fac; 
Aij2p(4,4)=Aij2p(4,4)*fac; 
Aij2p(4,5)=Aij2p(4,5)*fac; 
%823 0.75 et non 0.7
fac=0.758;
Aij2p(6,4)=Aij2p(6,4)*fac; 
%840 0.085 et non 0.053
fac=0.605;
Aij2p(7,4)=Aij2p(7,4)*fac; 
%904 0.38 et non 0.35
fac=0.9;
Aij2p(9,4)=Aij2p(9,4)*fac; 

            
            
Aij1s   =  1e6*[0   0.41e3    0          0.29e3           0];        
             
                %Data fu NIST
fji2p    =     [0   0.186    0     0.273       0;
                0   0.186    0     0.273       0;
                0   0.186    0     0.273       0;
                0   0.186    0     0.273       0;
                0   0        0     0.273       0;   
                0   0        0     0.273       0;
                0   0        0     0.273       0;
                0   0        0     0           0;
                0   0        0     0.273       0;
                0   0        0     0           0];   

fji1s    =     [0   0.186    0     0.273       0]; 

for j=1:5      
    for i=1:10
        Br2p(j,i) = Aij2p(i,j)/sum(Aij2p(i,:)); %Br2p(vers1s,du 2p)
    end
end

end
end