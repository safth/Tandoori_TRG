function [PopFond]= TePopulation_Metastable_TRG(ng,ne,rateGround_1s)

%On néglige le peuplement par désexcitation radiatives des 2Px puisqu'elle sont comprise en cascading.
%there is no FUDGE for those rates in donnelly's files
%pop par le fondamental par impact électronique
PopFond = ne*ng*[        0       ;  %RIEN
                 rateGround_1s(2);  %1s2 
                 rateGround_1s(3);  %1s3 
                 rateGround_1s(4);  %1s4
                 rateGround_1s(5)]; %1s5

