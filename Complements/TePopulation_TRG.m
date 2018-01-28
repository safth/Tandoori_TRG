function [PopFond,Pop1s,sig_Pop1s,PopAr]=TePopulation_TRG(i,ng,ne,n1s,rate1s_2p,rateGround_2p,gaz,nm_Ar)


%% Gain de population des niveaux 2p par impact électronique sur un neutre dans son niveau fondamental 

PopFond   =  ne*ng*[rateGround_2p(10,i);  %2p10
                    rateGround_2p(9,i);   %2p9
                    rateGround_2p(8,i);   %2p8
                    rateGround_2p(7,i);   %2p7
                    rateGround_2p(6,i);   %2p6
                    rateGround_2p(5,i);   %2p5
                    rateGround_2p(4,i);   %2p4
                    rateGround_2p(3,i);   %2p3
                    rateGround_2p(2,i);   %2p2
                    rateGround_2p(1,i)];  %2p1
                
                
%% Gain de population des niveaux 2p par impact électronique sur un 1s          
Pop1s = ne*[n1s(2)*rate1s_2p(2,10,i)+ n1s(3)*rate1s_2p(3,10,i)+ n1s(4)*rate1s_2p(4,10,i)+ n1s(5)*rate1s_2p(5,10,i); %2p10
            n1s(2)*rate1s_2p(2,9,i) + n1s(3)*rate1s_2p(3,9,i) + n1s(4)*rate1s_2p(4,9,i) + n1s(5)*rate1s_2p(5,9,i);  %2p9
            n1s(2)*rate1s_2p(2,8,i) + n1s(3)*rate1s_2p(3,8,i) + n1s(4)*rate1s_2p(4,8,i) + n1s(5)*rate1s_2p(5,8,i);  %2p8
            n1s(2)*rate1s_2p(2,7,i) + n1s(3)*rate1s_2p(3,7,i) + n1s(4)*rate1s_2p(4,7,i) + n1s(5)*rate1s_2p(5,7,i);  %2p7
            n1s(2)*rate1s_2p(2,6,i) + n1s(3)*rate1s_2p(3,6,i) + n1s(4)*rate1s_2p(4,6,i) + n1s(5)*rate1s_2p(5,6,i);  %2p6
            n1s(2)*rate1s_2p(2,5,i) + n1s(3)*rate1s_2p(3,5,i) + n1s(4)*rate1s_2p(4,5,i) + n1s(5)*rate1s_2p(5,5,i);  %2p5
            n1s(2)*rate1s_2p(2,4,i) + n1s(3)*rate1s_2p(3,4,i) + n1s(4)*rate1s_2p(4,4,i) + n1s(5)*rate1s_2p(5,4,i);  %2p4
            n1s(2)*rate1s_2p(2,3,i) + n1s(3)*rate1s_2p(3,3,i) + n1s(4)*rate1s_2p(4,3,i) + n1s(5)*rate1s_2p(5,3,i);  %2p3
            n1s(2)*rate1s_2p(2,2,i) + n1s(3)*rate1s_2p(3,2,i) + n1s(4)*rate1s_2p(4,2,i) + n1s(5)*rate1s_2p(5,2,i);  %2p2
            n1s(2)*rate1s_2p(2,1,i) + n1s(3)*rate1s_2p(3,1,i) + n1s(4)*rate1s_2p(4,1,i) + n1s(5)*rate1s_2p(5,1,i)]; %2p1

        
        
%% pour le calcul avec l'incertitude sur les 1s
sig_Pop1s = ne*[0 rate1s_2p(2,10,i) rate1s_2p(3,10,i) rate1s_2p(4,10,i) rate1s_2p(5,10,i); %2p10
                0 rate1s_2p(2,9,i)  rate1s_2p(3,9,i)  rate1s_2p(4,9,i)  rate1s_2p(5,9,i);  %2p9
                0 rate1s_2p(2,8,i)  rate1s_2p(3,8,i)  rate1s_2p(4,8,i)  rate1s_2p(5,8,i);  %2p8
                0 rate1s_2p(2,7,i)  rate1s_2p(3,7,i)  rate1s_2p(4,7,i)  rate1s_2p(5,7,i);  %2p7
                0 rate1s_2p(2,6,i)  rate1s_2p(3,6,i)  rate1s_2p(4,6,i)  rate1s_2p(5,6,i);  %2p6
                0 rate1s_2p(2,5,i)  rate1s_2p(3,5,i)  rate1s_2p(4,5,i)  rate1s_2p(5,5,i);  %2p5
                0 rate1s_2p(2,4,i)  rate1s_2p(3,4,i)  rate1s_2p(4,4,i)  rate1s_2p(5,4,i);  %2p4
                0 rate1s_2p(2,3,i)  rate1s_2p(3,3,i)  rate1s_2p(4,3,i)  rate1s_2p(5,3,i);  %2p3
                0 rate1s_2p(2,2,i)  rate1s_2p(3,2,i)  rate1s_2p(4,2,i)  rate1s_2p(5,2,i);  %2p2
                0 rate1s_2p(2,1,i)  rate1s_2p(3,1,i)  rate1s_2p(4,1,i)  rate1s_2p(5,1,i)]; %2p1

%% population par le transfert d'excitation de l'argon. % nm_Ar(1,1,i,5) i c'Est Te et 5 parceque c'est le 1s5 qui peuple
if gaz == 4
    PopAr = nm_Ar(1,1,i,5)*ng(gaz)*[0        ; %2p10
                                        0        ; %2p9 
                                        0        ; %2p8
                                        0.64e-18 ; %2p7
                                        5.6e-18  ; %2p6 
                                        0        ; %2p5 
                                        0        ; %2p4 
                                        0        ; %2p3 
                                        0        ; %2p2 
                                        0]       ; %2p1
else
    PopAr  =   [0        ; %2p10
                    0        ; %2p9 
                    0        ; %2p8
                    0        ; %2p7
                    0        ; %2p6 
                    0        ; %2p5 
                    0        ; %2p4 
                    0        ; %2p3 
                    0        ; %2p2 
                    0]       ; %2p1
end

end