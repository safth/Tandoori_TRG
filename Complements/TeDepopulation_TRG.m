function [DepopRad,PopAutoAbsRad,DepopColl2p2p,DepopColl1s,PopColl2p]=TeDepopulation_TRG(Aij,Thetaij,En2px_1s,En2px_2py,ng)

%% Depopulation due à une transition radiative 
%Aij change dans pour chaque gaz dans Teinitialisation
DepopRad=[(Aij(10,5)+Aij(10,4)+Aij(10,3)+Aij(10,2)),0,0,0,0,0,0,0,0,0; %2p10
           0,(Aij(9,5)+Aij(9,4)+Aij(9,3)+Aij(9,2)),0,0,0,0,0,0,0,0;    %2p9
           0,0,(Aij(8,5)+Aij(8,4)+Aij(8,3)+Aij(8,2)),0,0,0,0,0,0,0;    %2p8
           0,0,0,(Aij(7,5)+Aij(7,4)+Aij(7,3)+Aij(7,2)),0,0,0,0,0,0;    %2p7
           0,0,0,0,(Aij(6,5)+Aij(6,4)+Aij(6,3)+Aij(6,2)),0,0,0,0,0;    %2p6
           0,0,0,0,0,(Aij(5,5)+Aij(5,4)+Aij(5,3)+Aij(5,2)),0,0,0,0;    %2p5
           0,0,0,0,0,0,(Aij(4,5)+Aij(4,4)+Aij(4,3)+Aij(4,2)),0,0,0;    %2p4
           0,0,0,0,0,0,0,(Aij(3,5)+Aij(3,4)+Aij(3,3)+Aij(3,2)),0,0;    %2p3
           0,0,0,0,0,0,0,0,(Aij(2,5)+Aij(2,4)+Aij(2,3)+Aij(2,2)),0;    %2p2
           0,0,0,0,0,0,0,0,0,(Aij(1,5)+Aij(1,4)+Aij(1,3)+Aij(1,2))];   %2p1

%test
%% Population due au piégeage d'une transition radiative
AutoAbs=1-Thetaij; %Si Thetaij est la proportion de photons qui sortent, 1-Thetaij est la porportion d'autoabsorbés       

PopAutoAbsRad=[(Aij(10,5)*AutoAbs(10,5)+Aij(10,4)*AutoAbs(10,4)+Aij(10,3)*AutoAbs(10,3)+Aij(10,2)*AutoAbs(10,2)),0,0,0,0,0,0,0,0,0; %2p10
                0,(Aij(9,5)*AutoAbs(9,5)+Aij(9,4)*AutoAbs(9,4)+Aij(9,3)*AutoAbs(9,3)+Aij(9,2)*AutoAbs(9,2)),0,0,0,0,0,0,0,0;    %2p9
                0,0,(Aij(8,5)*AutoAbs(8,5)+Aij(8,4)*AutoAbs(8,4)+Aij(8,3)*AutoAbs(8,3)+Aij(8,2)*AutoAbs(8,2)),0,0,0,0,0,0,0;    %2p8
                0,0,0,(Aij(7,5)*AutoAbs(7,5)+Aij(7,4)*AutoAbs(7,4)+Aij(7,3)*AutoAbs(7,3)+Aij(7,2)*AutoAbs(7,2)),0,0,0,0,0,0;    %2p7
                0,0,0,0,(Aij(6,5)*AutoAbs(6,5)+Aij(6,4)*AutoAbs(6,4)+Aij(6,3)*AutoAbs(6,3)+Aij(6,2)*AutoAbs(6,2)),0,0,0,0,0;    %2p6
                0,0,0,0,0,(Aij(5,5)*AutoAbs(5,5)+Aij(5,4)*AutoAbs(5,4)+Aij(5,3)*AutoAbs(5,3)+Aij(5,2)*AutoAbs(5,2)),0,0,0,0;    %2p5
                0,0,0,0,0,0,(Aij(4,5)*AutoAbs(4,5)+Aij(4,4)*AutoAbs(4,4)+Aij(4,3)*AutoAbs(4,3)+Aij(4,2)*AutoAbs(4,2)),0,0,0;    %2p4
                0,0,0,0,0,0,0,(Aij(3,5)*AutoAbs(3,5)+Aij(3,4)*AutoAbs(3,4)+Aij(3,3)*AutoAbs(3,3)+Aij(3,2)*AutoAbs(3,2)),0,0;    %2p3
                0,0,0,0,0,0,0,0,(Aij(2,5)*AutoAbs(2,5)+Aij(2,4)*AutoAbs(2,4)+Aij(2,3)*AutoAbs(2,3)+Aij(2,2)*AutoAbs(2,2)),0;    %2p2
                0,0,0,0,0,0,0,0,0,(Aij(1,5)*AutoAbs(1,5)+Aij(1,4)*AutoAbs(1,4)+Aij(1,3)*AutoAbs(1,3)+Aij(1,2)*AutoAbs(1,2))];   %2p1
        
%% Depopulation vers un autre niveau 2p due à une collision avec un neutre fondamental
DepopColl2p2p=       ng*[0,0,0,0,0,0,0,0,0,0;                                                                 %2p10
                        0,(En2px_2py(9,8)+En2px_2py(9,10)),0,0,0,0,0,0,0,0;                                  %2p9
                        0,0,(En2px_2py(8,6)+En2px_2py(8,7)+En2px_2py(8,9)+En2px_2py(8,10)),0,0,0,0,0,0,0;    %2p8
                        0,0,0,(En2px_2py(7,6)+En2px_2py(7,8)+En2px_2py(7,9)),0,0,0,0,0,0;                    %2p7
                        0,0,0,0,(En2px_2py(6,7)+En2px_2py(6,8)+En2px_2py(6,9)),0,0,0,0,0;                    %2p6
                        0,0,0,0,0,(En2px_2py(5,4)+En2px_2py(5,6)+En2px_2py(5,8)),0,0,0,0;                    %2p5
                        0,0,0,0,0,0,(En2px_2py(4,3)+En2px_2py(4,5)+En2px_2py(4,6)+En2px_2py(4,7)+En2px_2py(4,8)+En2px_2py(4,9)),0,0,0;   %2p4
                        0,0,0,0,0,0,0,(En2px_2py(3,4)+En2px_2py(3,5)+En2px_2py(3,6)+En2px_2py(3,7)+En2px_2py(3,8)+En2px_2py(3,9)),0,0;   %2p3
                        0,0,0,0,0,0,0,0,(En2px_2py(2,3)),0;                                                  %2p2
                        0,0,0,0,0,0,0,0,0,0];                                                                %2p1

%% Depopulation vers un 1s due à une collision avec un neutre fondamental
DepopColl1s   =   ng*[(En2px_1s(10)),0,0,0,0,0,0,0,0,0;   %2p10
                      0,(En2px_1s(9)),0,0,0,0,0,0,0,0;   %2p9
                      0,0,(En2px_1s(8)),0,0,0,0,0,0,0;   %2p8
                      0,0,0,(En2px_1s(7)),0,0,0,0,0,0;   %2p7
                      0,0,0,0,(En2px_1s(6)),0,0,0,0,0;   %2p6
                      0,0,0,0,0,(En2px_1s(5)),0,0,0,0;   %2p5
                      0,0,0,0,0,0,(En2px_1s(4)),0,0,0;   %2p4
                      0,0,0,0,0,0,0,(En2px_1s(3)),0,0;   %2p3
                      0,0,0,0,0,0,0,0,(En2px_1s(2)),0;   %2p2
                      0,0,0,0,0,0,0,0,0,(En2px_1s(1))];  %2p1
              
%% Population due à une collision d'un 2p avec un neutre fondamental
PopColl2p=ng*[0,En2px_2py(9,10),En2px_2py(8,10),0,0,0,0,0,0,0;                                                     %2p10
              0,0,En2px_2py(8,9),En2px_2py(7,9),En2px_2py(6,9),0,En2px_2py(4,9),En2px_2py(3,9),0,0;                %2p9
              0,En2px_2py(9,8),0,En2px_2py(7,8),En2px_2py(6,8),En2px_2py(5,8),En2px_2py(4,8),En2px_2py(3,8),0,0;   %2p8
              0,0,En2px_2py(8,7),0,En2px_2py(6,7),0,En2px_2py(4,7),En2px_2py(3,7),0,0;                             %2p7
              0,0,En2px_2py(8,6),En2px_2py(7,6),0,En2px_2py(5,6),En2px_2py(4,6),En2px_2py(3,6),0,0;                %2p6
              0,0,0,0,0,0,En2px_2py(4,5),En2px_2py(3,5),0,0;                                                       %2p5
              0,0,0,0,0,En2px_2py(5,4),0,En2px_2py(3,4),0,0;                                                       %2p4
              0,0,0,0,0,0,En2px_2py(4,3),0,En2px_2py(2,3),0;                                                       %2p3
              0,0,0,0,0,0,0,0,0,0;                                                                                 %2p2
              0,0,0,0,0,0,0,0,0,0];                                                                                %2p1
end