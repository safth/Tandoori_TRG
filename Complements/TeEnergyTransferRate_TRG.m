%% ==============================================================================
%% =======Cette fonction sert à mettre en matrice les taux de transfert =========
%% ====================d'énergie entre niveaux 2px et 2py =======================
%% ==============================================================================

En2px_2py=zeros(10);                %En2px_2py(2,3) = 2p2 ->2p3

% ============== Energy transfer atom atom rate coeff. (atom collision) from PRA,18,(1978),1388===============
En2px_2py(2,3)=0.5e-12*1e-6;
En2px_2py(3,4)=27.0e-12*1e-6; %22
En2px_2py(3,5)=0.3e-12*1e-6;
En2px_2py(3,6)=44e-12*1e-6; %41
En2px_2py(3,7)=1.4e-12*1e-6;
En2px_2py(3,8)=1.9e-12*1e-6;
En2px_2py(3,9)=0.8e-12*1e-6;
En2px_2py(4,5)=0.7e-12*1e-6;
En2px_2py(4,6)=4.8e-12*1e-6;
En2px_2py(4,7)=3.2e-12*1e-6;
En2px_2py(4,8)=1.4e-12*1e-6;
En2px_2py(4,9)=3.3e-12*1e-6;
En2px_2py(5,6)=11.3e-12*1e-6;
En2px_2py(5,8)=9.5e-12*1e-6;
En2px_2py(6,7)=4.1e-12*1e-6; %8
En2px_2py(6,8)=6e-12*1e-6;  %8
En2px_2py(6,9)=1e-12*1e-6; %2
En2px_2py(7,8)=14.3e-12*1e-6; %16
En2px_2py(7,9)=23.3e-12*1e-6;
En2px_2py(8,9)=18.2e-12*1e-6;  %17
En2px_2py(8,10)=1e-12*1e-6;
En2px_2py(9,10)=5.1e-12*1e-6;  %18

% ============== Energy transfer atom atom rate coeff. reported PRA,18,(1978),1388 using detailed balance ===============
En2px_2py(4,3)=23e-12*1e-6; %1.17
En2px_2py(5,4)=1.7e-12*1e-6;
En2px_2py(7,6)=2.5e-12*1e-6;
En2px_2py(8,6)=0.3e-12*1e-6;
En2px_2py(8,7)=0.8e-12*1e-6;
En2px_2py(9,8)=6.8e-12*1e-6;

% =========== Additional Energy transfer atom atom rate coeff. obtained by intrapolating the rates reported in from PRA,18,(1978),1388 ==============
En2px_2py(1,2)=3.22679*1e-12*1e-6;
En2px_2py(1,3)=2.18211*1e-12*1e-6;
En2px_2py(1,4)=1.63954*1e-12*1e-6;
En2px_2py(1,5)=1.41052*1e-12*1e-6;
En2px_2py(1,6)=3.08606*1e-13*1e-6;
En2px_2py(1,7)=2.31874*1e-13*1e-6;
En2px_2py(1,8)=9.68849*1e-14*1e-6;
En2px_2py(1,9)=7.27951*1e-14*1e-6;
En2px_2py(1,10)=5.72518*1e-15*1e-6;
En2px_2py(2,4)=1.61418*1e-11*1e-6;
En2px_2py(2,5)=1.3887*1e-11*1e-6;
En2px_2py(2,6)=3.03832*1e-12*1e-6;
En2px_2py(2,7)=2.28286*1e-12*1e-6;
En2px_2py(2,8)=9.5386*1e-13*1e-6;
En2px_2py(2,9)=7.1669*1e-13*1e-6;
En2px_2py(2,10)=5.6366*1e-14*1e-6;
En2px_2py(3,10)=8.33512*1e-14*1e-6;
En2px_2py(4,10)=1.10934*1e-13*1e-6;
En2px_2py(5,7)=5.22242*1e-12*1e-6;
En2px_2py(5,9)=1.63954*1e-12*1e-6;
En2px_2py(5,10)=1.28947*1e-13*1e-6;
En2px_2py(6,10)=5.89365*1e-13*1e-6;
En2px_2py(7,10)=7.844*1e-13*1e-6;

save('En2px_2py.mat','En2px_2py')