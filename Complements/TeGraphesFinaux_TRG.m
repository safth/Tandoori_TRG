function TeGraphesFinaux(NumFichier,NmOptimal,NmMin,NmMax,TeOptimal,TeMin,TeMax,Iexp)

%% Évolution de l'intensité normalisée de certaines raies
% figure
% %% 585 Ne
%     subplot(2,4,1)
%     subplot('position',[0.04 0.55 0.19 0.34])
%     scatter(NumFichier,Iexp(:,1),'filled')
%     title('585 Ne')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
% 
% %% 640 Ne
%     subplot(2,4,2)
%     subplot('position',[0.29 0.55 0.19 0.34])
%     scatter(NumFichier,Iexp(:,2),'filled')
%     title('640 Ne')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
%         
% %% 750 Ar
%     subplot(2,4,3)
%     subplot('position',[0.55 0.55 0.19 0.34])
%     scatter(NumFichier,Iexp(:,9),'filled')
%     title('750 Ar')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
%     
% %% 751 Ar
%     subplot(2,4,4)
%     subplot('position',[0.79 0.55 0.19 0.34])
%     scatter(NumFichier,Iexp(:,10),'filled')
%     title('751 Ar')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
%     
% %% 826 Ar
%     subplot(2,4,5)
%     subplot('position',[0.04 0.07 0.19 0.34])
%     scatter(NumFichier,Iexp(:,17),'filled')
%     title('801 Ar')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
%     
% %% 760 Kr
%     subplot(2,4,6)
%     subplot('position',[0.29 0.07 0.19 0.34])
%     scatter(NumFichier,Iexp(:,23),'filled')
%     title('760 Kr')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
%     
% %% 819 Kr
%     subplot(2,4,7)
%     subplot('position',[0.55 0.07 0.19 0.34])
%     scatter(NumFichier,Iexp(:,29),'filled')
%     title('819 Kr')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')
%     ylim([0 1]);
%     
% %% 828 Xe
%     subplot(2,4,8)
%     subplot('position',[0.79 0.07 0.19 0.34])
%     scatter(NumFichier,Iexp(:,37),'filled')
%     title('828 Xe')
%     xlabel('#Fichier') 
%     ylabel('Intensité normalisée')        
%     ylim([0 1]);
%     
%     
figure

%% Évolution des métastables en fonction du fichier
subplot(2,1,1)
subplot('position',[0.04 0.55 0.9 0.36])
scatter(NumFichier,NmOptimal,'filled')
set(gca, 'yscale', 'log')
title('Densite électronique')
xlabel('#Fichier') 
ylabel('Densité électronique')
hold on
scatter(NumFichier,NmMin,2,'r')
scatter(NumFichier,NmMax,2,'r')
hold off

%% Évolution de la température en fonction du fichier
subplot(2,1,2)
subplot('position',[0.04 0.07 0.9 0.36])
scatter(NumFichier,TeOptimal,'filled')
title('Temperature electronique')
xlabel('#Fichier')
ylabel('Temprature electronique')
hold on
scatter(NumFichier,TeMin,2,'r')
scatter(NumFichier,TeMax,2,'r')
hold off