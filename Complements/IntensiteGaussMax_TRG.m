function [I_exp,sig_I_exp,Fit] = IntensiteGaussMax_TRG(E,Lambda,Int,GraphAll,Doublet,Fit)
%% ========================================================================================
%% ======= Cette fonction calcule le FWHM à partir de la raie à 763nm et identifie ========
%% ========== l'intensité des raies à analyser en leur fittant une/des gaussienne =========
%% ========================================================================================
GraphExp=0; %Variable pour voir les Fits sur chaques raies
%% ====== INFOS SUR LA FONCTION ======
%La variable E contient les longueurs d'onde des raies potentiellement identifiables
%Les variables Lambda et Int sont les x et y du spectre expérimental
%La variable GraphExp est une variable binaire (0 ou 1) demandant ou non
    %d'afficher les graphiques des fits
%La variable Doublet est un vecteur de valeurs (0, 1 ou 2) disant si il
    %faut fiter qu'une Gaussienne, 2 ou 3 gausiennes dans l'intervalle de
    %plage_initiale = 2 nm
%------------------------------------------------------------------------------------

%% ====== DÉBUT DE LA FONCTION =======
%% Préallocation d'espace
I_exp          = zeros(1,length(E))'; 
sig_I_exp      = zeros(1,length(E))'; 
plage_initiale = 2;
%% À priori, on s'attend à trouver toutes les raies demandées
%% Premiere boucle sur les raies pour voir si certaines sont au delà du spectre mesuré
for i=1:length(E)
    if Fit(i)==1 %ne pas regarder si elle est déjà exclue
    % Si la longueur d'onde de la raie est inférieure à la plus petite longueur d'onde mesurée
    if E(i)<min(Lambda)
        I_exp(i)=0;     %et on met son intensité à 0
        Fit(i)=0; 
        E(i)
        disp('Raie hors des limites du spectre');
    % Ou si elle est plus grand que la plus grande longueur d'onde mesurée
    elseif E(i)>max(Lambda)
        I_exp(i)=0;     %et on met son intensité à 0
        Fit(i)=0;
        E(i)
        disp('Raie hors des limites du spectre')
    end
    end
end       
    
%% Ensuite on établie la valeur du bruit de fond pour fin de comparaison ultérieure
bkgd=mean(Int) ;            %Valeur d'un premier background considérant même les raies
IntBkgd=Int((Int)<bkgd);    %On exclue du spectre les raies
 %plot(IntBkgd)             %Pour voir si les raies ont bien été exclues
bkgd=mean(IntBkgd);      %Valeur d'un second background excluant les raies
   
clear IntBkgd


%% Puis on démarre la boucle pour trouver chacune de ces raies
plage = plage_initiale;
for q=1:length(E)  
 if Fit(q) == 1
        %pouvoir ploter le smooth en ayant un q=0. sinon il plot 2 fois
        %pour les q differents de zeros.
          %========= Prend juste les points autour du pic===========
        [IntTemp,LambdaTemp,moyenne] = PlageIntegration(E(q),Lambda,Int,plage,bkgd);
          %========================== fin ==========================

    
        %Si le pic détecté dans l'interval est au-dessus de 0, on continue les démarches
        if max(IntTemp)>0
            %La procédure du fit ne s'entame que si la raie se distingue du bruit de fond
            if max(IntTemp)>bkgd
                %% Puis on fit soit une gaussienne si la raie est isolée
                if Doublet(q) == 0 
                    k=0; %boucle tant que k négale pas 1
                 while k==0
                     %========= Prend juste les points autour du pic===========
                  [IntTemp,LambdaTemp,moyenne] = PlageIntegration(E(q),Lambda,Int,plage,bkgd);
                     %========================== fin ==========================
                     try %si le fit calcul un NAN, le code continu
                    [coeff] = fit(LambdaTemp',IntTemp','gauss1'); % FIT
                    [intervals] = confint(coeff); %prend l'incertitude
                     catch
                        disp('NAN compute, intensité mise à zero') 
                        I_exp(q) = 0;   %Elle est trop petite et on met son intensité à 0
                        sig_I_exp(q) = 0;
                        Fit(q) = 0;
                        break %sort de la boucle While
                     end
                    coeff = coeffvalues(coeff);  
                   
                    %Si l'amplitude de la gausienne fittée se distingue du bruit de fond
                    if coeff(1) > abs(moyenne) & coeff(1) > bkgd
                            %si la distance entre ce qu'on demande et ce quon fit est ok
                            if abs(coeff(2)-E(q))<0.4
                                %La valeur maximale de la gaussienne est l'intensité de la raie
                                I_exp(q) = coeff(1);
                                sig_I_exp(q) = sqrt((max(abs(intervals(:,1)-coeff(1))))^2 + (bkgd/coeff(1))^2  );  %prend le max de l'incertitude sur l'intensité de la raie
                                Fit(q) = 1;
                                k=1; % fin de la boucle

                            else % sinon, on baisse la plage
                               plage = plage-0.4;
                            end
                    %Si l'amplitude de la gausienne fittée ne se distingue pas du bruit de fond
                    else
                        k=1;
                        I_exp(q) = 0;   %Elle est trop petite et on met son intensité à 0
                        sig_I_exp(q) = 0;
                        Fit(q) = 0;
                         %marque l'erreur
                        disp('Intensité de raie mise à 0: Intensité < Background')
                    
                    end     
                    %====== si ca loop depuis trop longtemps=====
                  
                        if plage <0.9
                           k=1 ;
                           E(q)
                           Fit(q) =0;
                           I_exp(q) =  0;   %Et on met l'intensité de la raie à 0 
                           sig_I_exp(q) = 0;
                           %marque l'erreur
                           disp('Intensité de raie mise à 0: trouve pas la raie apres plusieurs essaies')                  
                        end
                             
                 end %fin boucle while
                    % Vérification que les gaussiennes s'ajustent bien sur chaque raie
                    if GraphExp==1
                        xfit = (coeff(2)-2):0.01:(coeff(2)+2);    %On construit un vecteur longueurs d'ondes fictif avec plus de précision que le spectre et centré sur la valeur centrale trouvé à l'aide du fit (antoine(2))
                        yfit = coeff(1).*exp(-((xfit-coeff(2))/(coeff(3))).^2);     %Produit une fonction gaussienne à l'aide des paramètres déduits du fit

                        figure
                        plot(xfit,yfit)
                        hold on
                        plot(LambdaTemp,IntTemp,'r')
                        hold off
                        clear xfit yfit
                    end
                    
                    clear coeff Ampli_fit FWHM_fit ecart IntTemp LambdaTemp
                    
                end %fin du fit sur une gausienne
                  clear k
                 
                 
                %% Soit deux gausiennes s'il s'agit d'une raie d'un doublet
                plage = plage_initiale;
                if Doublet(q) == 1 %si on a indiqué que c'est un doublet
                    k=0; %loop tant que ca ne vaut pas 1
                    plage = plage_initiale+0.8; %plus grande plage pour un doublet
                 while k==0 %bouléen, arrêter si on a la bonne raie.
                     
                      %========= Prend juste les points autour du pic===========
                    [IntTemp,LambdaTemp,moyenne] = PlageIntegration(E(q),Lambda,Int,plage,bkgd);
                      %========================== fin ==========================
                    
                      %Fit 2 gaussiennes
                    [coeff] = fit(LambdaTemp',IntTemp','gauss2');
                    [intervals] = confint(coeff); %prend l'incertitude
                    coeff = coeffvalues(coeff);

                    %Extraction des coefficients des fits pour les tests conditionnels (FWHM et Amplitude)
                    Ampli_fit(1) = coeff(1); Ampli_fit(2) = coeff(4);
                    FWHM_fit(1)=sqrt(8)*(log(2))^(3/2)*coeff(6); FWHM_fit(2)=sqrt(8)*(log(2))^(3/2)*coeff(3);
                     
                    %Extraction de l'eccart entre le fit et le lambda demandé
                    ecart(1)=abs(E(q)-coeff(2));
                    ecart(2)=abs(E(q)-coeff(5));
                   
                    %===on check que la raie la plus proche est au dessus du bruit===  
                      [mini pmin] = min(ecart);%pointe le minimum de ecart
                      if Ampli_fit(pmin) > abs(moyenne) & Ampli_fit(pmin) > bkgd
                        %Si les deux fit ont des FWHM qui ne sont pas
                        %démesurément différents (car en principe fixés par le spectro) (ATTENTION CAR AVEC ISOPLANE LES FWHH DES DOUBLETS SONT QUAND MËME DISPARATES) 
                       
                        if max(FWHM_fit)<5*min(FWHM_fit)     
                            %Et si les longueurs d'onde centrales sont assez distantes pour bien correspondre a deux raies
                            if abs(coeff(2)-coeff(5))>0.1
                               %On vérifie si la longueur d'onde recherchée est plus près de la longueur d'onde centrale fittée de la 1ere ou la 2eme Gaussienne                          
                                if ecart(2) < ecart(1) 
                                    I_exp(q) = coeff(4);
                                    sig_I_exp(q) = sqrt((max(abs(intervals(:,4)-coeff(4))))^2 + (bkgd/coeff(4))^2);
                                    Fit(q) = 1;
                                else
                                    I_exp(q) = coeff(1);
                                    sig_I_exp(q) = sqrt((max(abs(intervals(:,1)-coeff(1))))^2 + (bkgd/coeff(1))^2);
                                    Fit(q) = 1;
                                end      
                            %Par contre si les longueurs d'onde centrales sont trop proches, le fit echoue    
                            else
                                I_exp(q) = 0;   %Et on met l'intensité de la raie à 0
                                sig_I_exp(q) = 0;
                                Fit(q) = 0;    %Ainsi que Fit à 0 car on n'a pas réussi à la trouver
                                %écrit l'erreur
                                disp('Intensité de raie mise à 0: échec du fit, longueurs d onde centrales trop proches')   
                            end
                        %Par contre si le FWHM d'une des raies est démesurément tros large p/r à l'autre, le fit échoue
                        else
                            Fit(q) = 0;
                            I_exp(q) = 0;   %Et on met l'intensité de la raie à 0 
                            sig_I_exp(q) = 0;
                            %écrit l'erreur
                            disp('Intensité de raie mise à 0: échec du fit, FWHM trop différents')
                        end
                        
                      else
                        I_exp(q) = 0;   %On met l'intensité de la raie à 0
                        sig_I_exp(q) = 0;
                        Fit(q) = 0;
                        %écrit l'erreur
                        disp('Intensité de raie mise à 0: échec du fit, amplitude trop petite')
                      end
                    
                    %regarde si l'eccart entre la position du fit et la raie est trop
                    %grand, si oui, ont prend moins de points.
                    if min(ecart) < 0.5
                        k=1;     %ont ne refait pas une boucle 
                    else
                       plage = plage-0.5;
                    end
                   %===arret de la boucle si ca fait trop longtemps===
                   if plage <1.5
                      k=1 ;
                      Fit(q) = 0;
                      I_exp(q) = 0;   %Et on met l'intensité de la raie à 0 
                      sig_I_exp(q) = 0;
                      %écrit l'erreur
                      disp('Intensité de raie mise à 0: doublet introuvable')                  
                   end

                 end %fin de la boucle while
                  % Vérification que les gaussiennes s'ajustent bien sur chaque raie
                   if GraphExp==1
                      xfit = (coeff(5)-3):0.01:(coeff(5)+5);    %On construit un vecteur longueurs d'ondes fictif avec plus de précision que le spectre et centré sur la valeur centrale trouvé à l'aide du fit (antoine(2))
                      yfit = coeff(1).*exp(-((xfit-coeff(2))/(coeff(3))).^2)+coeff(4).*exp(-((xfit-coeff(5))/(coeff(6))).^2);     %Produit une fonction gaussienne à l'aide des paramètres déduits du fit

                      figure
                      plot(xfit,yfit)
                      hold on
                      plot(LambdaTemp,IntTemp,'r')
                      hold off
                      clear xfit yfit
                   end
                   
                  clear coeff Ampli_fit FWHM_fit ecart IntTemp LambdaTemp
                 
                end %fin boucle doublet=1  
                k=0;
                 clear k
                              
                 %% Soit trois gausiennes s'il s'agit d'une raie d'un triplet
                  plage = plage_initiale;   
                if Doublet(q) == 2
                    k=0; %Bouléen, arrete la boucle si condition remplis
                    plage = plage_initiale+0.8; %plus grande plage pour un triplet!
                while k==0
                         %========= Prend juste les points autour du pic===========
                        [IntTemp,LambdaTemp,moyenne] = PlageIntegration(E(q),Lambda,Int,plage,bkgd);
                         %========================== fin ==========================
                    
                    [coeff rmse] = fit(LambdaTemp',IntTemp','gauss3');
                    [intervals] = confint(coeff); %prend l'incertitude
                    coeff = coeffvalues(coeff);
                    %Extraction des coefficients des fits pour les tests conditionnels (FWHM et Amplitude)
                    Ampli_fit(1) = coeff(1); Ampli_fit(2) = coeff(4); Ampli_fit(3) = coeff(7);
                    FWHM_fit(1)=sqrt(8)*(log(2))^(3/2)*coeff(6); FWHM_fit(2)=sqrt(8)*(log(2))^(3/2)*coeff(3); FWHM_fit(3)=sqrt(8)*(log(2))^(3/2)*coeff(9);
                    %Extraction de l'eccart entre le fit et le lambda demandé
                                ecart(1)=abs(E(q)-coeff(2));
                                ecart(2)=abs(E(q)-coeff(5));
                                ecart(3)=abs(E(q)-coeff(8));                    
                    
                    %Si l'amplitude de la plus petite gaussienne se distingue du bruit de fond
                    %===on check que la raie la plus proche est au dessus du bruit===  
                     [mini pmin] = min(ecart);%pointe le minimum de ecart
                     if Ampli_fit(pmin) > abs(moyenne) & Ampli_fit(pmin) > bkgd
                        %Si les deux fit ont des FWHM qui ne sont pas
                        %démesurément différents (car en principe fixés par le spectro) (ATTENTION CAR AVEC ISOPLANE LES FWHH DES DOUBLETS SONT QUAND MËME DISPARATES) ++++++++++++++++++++++++++++++++++++++++++
                        if (FWHM_fit(pmin)==min(FWHM_fit) & max(FWHM_fit)>5*FWHM_fit(pmin)) | FWHM_fit(pmin)==max(FWHM_fit)  &  FWHM_fit(pmin)>5*min(FWHM_fit)        
                            I_exp(q) = 0;   %Et on met l'intensité de la raie à 0
                            sig_I_exp(q) = 0;
                            Fit(q) = 0;
                            %écrit l'erreur
                            disp('Intensité de raie mise à 0: échec du fit, FWHM trop différents')              
                        
                        %Par contre si le FWHM d'une des raies est démesurément tros large p/r à l'autre, le fit échoue
                        else   
                               %check la gaussienne la plus prete
                                if find(ecart == min(ecart)) ==1
                                    I_exp(q) = coeff(1);
                                    sig_I_exp(q) = sqrt((max(abs(intervals(:,1)-coeff(1))))^2+ (bkgd/coeff(1))^2);
                                    Fit(q) = 1;
                                elseif find(ecart == min(ecart)) ==2
                                    I_exp(q) = coeff(4);
                                    sig_I_exp(q) = sqrt((max(abs(intervals(:,4)-coeff(4))))^2+ (bkgd/coeff(4))^2);
                                    Fit(q) = 1;
                                else
                                    I_exp(q) = coeff(7); 
                                    sig_I_exp(q) = sqrt((max(abs(intervals(:,7)-coeff(7))))^2 + (bkgd/coeff(7))^2);
                                    Fit(q) = 1;
                                end   
                        end
                        
                    %Si l'amplitude de la plus petite gaussienne est trop petite
                     else
                        I_exp(q) = 0;   %On met l'intensité de la raie à 0
                        sig_I_exp(q) = 0;
                        Fit(q) = 0;
                         %écrit l'erreur
                         disp('Intensité de raie mise à 0: échec du fit, amplitude trop petite')   
                    end
                    
                    %regarde si l'eccart entre le fit et la raie est trop
                    %grand, si oui, ont prend moins de points.
                    
                    if min(ecart) < 0.4
                        k=1; %ont ne refait pas une boucle
                    else
                       plage = plage-0.5;
                    end
                      %===arret de la boucle si ca fait trop longtemps===
                     if plage <1.5
                      k=1 ;
                      Fit(q) = 0;
                      I_exp(q) = 0;   %Et on met l'intensité de la raie à 0 
                      sig_I_exp(q) = 0;
                      %écrit l'erreur
                      disp('Intensité de raie mise à 0: doublet introuvable')   
                   end
                  
                    
                 end %fin de la boucle while
                           % Vérification que les gaussiennes s'ajustent bien sur chaque raie
                           if GraphExp==1
                               xfit = (coeff(5)-3):0.01:(coeff(5)+5);    %On construit un vecteur longueurs d'ondes fictif avec plus de précision que le spectre et centré sur la valeur centrale trouvé à l'aide du fit (antoine(2))
                               %yfit = coeff(1).*exp(-((xfit-coeff(2))/(coeff(3))).^2)+coeff(4).*exp(-((xfit-coeff(5))/(coeff(6))).^2)+coeff(7).*exp(-((xfit-coeff(8))/(coeff(9))).^2);     %Produit une fonction gaussienne à l'aide des paramètres déduits du fit
                               yfit1 = coeff(1).*exp(-((xfit-coeff(2))/(coeff(3))).^2);
                                    yfit2 = coeff(4).*exp(-((xfit-coeff(5))/(coeff(6))).^2);
                                    yfit3 = coeff(7).*exp(-((xfit-coeff(8))/(coeff(9))).^2);
                               figure
                               plot(xfit,yfit1)
                               hold on
                               plot(xfit,yfit2)
                               plot(xfit,yfit3)
                               plot(LambdaTemp,IntTemp,'r')
                               hold off
                               clear xfit yfit
                           end      
                           
                   clear coeff Ampli_fit FWHM_fit ecart IntTemp LambdaTemp
                 
                end %fin boucle doublet=1  
                 k=0;
             clear k;
             plage = plage_initiale;
            %Si par contre le max dans l'interval regardé se distingue trop peu du bruit de fond
            else
                I_exp(q)=0; %On met automatiquement son intensité à 0
                sig_I_exp(q) = 0;
                Fit(q) = 0;
                %écrit l'erreur
                disp('Intensité de raie mise à 0: plus petite que le Background')            
         end
            
        %Si le pic est sous 0, on abandonne tout de suite    
        else 
            I_exp(q)=0;         %On met automatiquement son intensité à 0
            sig_I_exp(q) = 0;
            Fit(q) = 0;
            %écrit l'erreur
            disp('Intensité de raie mise à 0: intensité négative')  
        end
        clear IntTemp LambdaTemp
 end  %end sur le if fit(i)=0   
end

%% Visualisation de l'obtention des intensités fittées pour comparaison 
if GraphExp==1
    figure
    set(gcf,'color','w');
    plot(Lambda,(Int-bkgd))
    hold on
    scatter(E,I_exp)
    hold off
end

if GraphAll==1
    figure
    set(gcf,'color','w');
    plot(Lambda,(Int),'k')
    hold on
    %met les zeros en NAN donc pas affichés
    temp_I_exp=I_exp;
    temp_I_exp(temp_I_exp==0) = nan;
    scatter(E(1:2),temp_I_exp(1:2),'MarkerEdgeColor',[1 0.5 0.2])
    scatter(E(3:21),temp_I_exp(3:21),'MarkerEdgeColor' ,[0 0.8 0])
    scatter(E(22:33),temp_I_exp(22:33),'r')
    scatter(E(34:41),temp_I_exp(34:41),'MarkerEdgeColor' ,[0 1 1])
    hold off
end
clear q
end