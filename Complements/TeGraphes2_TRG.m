function TeGraphes2_TRG(Ne,Te,SEGD,a,b,c,I_exp,sig_I_exp,FitCorrige,ChoixErreur,moyenne,P,E,ChoixAutoabs)

% %% =====================================================================================================
% % ======== Cette fonction trace le graphiques 3D de SD en fonction de nm et Te ainsi que ceux de =======
% % ==== divers quantité (Densité des niveaux 2p, % des contributions aux gains/pertes de population, ====
% % ============ intensité et auto-absorption des raies) pour les paramètre Te et nm optimaux ============
% % ======================================================================================================

%% Préallocation d'espace
 w1=zeros(1,41); % les poid de la weighted STD
 w2=zeros(41,2);
% Population=zeros(10,4);
% Depopulation=zeros(10,3);
% PopFond=zeros(10,2);

%% Extraction des données des matrices
I_theo=zeros(1,size(a,4));
sig_I_theo=zeros(1,size(b,4));
Thetaij=zeros(1,size(c,4));
I_theo(:)=a(1,1,1,:);
sig_I_theo(:)=b(1,1,1,:);
Thetaij(:)=c(1,1,1,:);
AutoAbsTheo=1-Thetaij;

  %rejete les raies avec une incertitude plus grande que la
  %valeur de l'intensité (donne des poids < 1 ...caca) 
if ChoixErreur ==1  
    for m=1:length(I_theo)
        if 1- sqrt((sig_I_theo(m)/I_theo(m))^2 + (sig_I_exp(m)/I_exp(m))^2)<0
           I_theo(m) = 0;
           sig_I_theo(m) = 0;
           I_exp(m)=0;
           sig_I_exp(m)=0;
           E(m)
           disp('intensité plus petite que l incertitude')
        end
    end
end


clear k gaz m
%% Normalisation des intensités
I_theo=I_theo.*FitCorrige;
sig_I_theo=sig_I_theo.*FitCorrige;

%       Ratio = (I_exp')./(I_theo');
%       Ratio(isnan(Ratio)) = 0 ;
%       fopen( 'Ratio_Tandoori.txt','wt')
%       for i=1: length(Ratio)
%      fileID = fopen( 'Ratio_Tandoori.txt','at'); 
%       formatSpec = '%e \r\n';
%       fprintf(fileID,formatSpec,Ratio(i));
%       fclose(fileID);
%       end 
%       clear Ratio
% plot
Ratio = (I_exp')./(I_theo');
Ratio(isnan(Ratio)) = 0 ;   
Dat=(Ratio/mean(Ratio(Ratio~=0)));
Dat(Dat==1) =nan;

clear tempData

        I_theo=I_theo*moyenne;
    sig_I_theo=sig_I_theo*moyenne;


% ############# fig temporaire de I exp/I theo

    figure2=figure;
    set(figure2,'name','Rapport Iexp/Itheo','numbertitle','off','Position',[200 600 1000 600])
    %subplot('position',[0.05 0.55 0.90 0.36])    
    bar3=bar(1:length(Ratio),Dat , 0.75);
    %bar3=bar(1:length(Ratio),Ratio/mean(Ratio(Ratio~=0)) , 0.75);
    set(gca,'XTick',1:41)
    linename=['585';'640'; '667'; '696'; '706'; '714'; '727'; '738'; '750'; '751'; '763'; '794'; '800'; '801'; '810'; '811'; '826'; '840'; '842'; '852'; '866'; '758'; '760'; '768'; '769'; '785'; '810'; '811'; '819'; '829'; '850'; '877'; '892'; '788'; '820'; '823'; '828'; '834'; '881'; '895'; '904'];
    title('Ratio','fontweight','bold')
    legend('Iexp/Icalc','Location','northwest');
    set(gca,'XTickLabel',linename)


%##########################################################################
%##########################################################################    
%########################## Figure 1 ######################################    
%##########################################################################    
%##########################################################################  
 

%% Graphe de comparaison des intensités des raies
    %dimension et nom de la figure
    figure2=figure;
    set(figure2,'name','Intensité de raie et Auto absorption','numbertitle','off','Position',[200 600 1000 800])
    if ChoixAutoabs==1
        subplot(3,1,1)
       % subplot('position',[0.05 0.55 0.90 0.36])    
    else
        subplot(2,1,1)
        subplot('position',[0.05 0.55 0.90 0.36])    
    end
    bar3=bar(1:length(I_exp), [I_theo' I_exp'], 0.75, 'grouped');
    set(bar3(1),'FaceColor',[1 0 0])
    set(bar3(2),'FaceColor',[0 0 1])
    set(gca,'XTick',1:41)
    linename=['585';'640'; '667'; '696'; '706'; '714'; '727'; '738'; '750'; '751'; '763'; '794'; '800'; '801'; '810'; '811'; '826'; '840'; '842'; '852'; '866'; '758'; '760'; '768'; '769'; '785'; '810'; '811'; '819'; '829'; '850'; '877'; '892'; '788'; '820'; '823'; '828'; '834'; '881'; '895'; '904'];
    title('Comparaison intensité prédite vs mesurée','fontweight','bold')
    legend('Intensité théorique','Intensité expérimentale','Location','northwest');
    set(gca,'XTickLabel',linename)
    
     
    %pogne les max pour afficher les "2px" sur le graphique 
    %qui compare Ith et Iobs
    for i=1:41
        if I_theo(i)>=I_exp(i)
        position(i)=I_theo(i);
        elseif I_theo(i)<I_exp(i)
        position(i)=I_exp(i);
        else
        position(i)=200;
        end
    end
    %ajout du nom des niveaux
    name = {'2P1' '2P9'... %Ne
            '2P1' '2P2' '2P3' '2P4' '2P2' '2P3' '2P1' '2P5' '2P6' '2P4' '2P6' '2P8' '2P7' '2P9' '2P2' '2P3' '2P8' '2P4' '2P7' ... %Ar
            '2P5' '2P6' '2P1' '2P7' '2P3' '2P8' '2P9' '2P6' '2P7' '2P4' '2P8' '2P10' ... %Kr
            '2P1' '2P4' '2P6' '2P5' '2P3' '2P8' '2P6' '2P9'}; %Xe
    text(1:41,position,name,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

      set(gca, 'yscale', 'log')
    
      
%% Graphe de l'auto-absorption sur chaque raie analysée
if ChoixAutoabs==1
    subplot(3,1,2)
   % subplot('position',[0.05 0.07 0.90 0.36])
    bar(1:length(I_exp), 100*AutoAbsTheo, 0.6, 'grouped')
    ylabel('%','FontSize',11);
    set(gca,'XTick',1:41)
    linename=['585';'640'; '667'; '696'; '706'; '714'; '727'; '738'; '750'; '751'; '763'; '794'; '800'; '801'; '810'; '811'; '826'; '840'; '842'; '852'; '866'; '758'; '760'; '768'; '769'; '785'; '810'; '811'; '819'; '829'; '850'; '877'; '892'; '788'; '820'; '823'; '828'; '834'; '881'; '895'; '904'];
    title('Auto-absorption des raies','fontweight','bold');
    set(gca,'XTickLabel',linename)
end

%Graphe des Weights
if ChoixAutoabs == 1 
   subplot(3,1,3)
else
    subplot(2,1,2)
    subplot('position',[0.05 0.07 0.90 0.36])
end
if ChoixErreur == 1
    w1 = 1-sqrt((sig_I_theo./I_theo).^2 + (sig_I_exp./I_exp).^2) ; %Poid pour la STD
       w2(:,1)= w1.*(sig_I_theo./I_theo)./((sig_I_exp./I_exp+sig_I_theo./I_theo));
    w2(:,2)= w1.*(sig_I_exp./I_exp)./((sig_I_exp./I_exp+sig_I_theo./I_theo));
    

    bar11=bar(1:length(I_exp),w2,0.70,'stacked','k');
    ylabel('Poid statistique','FontSize',11);
    set(gca,'XTick',1:41)
    ylim([0 1]);
    linename=['585';'640'; '667'; '696'; '706'; '714'; '727'; '738'; '750'; '751'; '763'; '794'; '800'; '801'; '810'; '811'; '826'; '840'; '842'; '852'; '866'; '758'; '760'; '768'; '769'; '785'; '810'; '811'; '819'; '829'; '850'; '877'; '892'; '788'; '820'; '823'; '828'; '834'; '881'; '895'; '904'];
    title('Poids statistiques','fontweight','bold');
    set(gca,'XTickLabel',linename)
     hold on
    bar(1:length(I_exp), w1, 0.35, 'Grouped')
    set(bar11(1),'FaceColor',[1 1 0])
    set(bar11(2),'FaceColor',[0 1 0])
        lgnd=legend('% sigma Théo','% sigma Exp','poids','Orientation','horizontal','Location','south');
    hold off   
elseif ChoixErreur == 2
    w1 = 1 ./ sqrt((sig_I_theo./I_theo).^2 + (sig_I_exp./I_exp).^2) ; %Poid pour la STD
    w2(:,1)= w1.*(sig_I_theo./I_theo)./((sig_I_exp./I_exp+sig_I_theo./I_theo));
    w2(:,2)= w1.*(sig_I_exp./I_exp)./((sig_I_exp./I_exp+sig_I_theo./I_theo));
    

    bar11=bar(1:length(I_exp),w2,0.70,'stacked','k');
    ylabel('Poid statistique','FontSize',11);
    set(gca,'XTick',1:41)
    if min(w1)<0
        ylim([0.5*min(w1) 1.1*max(w1)]);
    else
        ylim([0 1.1*max(w1)]);
    end
    linename=['585';'640'; '667'; '696'; '706'; '714'; '727'; '738'; '750'; '751'; '763'; '794'; '800'; '801'; '810'; '811'; '826'; '840'; '842'; '852'; '866'; '758'; '760'; '768'; '769'; '785'; '810'; '811'; '819'; '829'; '850'; '877'; '892'; '788'; '820'; '823'; '828'; '834'; '881'; '895'; '904'];
    title('Poids statistiques','fontweight','bold');
    set(gca,'XTickLabel',linename)
     hold on
    bar(1:length(I_exp), w1, 0.35, 'Grouped')
    set(bar11(1),'FaceColor',[1 1 0])
    set(bar11(2),'FaceColor',[0 1 0])
        lgnd=legend('% sigma Théo','% sigma Exp','poids','Orientation','horizontal','Location','south');
    hold off
elseif ChoixErreur == 3
    w1 = sig_I_theo./sig_I_theo; %donc 1
     bar(1:length(I_exp), w1, 0.35, 'Grouped')
         ylabel('Poid statistique','FontSize',11);
    set(gca,'XTick',1:41)
    linename=['585';'640'; '667'; '696'; '706'; '714'; '727'; '738'; '750'; '751'; '763'; '794'; '800'; '801'; '810'; '811'; '826'; '840'; '842'; '852'; '866'; '758'; '760'; '768'; '769'; '785'; '810'; '811'; '819'; '829'; '850'; '877'; '892'; '788'; '820'; '823'; '828'; '834'; '881'; '895'; '904'];
    title('Poids statistiques','fontweight','bold');
    set(gca,'XTickLabel',linename)
        lgnd=legend('poids','Orientation','horizontal','Location','south');

end

    

%##########################################################################
%##########################################################################    
%########################## Figure 2 ######################################    
%##########################################################################    
%########################################################################## 
 %% Graphe 3D de SD en fonction de nm et Te
     figure
     surf(Ne,Te,SEGD','edgecolor','none') 
     colormap(hsv)
     set(gca, 'xscale', 'log')
     %set(gca, 'zscale', 'log')
     xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
     if ChoixErreur==1
         zlabel({'Squared Eucledian Distance'},'FontSize',10,'fontweight','bold')
         title('Squared Eucledian Distance','fontweight','bold')
     elseif ChoixErreur==2
         zlabel({'Percentage standard error'},'FontSize',10,'fontweight','bold')
         title('Percentage standard error','fontweight','bold')
     end
   