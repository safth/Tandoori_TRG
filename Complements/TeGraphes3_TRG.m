function TeGraphes3_TRG(densite2p,Gains2p,Pertes2p,ContributionFond,energie1s,energie2p,densite1s,Gains1s,Pertes1s,ng)
% %% =====================================================================================================
% % ======== Cette fonction trace le graphiques 3D de SD en fonction de nm et Te ainsi que ceux de =======
% % ==== divers quantité (Densité des niveaux 2p, % des contributions aux gains/pertes de population, ====
% % ============ intensité et auto-absorption des raies) pour les paramètre Te et nm optimaux ============
% % ======================================================================================================

%% Préallocation d'espace
 I_theo=zeros(1,36);
 AutoAbsTheo=zeros(1,36);





global gaz_i gaz_f
%% Extraction des données des matrices
for gaz=gaz_i:gaz_f 
    %%
    for l=1:10
   Population2p(gaz,l,1)=Gains2p(gaz,1,1,l,1); %Population impact électronique sur le fondamental
   Population2p(gaz,l,2)=Gains2p(gaz,1,1,l,2); %Population impact électronique sur les 1s
   Population2p(gaz,l,3)=Gains2p(gaz,1,1,l,4); %Population ré-absorption de photon
   Population2p(gaz,l,4)=Gains2p(gaz,1,1,l,5); %Transfert radiatif du Ar(1s5)
   Depopulation2p(gaz,l,1)=Pertes2p(gaz,1,1,l,1);     %Dépopulation desexcitation radiative
   Depopulation2p(gaz,l,2)=Pertes2p(gaz,1,1,l,2);%/(Pertes(1,l)+Pertes(2,l)+Pertes(3,l));      %Dépopulation transfert collisionnel entre les 2p
   Depopulation2p(gaz,l,3)=Pertes2p(gaz,1,1,l,3);%/(Pertes(1,l)+Pertes(2,l)+Pertes(3,l));      %Dépopulation transfert collisionnel vers les 1s
   
   PopFond2p(gaz,l,1)=ContributionFond(gaz,1,1,l);
   PopFond2p(gaz,l,2)=100-ContributionFond(gaz,1,1,l);
    end
    for l=1:5
   Population1s(gaz,l,1)=Gains1s(gaz,1,1,l,1);
   Population1s(gaz,l,2)=Gains1s(gaz,1,1,l,2);
   Population1s(gaz,l,3)=Gains1s(gaz,1,1,l,3)+Gains1s(gaz,1,1,l,4); %mix pis updown c'est le meme mécanisme
   Population1s(gaz,l,4)=Gains1s(gaz,1,1,l,5);
   Depopulation1s(gaz,l,1)=Pertes1s(gaz,1,1,l,1);
   Depopulation1s(gaz,l,2)=Pertes1s(gaz,1,1,l,2);
   Depopulation1s(gaz,l,3)=Pertes1s(gaz,1,1,l,3);
   Depopulation1s(gaz,l,4)=Pertes1s(gaz,1,1,l,4);
   Depopulation1s(gaz,l,5)=Pertes1s(gaz,1,1,l,5);
   Depopulation1s(gaz,l,6)=Pertes1s(gaz,1,1,l,6);
    end   
    
   
   
end



%##########################################################################
%##########################################################################    
%########################## Figure 1 ######################################    
%##########################################################################    
%##########################################################################  
for gaz=gaz_i:gaz_f 
    figure1=figure;
    name = {'Pop Depop He' 'Pop Depop Ne' 'Pop Depop Ar' 'Pop Depop Kr' 'Pop Depop Xe'};
    set(figure1,'name',name{gaz},'numbertitle','off','Position',[200 600 1000 600])
%% Graphe de la densité des niveaux 2p

      subplot(2,3,4)
      subplot('position',[0.04 0.07 0.29 0.36])
      scatter(energie2p(gaz,:)',densite2p(gaz,1,1,:)/ng(gaz),'filled','d')
      set(gca, 'yscale', 'log')
      xlabel('Énergies Niveaux 2p_x','FontSize',11,'fontweight','bold');
      title('Densite des niveaux 2p_x par rapport au Fondamental','fontweight','bold')

%% Graphe des mécanismes de population 2P
    pop1(:,:) = Population2p(gaz,:,:);
   % pop2(:,:) = PopFond2p(gaz,:,:);
    
    subplot(2,3,5)
    subplot('position',[0.38 0.07 0.29 0.36])    
    bar1=bar(1:10,pop1,0.75,'stacked');
    set(bar1(4),'FaceColor',[1 0 0],'LineStyle','none')
    set(bar1(2),'FaceColor',[0 1 0],'LineStyle','none')
    set(bar1(3),'FaceColor',[0 0 1],'LineStyle','none')
    set(bar1(1),'FaceColor',[1 0 1],'LineStyle','none')
    set(gca,'XTick',1:10,'Xdir','reverse')
    ylim([0 120]);
    xlim([0.4 10.6]);
    title('Contribution des mécanismes de population des 2p','fontweight','bold')
    xlabel('Niveau 2p_x','FontSize',11,'fontweight','bold');
    lgnd=legend('Coll e^- - fond','Coll e^- - 1s','Auto-abs','Transfert Ar(1s_5)','Orientation','horizontal','Location','north');
    set(lgnd,'FontSize',8)
%     hold on
%     bar11=bar(1:10,pop2,0.35,'stacked');
%     set(bar11(1),'FaceColor',[1 0 1])
%     set(bar11(2),'FaceColor',[0 1 0])
%     set(gca,'XTick',1:10,'Xdir','reverse')
%     ylim([0 120]);
%     xlim([0.4 10.6]);
%     hold off
     clear pop1 

    
    
%% Graphe des mécanismes de dépopulation des 2p
    depop(:,:)=Depopulation2p(gaz,:,:);
    subplot(2,3,6)
    subplot('position',[0.70 0.07 0.29 0.36])
    bar2=bar(1:10,depop,0.75,'stacked');
    set(bar2(1),'FaceColor',[1 0 0],'LineStyle','none')
    set(bar2(2),'FaceColor',[0 1 0],'LineStyle','none')
    set(bar2(3),'FaceColor',[0 0 1],'LineStyle','none')
    set(gca,'XTick',1:10,'Xdir','reverse')
    ylim([0 120])
    xlim([0.4 10.6]);
    title('Contribution des mécanismes de dépopulation des 2p','fontweight','bold')
    xlabel('Niveau 2p_x','FontSize',11,'fontweight','bold');
    lgnd=legend('Trans rad','Coll 2p_x->2p_y','Coll 2p_x->1s','Orientation','horizontal','Location','north');
    set(lgnd,'FontSize',8)    
    clear depop
    
%% Graphe de la densité des niveaux 1s
      subplot(2,3,1)
      subplot('position',[0.04 0.55 0.29 0.4])
      scatter(energie1s(gaz,:)',densite1s(gaz,1,1,:),'filled','d')
      %errorbar(energie1s(gaz,:)',densite1s(gaz,:)/ng(gaz),sig_densite1s(gaz,:)/ng(gaz),'o','CapSize',18)
      set(gca, 'yscale', 'log')
      xlabel('Énergies Niveaux 1s_5 à 1s_2','FontSize',11,'fontweight','bold');
      title('Densite des niveaux 1s_x m^-3','fontweight','bold')
      if gaz==2
        xlim([16.5 17]);
      elseif gaz==3
        xlim([11.45 11.9]);
      elseif gaz==4
        xlim([9.8 10.7]);
      elseif gaz==5
        xlim([8.2 9.7]);
      end
      
      
          
      %% Graphe des mécanismes de population 1s
    pop1(:,:) = Population1s(gaz,:,:);
    
    subplot(2,3,2)
    subplot('position',[0.38 0.55 0.29 0.36])  
    bar1=bar(1:5,pop1,0.75,'stacked');
    set(bar1(2),'FaceColor',[0 1 0],'LineStyle','none')
    set(bar1(3),'FaceColor',[0 0 1],'LineStyle','none')
    set(bar1(1),'FaceColor',[1 0 1],'LineStyle','none')
    set(gca,'XTick',2:5,'Xdir','reverse')
    ylim([0 120]);
    xlim([1 6]);
    title('Contribution des mécanismes de population des 1s','fontweight','bold')
    xlabel('Niveau 1s_x','FontSize',11,'fontweight','bold');
    lgnd=legend('Coll e^- - fond','Auto-abs','Mixing','Transfert Ar(1s_5)','Orientation','horizontal','Location','north');
    set(lgnd,'FontSize',8)

    clear pop1 

          %% Graphe des mécanismes de depopulation des 1s
    depop1(:,:) = Depopulation1s(gaz,:,:);
    subplot(2,3,3)
    subplot('position',[0.70 0.55 0.29 0.43])
    bar1=bar(1:5,depop1,0.75,'stacked');
    set(bar1(2),'FaceColor',[0 1 0],'LineStyle','none')
    set(bar1(3),'FaceColor',[0 0 1],'LineStyle','none')
    set(bar1(1),'FaceColor',[1 0 1],'LineStyle','none')
    set(bar1(4),'FaceColor',[1 0 0],'LineStyle','none')
    set(bar1(5),'FaceColor',[0 1 1],'LineStyle','none')
    set(bar1(6),'FaceColor',[1 1 0],'LineStyle','none')
    set(gca,'XTick',2:5,'Xdir','reverse')
    ylim([0 120]);
    xlim([1 6]);
    title('Contribution des mécanismes de depopulation des 1s','fontweight','bold')
    xlabel('Niveau 1s_x','FontSize',11,'fontweight','bold');
    lgnd=legend('Radiative','vers 2P','Mixing','Superelastique','Ionization','Neutre','Orientation','horizontal','Location','north');
    set(lgnd,'FontSize',8,'Location','NorthOutside')

    clear depop1 

    
end




end

