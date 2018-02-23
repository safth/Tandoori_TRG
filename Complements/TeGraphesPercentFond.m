%% fonction pas utilisé, permet d'avoir le pourcentage du fondamental pour tous les nievaux de tous les gaz!!

function  TeGraphesTEMP(gaz,Ne,Te,ContributionFond)
% if gaz==2
%     figure('name','Neon')
% elseif gaz==3
%     figure('name','Argon')
% elseif gaz==4
%     figure('name','Krypton') 
% elseif gaz==5
%     figure('name','Xénon')    
% end
%% graphes 3D
% subplot(2,5,1)
%      surf(Ne,Te,ContributionFond(:,:,1)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P1','fontweight','bold')
%      zlim([20 100])
% subplot(2,5,2)
%      surf(Ne,Te,ContributionFond(:,:,2)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P2','fontweight','bold')
%      zlim([20 100])
% subplot(2,5,3)
%      surf(Ne,Te,ContributionFond(:,:,3)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P3','fontweight','bold')
%      zlim([20 100])
% subplot(2,5,4)
%      surf(Ne,Te,ContributionFond(:,:,4)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P4','fontweight','bold')   
%      zlim([20 100])
% subplot(2,5,5)
%      surf(Ne,Te,ContributionFond(:,:,5)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P5','fontweight','bold')   
%      zlim([20 100])
% subplot(2,5,6)
%      surf(Ne,Te,ContributionFond(:,:,6)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P6','fontweight','bold')  
%      zlim([20 100])
% subplot(2,5,7)
%      surf(Ne,Te,ContributionFond(:,:,7)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P7','fontweight','bold')   
%      zlim([20 100])
% subplot(2,5,8)
%      surf(Ne,Te,ContributionFond(:,:,8)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P8','fontweight','bold')  
%      zlim([20 100])
% subplot(2,5,9)
%      surf(Ne,Te,ContributionFond(:,:,9)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P9','fontweight','bold')
%      zlim([20 100])
% subplot(2,5,10)
%      surf(Ne,Te,ContributionFond(:,:,10)','edgecolor','none') 
%      colormap(hsv)
%      set(gca, 'xscale', 'log')
%      xlabel('Densité électronique(m-3)','FontSize',11,'fontweight','bold'); ylabel('Te(eV)','FontSize',11,'fontweight','bold');
%      zlabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%      title('2P10','fontweight','bold')   
%      zlim([20 100])
     

        %% plot 2D en fonction de Te
% 
%       plot(Te,ContributionFond(1,:,1),'LineWidth',2)
%       xlabel('Te(eV)','FontSize',11,'fontweight','bold');
%       ylabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
%       hold on
%       plot(Te,ContributionFond(1,:,2),'LineWidth',2)
%       plot(Te,ContributionFond(1,:,3),'LineWidth',2)
%       plot(Te,ContributionFond(1,:,4),'LineWidth',2)
%       plot(Te,ContributionFond(1,:,5),'LineWidth',2)
%       plot(Te,ContributionFond(1,:,6),'LineWidth',2,'Linestyle','--')
%       plot(Te,ContributionFond(1,:,7),'LineWidth',2,'Linestyle','--')
%       plot(Te,ContributionFond(1,:,8),'LineWidth',2,'Linestyle','--')
%       plot(Te,ContributionFond(1,:,9),'LineWidth',2,'Linestyle','--')
%       plot(Te,ContributionFond(1,:,10),'LineWidth',2,'Linestyle','--')
%       lgnd=legend('2P1','2P2','2P3','2P4','2P5','2P6','2P7','2P8','2P9','2P10','Orientation','horizontal','Location','north');
%     set(lgnd,'FontSize',8)
%       hold off


        %% plot 2D en fonction de Ne
        %AVEC LE PREMIER TE
        name = {'Néon','Argon','Krypton','Xénon'}
fig_temp = figure;
for gaz=2:5
      subplot(2,2,gaz-1)
      plot(Ne,ContributionFond(gaz,:,1,1),'LineWidth',2)
      xlabel('n_e(m^{-3})','FontSize',11,'fontweight','bold');
      ylabel({'% du fondamental'},'FontSize',10,'fontweight','bold')
      title(name{gaz-1})
      hold on
      plot(Ne,ContributionFond(gaz,:,1,2),'LineWidth',2)
      plot(Ne,ContributionFond(gaz,:,1,3),'LineWidth',2)
      plot(Ne,ContributionFond(gaz,:,1,4),'LineWidth',2)
      plot(Ne,ContributionFond(gaz,:,1,5),'LineWidth',2)
      plot(Ne,ContributionFond(gaz,:,1,6),'LineWidth',2,'Linestyle','--')
      plot(Ne,ContributionFond(gaz,:,1,7),'LineWidth',2,'Linestyle','--')
      plot(Ne,ContributionFond(gaz,:,1,8),'LineWidth',2,'Linestyle','--')
      plot(Ne,ContributionFond(gaz,:,1,9),'LineWidth',2,'Linestyle','--')
      plot(Ne,ContributionFond(gaz,:,1,10),'LineWidth',2,'Linestyle','--')
      lgnd=legend('2P1','2P2','2P3','2P4','2P5','2P6','2P7','2P8','2P9','2P10','Orientation','horizontal','Location','north');
    set(lgnd,'FontSize',8)
    set(gca, 'xscale', 'log') 
      hold off
end    
      
      %%
% if gaz==2
%     a(:,:)=ContributionFond(gaz,:,1,:);
%     save('Neon.mat','a')
%     b(:) = Ne;
%     save('Ne.mat','b')
% elseif gaz==3
%     a(:,:)=ContributionFond(gaz,:,1,:);
%     save('Argon.mat','a')
% elseif gaz==4
%     a(:,:)=ContributionFond(gaz,:,1,:);
%     save('Krypton.mat','a')
% elseif gaz==5
%     a(:,:)=ContributionFond(gaz,:,1,:);
%     save('Xenon.mat','a')   
% end


      