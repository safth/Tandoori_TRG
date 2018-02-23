%% Cette fonction n'est pas utilisé

function TeGraphes1(ne,Te,I_theo,I_exp,FitCorrige,moyenne,Neopt)

I_exp=I_exp(I_exp~=0);
tempI_exp(:)=I_exp(1,:);
% 
% 
% I_simul=zeros(1,41);
% for m=1:41
%     I_simul(m)=I_theo(j,i,m)*FitCorrige(m);
% end
% I_simul=I_simul(I_simul~=0);
%% ==============================================================================
%% ====== Cette fonction trace les graphiques des rapports Iobs/Iexp ======
%% ======================= avec Ne ou Te de fixe  =====================
%% ==============================================================================
 %I_theo(Neopt,40,:)
 for i=1:length(Te)
     
     tempI_theo(:)=I_theo(Neopt,i,:);
     tempI_theo(:)= tempI_theo(:)*moyenne;
     tempI_theo(:)= tempI_theo.*FitCorrige;
     tempI_theo=tempI_theo(tempI_theo~=0);

     Ratio(i,:)= I_exp./tempI_theo;
     clear tempI_theo
     
     for j=1:length(I_exp) % on veut un minimum en Ratio=1
        if Ratio(i,j) < 1 
         Ratio(i,j)=1/Ratio(i,j);
        end
     end
     
 end
 
  figure
    plot(Te,Ratio)
    set(gca, 'yscale', 'log')
    set(gcf,'color','w');
    xlabel('Te','FontSize',11,'fontweight','bold'); ylabel('I_{exp}/I_{theo}','FontSize',11,'fontweight','bold');
    title('I_{exp}/I_{theo} pour le premier Ne donnée','FontSize',12,'fontweight','bold')
    ylim([-1 100]);
     legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'...
        ,'21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41')
%   for i=1:length(ne)
%      
%      tempI_theo(:)=I_theo(i,1,:);
%      tempI_theo(:)= tempI_theo(:)*moyenne;
%      tempI_theo(:)= tempI_theo.*FitCorrige;
%      tempI_theo=tempI_theo(tempI_theo~=0);
% 
%      Ratio(i,:)= I_exp./tempI_theo;
%      clear tempI_theo
%      
%      for j=1:length(I_exp) % on veut un minimum en Ratio=1
%         if Ratio(i,j) < 1 
%          Ratio(i,j)=1/Ratio(i,j);
%         end
%      end
%      
%  end
%     figure
%     plot(ne,Ratio)
%     set(gca, 'yscale', 'log')
%     set(gca, 'xscale', 'log')
%     set(gcf,'color','w');
%     xlabel('Ne','FontSize',11,'fontweight','bold'); ylabel('I_{exp}/I_{theo}','FontSize',11,'fontweight','bold');
%     title('I_{exp}/I_{theo} pour le premier Ne donnée','FontSize',12,'fontweight','bold')

end
