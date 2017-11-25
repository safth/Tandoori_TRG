function[Int_time_correction] = IntTime_trtx(fileNames)         

%Cette fonction va chercher les temps d'integrations des deux spectromètres
%Avantes 500-700 et 700-900 synchronisé. elle fait ensuite un facteur de
%correction pour les raies du 500-700. pour les fichiers .trtx
  
        tempfile=fopen(fileNames);
           tempdata = textscan(tempfile,'%s %s','delimiter',';');
           int_time1 = tempdata{1,1};
           int_time1 = int_time1{6,1}; %peut etre 6 et 7 ou [5 6]
           int_time_selector = str2double(int_time1(1)); %regarde si c'est le 700-99 ou le 500-700
           int_time1 = int_time1(54:strfind(int_time1,'.00')+2); %prend le temps d'integration
           time(1) = str2double(int_time1);
           int_time2 = tempdata{1,1};
           int_time2 = int_time2{7,1}; %peut etre 6 et 7 ou  5 6 
           int_time2 = int_time2(54:strfind(int_time2,'.00')+2); %prend le temps d'integration
           time(2) = str2double(int_time2);
                %si ce n'est pas les bonnes lignes (Arrive parfois)
           if isnan(time(1)) | isnan(time(2))
                clear int_time1 int_time2 time
                int_time1 = tempdata{1,1};
                int_time1 = int_time1{5,1}; %peut etre 6 et 7 ou [5 6]
                int_time_selector = str2double(int_time1(1)); %regarde si c'est le 700-99 ou le 500-700
                int_time1 = int_time1(54:strfind(int_time1,'.00')+2); %prend le temps d'integration
                time(1) = str2double(int_time1);
                int_time2 = tempdata{1,1};
                int_time2 = int_time2{6,1}; %peut etre 6 et 7 ou  5 6 
                int_time2 = int_time2(54:strfind(int_time2,'.00')+2); %prend le temps d'integration
                time(2) = str2double(int_time2);
           end
           
          %====faceur de correction pour les temps d'intégrations diff des 2 spectros
          if int_time_selector == 7 %si le premier est le 700=900
             Int_time_correction = time(1)/time(2); %I_exp = I_exp*correction
          else
             Int_time_correction = time(2)/time(1); %I_exp = I_exp*correction
          end
          fclose(tempfile);
          clear tempdata int_time_selector int_time1 int_time2 time