clear all 

clc

        file1=fopen('NEW0001.TXT');
        tempdata=textscan(file1, '%s\t%s\t%s\t%s','delimiter',';', 'headerlines',8);  % Extraction des datas 
        lambdaTe=tempdata{1}  ;                 % Longueurs d'ondes des mesures
        lambdaTe=str2double(strrep(lambdaTe,',','.'));
        intTe=tempdata{2};                      % Intensité des mesures
        intTe=str2double(strrep(intTe,',','.'));
        clear tempdata                          % On a plus besoin de tempdata 