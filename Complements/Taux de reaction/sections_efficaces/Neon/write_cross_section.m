%Write down cross section fichiers

%%1si - 2pj
for i=3:2:5 %écrit 1s3 et 1s5
    for j=1:10 %écrit les 2P1 à 2P10
        fileID = fopen( sprintf('1s%01d-2p%01d.txt',i,j),'at'); %ouvre ou crée un fichier en mode txt avec comme permission décrire un a la suite de lautre
        formatSpec = '%s \r\n'; %string (/r/n=saut de page sur windows)
        fclose(fileID);
    end
end
 %% Ground - 2pj
for j=1:10 %écrit les 2P1 à 2P10
    fileID = fopen( sprintf('gs-2p%01d.txt',j),'at'); %ouvre ou crée un fichier en mode txt avec comme permission décrire un a la suite de lautre
    formatSpec = '%s \r\n'; %string (/r/n=saut de page sur windows)
    fclose(fileID);
end

%% Quenching 1si=1sj
for i=3:2:5 %écrit 1s3 et 1s5
    for j=2:5 %écrit les 2P1 à 2P10
        if j ~=i
            fileID = fopen( sprintf('1s%01d-1s%01d.txt',i,j),'at'); %ouvre ou crée un fichier en mode txt avec comme permission décrire un a la suite de lautre
            formatSpec = '%s \r\n'; %string (/r/n=saut de page sur windows)
            fclose(fileID);
        elseif j==i
            %1si -> ground state
            fileID = fopen( sprintf('1s%01d-gs.txt',i),'at'); %ouvre ou crée un fichier en mode txt avec comme permission décrire un a la suite de lautre
            formatSpec = '%s \r\n'; %string (/r/n=saut de page sur windows)
            fclose(fileID);
            
            %ground->1si
            fileID = fopen( sprintf('gs-1s%01d.txt',i),'at'); %ouvre ou crée un fichier en mode txt avec comme permission décrire un a la suite de lautre
            formatSpec = '%s \r\n'; %string (/r/n=saut de page sur windows)
            fclose(fileID);
            
            %ionisation
            fileID = fopen( sprintf('1s%01d-Ionization.txt',i),'at'); %ouvre ou crée un fichier en mode txt avec comme permission décrire un a la suite de lautre
            formatSpec = '%s \r\n'; %string (/r/n=saut de page sur windows)
            fclose(fileID);
        end
    end
    
end