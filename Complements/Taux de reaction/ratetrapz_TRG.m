function [result]= ratetrapz(Te, file, exposant,HighPressure,Niveau)

%% intèegre sur une fonction de distribution généralisé la section efficace choisis.

%Definition des constantes
e=1.602176462*10^(-19);
me=9.10938188*10^(-31);

%Chargement des sections efficaces
Data=load(file);
energie=Data(:,1);           %Energy is 1st column of the file
sectioneff=Data(:,2);        %Cross section is 2nd column

%On regarde si les sections efficaces ont besoin d'être ramenées en m^2
if max(sectioneff)>1
    sectioneff=10^(-23)*sectioneff;
end

%On regarde si on doit corriger pour la haute pression
a=[1 1.2 0.8 1.4 2.92 1.05 1 0.65 1 1.23];
b=[0 0.06 0.04 0.06 0.0785 0.035 0.785 0.06 0 0.08];
if HighPressure==1
    for i=1:length(energie)
        if energie(i)>14
            frel=1+a(Niveau)*(1-exp(b(Niveau)*(14-energie(i))));
            sectioneff(i)=frel*sectioneff(i);
        end
    end
end

%On définit la EEDF
B1=exposant.*(2/3).^(3/2).*(gamma(5./(2.*exposant))).^(3/2)./(gamma(3./(2.*exposant))).^(5/2);
B2=(2/3).*gamma(5./(2.*exposant))./gamma(3./(2.*exposant));


sigmaE = @(E) (interp1(energie,sectioneff,E).*E.*B1*Te^(-3/2).*exp(-(B2*E./Te).^exposant));

result = (sqrt(2*e/me)).*(integral(sigmaE,min(energie),max(energie)));
end