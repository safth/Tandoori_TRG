% from paper Donnelly 2004 and mullen program paper 
%ve(E)=sqrt(2*e*E/me), and k= integration (cross section*EEDF*ve(E))
% EEDF=2*(pi^-0.5)*(1/k*Te^1.5)*(E^.5)exp(-E/k*Te)

function [result]= ratetrapz(Te, file, norm)
He=load(file);
energie=He(:,1);                    %Energy is 1st column of the file
sectionHe=10^(-23)*He(:,2);                  %Cross section is 2nd column
e=1.602176462*10^(-19);
me=9.10938188*10^(-31);             %Mass of electron       

%fHe=2*sqrt(energie/pi)*(1/Te)^(3/2).*exp(-energie./Te);
fHe=exp(-energie./Te);%sqrt((e)/(pi*me))*(2/Te)^(3/2).*
YHe=energie.*fHe.*sectionHe;

 %result=sqrt(2/me)*trapz(energie,YHe)/trapz(E_norm, ZHe);
 result=(sqrt(2*e/me)).*trapz(energie,YHe)/norm;
 
