%fonction copié collé de PumpSub de Donnelly Teoes 2004.
%calcul la variable Pump qui donne le flow de gaz effectif pour un mélange.
%
function [Pump,Pression_calc]= TePumpingSpeed_TRG(pmtorr,flow,name_Dimension)
%% récupération des dimension du réacteur dans Dimension.txt
    tempfile=fopen(name_Dimension); %dimension du réacteur + pompe
    data = textscan(tempfile,'%s','delimiter',' ','headerlines',9);
    Dimension=str2double(data{1,1});  % Variable avec les dimension du réacteur.
    fclose(tempfile);
    clear data

%% adaptation des paramètres pour avoir les mêmes que Donnelly
igas = find(flow~=0); %donc les numéro de gaz utilisés 13 pour O2, 2 Ne, 3 Ar, etc...
ngases=length(igas); % nombre de gaz
frot = Dimension(1);
rdisk = Dimension(2);
rbase = Dimension(3);
ltube = Dimension(4);
dtube = 2*Dimension(5);
dvalve = Dimension(6);

%% initialisation de tableaux

ktrap      =ones(1,100);
mfpav      =zeros(1,25);
mfpc       =zeros(1,25);
knav       =zeros(1,25);
ndav       =zeros(1,25);
Pump       =zeros(1,25);

        
%% constant
k= 1.38e-16;
na=6.023e23;
r=1.07e-16;
% !	molecular diameters (cm) (from "Vac. Sci. and Tech",p.34., for He-Xe,
% !	N2 (air) and O2. Cl2 is an estimate from the sigma L-J parameter (4.115)
% !	increased some to match the slightly larger mol-d of O2, relative to
% !	its sigma L-J, as well as for CO2.  Likewise, md for Cl is assumed to
% !	be L-J sigma = 3.5.  C2F6 and c-C4F8 are from Sematech report #4. CF4 is
% !	a guess (=CHF3, Sematech #4).  F, BCl3, HBr, CF, CF2, C2F4 and C4F6 are guesses.
% !
	 md=[2.18e-8, 2.60e-8, 3.67e-8, 4.15e-8, 4.91e-8, 5.44e-8, 3.5e-8, 2.6e-8,0,0,...
          5.0e-8, 4.0e-8, 3.64e-8, 3.74e-8,0,0,0,3.7e-8, 3.8e-8, 4.0e-8, ...
          4.1e-8, 4.3e-8, 4.65e-8, 4.8e-8, 4.90e-8];


	%sigg=1.0e-8*sigg	!convert the L-J sig values to cm.

	 sigg=[2.576e-8, 2.789e-8, 3.418e-8, 3.61e-8, 4.055e-8, 4.115e-8,...
         3.548e-8, 2.7e-8, 0, 0, 5.0e-8, 3.7e-8, 3.43e-8, 3.681e-8, 0, 0, 0,...
         3.47e-8, 3.47e-8,4.0e-8, 4.7e-8, 4.67e-8, 4.65e-8,4.8e-8,4.90e-8];
     	
     %	molecular weights

	m=[4.00, 20.18, 39.95, 83.8, 131.8, 70.8, 35.4, 19.0, 0, 0, 117.16, ...
        80.91, 32.0, 28.0, 0, 0, 0, 46.2, 31.01, 50.01, 88.00, 100.00, ...
        138.01, 162.03, 200.02];
	
   	%gamma (Cp/Cv).  C4F6 (igas=24) and C4F8 (25) are guesses

	 gamma=[1.66667,1.66667,1.66667,1.66667,1.66667,1.325,1.66667,1.66667,...
         0, 0, 1.135,1.287,1.396,1.405, 0, 0, 0, 1.266,1.279,1.216 ,...
         1.14,1.104,1.078,1.05,1.05];
     
%     	trapping probabilities (from JVST A 1, 136 (1983).)
% 	First set all trapping probabilities to 1.0, then set those known
% 	to be less than 1.
    ktrap(1:5)=[0.689,0.918,1.0,1.0,1.0];

    
    
%%

	vblade=(4.0/3.0)*pi*frot*(rdisk^3-rbase^3)/(rdisk^2-rbase^2);
	fblade=pi*(rdisk^2-rbase^2);
	tk=300.0;
    
   
	reduce=1.01;
	lod=ltube/dtube ;
    aprimelt=(4.0*dtube)/(3.0*ltube);
    
    
    
 	for i=1:ngases
%	   Inside this loop are things that do not depend on pressure
	   v(i)=sqrt(8.0*k*tk*na/(pi*m(igas(i))));
	   s(i)=0.8*vblade*fblade/(4.0*(1.0/ktrap(igas(i))+vblade/v(i)));
	   k2(i)=(0.5/sqrt(2.0*pi))*((gamma(igas(i))+1.0)/2.0)^(1/(gamma(igas(i))-1))*sqrt((gamma(igas(i))+1.0)/(2.0*gamma(igas(i))));
	   eta(i)=0.499*m(igas(i))*v(i)/(sqrt(2.0)*pi*na*md(igas(i))^2);
	   clongtubem(i)=aprimelt*v(i)*pi*dtube^2/16.0;
	   cvalvem(i)=pi*v(i)*dvalve^2/16.0;
	   capertubem(i)=pi*v(i)*dtube^2/16.0;
    end
    
 for i=1:ngases
	   for j=1:ngases		
		   phi(i,j)=1.0/sqrt(8.0+8.0*m(igas(i))/m(igas(j)))*(1.0+sqrt(eta(i)/eta(j))*(m(igas(j))/m(igas(i)))^0.25)^2;		
       end
 end  
 totalflow=sum(flow);
% !	Make a crude estimation of the pressure of each species in the chamber 
% !	and above the valve to begin iterations below.
% !
	for i=1:ngases
	   pc(i)=(flow(igas(i))/totalflow).*pmtorr;
	   pv(i)=pc(i)/2.0;
	   ndc(i)=3.106e13*(298.16/tk)*pc(i);
	   ndv(i)=ndc(i)/2.0;
    end

    %%
    %	ITERATE FROM HERE TO GET REDUCE
    iter=0; %compte le nombre d'itération
	ptot=0.0;
 while (ptot/pmtorr < 0.99 && reduce > 1.0)
%	do ireduce=1,10

%	 Iterate 5 times to compute pressure above the valve (pv)
%	 Mean-free-paths for all gases in system

	 for j=1:5
	   ndvtot=0.0;
	   for i=1:ngases
		ndv(i)=3.106e13*(298.16/tk)*pv(i);
		ndvtot=ndvtot+ndv(i);
        end

	   sigmix3=0.0;
	   avermass=zeros(1,25); %nombre de gaz possibles
	   for i=1:ngases
		for i2=1:ngases
		   if i2~=i
           avermass(igas(i)) = avermass(igas(i)) + m(igas(i2))*ndv(i2)/(ndvtot-ndv(i));
           sigpair3=(0.5*(sigg(igas(i))+sigg(igas(i2))))^3;
		   fracts=(ndv(i)/ndvtot)*(ndv(i2)/ndvtot);
		   sigmix3=sigmix3+fracts*sigpair3;
           end
         end
       end
	   sigmix=(sigmix3)^0.333333;

	   for i=1:ngases
		mfpv(i)=1.0/(sqrt(2.0)*pi*ndv(i)*sigg(igas(i))^2 + pi*(ndvtot-ndv(i))*((sigg(igas(i))+sigmix)/2.0)^2*(1.0+m(igas(i))/avermass(igas(i))));
	 	knv(i)=mfpv(i)/dvalve;
	 	cvalvet(i)=(cvalvem(i)/reduce)*(10.0+0.5*(dvalve/(reduce*mfpv(i)))^1.5)/(10.0+k2(i)*(dvalve/(reduce*mfpv(i)))^1.5);
	 	pv(i)=r*tk*4.482e17*flow(igas(i))*(1.0/s(i)+1.0/cvalvet(i));
        end
     end
	   pvtot=sum(pv);


%	   Iterate 5 times to compute pressure in chamber (pc)

	 for j=1:5
	   pctot=0.0;
	   ndctot=0.0;
	   for i=1:ngases
		ndc(i)=3.106e13*(298.16/tk)*pc(i);
	 	pctot=pctot+pc(i);
		ndctot=ndctot+ndc(i);
       end

	   sigmix3=0.0;
	   avermass=zeros(1,25); %nombre de gaz possibles
	   for i=1:ngases
		for i2=1:ngases
		   if(i2 ~= i) 
           avermass(igas(i)) = avermass(igas(i)) + m(igas(i2))*ndc(i2)/(ndctot-ndc(i));
		   sigpair3=(0.5*(sigg(igas(i))+sigg(igas(i2))))^3;
		   fracts=(ndc(i)/ndctot)*(ndc(i2)/ndctot);
		   sigmix3=sigmix3+fracts*sigpair3;
           end
         end
        end
	   sigmix=(sigmix3)^0.333333;
	   ndavtot=0.0;
       
       
       
	   for i=1:ngases

		mfpc(i)=1.0/(sqrt(2.0)*pi*ndc(i)*sigg(igas(i))^2 + pi*(ndctot-ndc(i))*((sigg(igas(i))+sigmix)/2.0)^2*(1.0+m(igas(i))/avermass(igas(i))));
		
	   	knc(i)=mfpc(i)/dtube;
%		Use mfp and Kn for average of pressures in chamber and
%		above valve:
		mfpav(i)=(mfpc(i)+mfpv(i))/2.0;
		knav(i)=(knc(i)+knv(i))/2.0;
		ndav(i)=(ndc(i)+ndv(i))/2.0;
		ndavtot=ndavtot+ndav(i);
       end
%	   Compute viscosity for the mixture

	   sumdenom=0.0;
	   etamix=0.0;
	   for ieta=1:ngases
		for jeta=1:ngases
		   sumdenom=sumdenom+(ndav(jeta)/ndavtot)*phi(ieta,jeta);
        end
		etamix=etamix+(ndav(ieta)/ndavtot)*eta(ieta);
       end
	   for i=1:ngases
		capertubet=capertubem(i)*(10.0+0.5*(dtube/mfpav(i))^1.5)/(10.0+k2(i)*(dtube/mfpav(i))^1.5);	%???should it be k2(i) or k2(6)???
%		1.33 converts pressure to dynes/cm2
		clongtubev=1.33*(pctot+pvtot)*pi*dtube^4/(256.0*etamix*ltube);
		clongtubet=clongtubev+clongtubem(i)*(knav(i)+sqrt(pi/2.0))/(knav(i)+1.235*sqrt(pi/2.0));
		ctubet=1.0/(1.0/clongtubet+1.0/capertubet);

%		BOTTOM LINE:

		pc(i)=r*tk*4.482e17*flow(igas(i))*(1.0/s(i)+1.0/cvalvet(i)+1.0/ctubet);        
        end
     end
	   ptot=sum(pc);
	 %pump(99)=reduce	%Sneak reduce back to main program as pump(99)
	 reduce=reduce*pmtorr/ptot;
     iter=iter+1; %compte le nb d'itérations
 end
 %iter
 Pression_calc=ptot; % pression effective calculer par le programme
% 
% pump(igas(i)) is the inverse of the efficiency of pumping of species i.
% !	Number densities should be multiplied by pump(igas(i))
% !	to correct for pumping speed effects.
% !	First compute relative to Xe
% !
	for i=1:ngases
        if ngases>4
	   pumprel(igas(i))=(pc(i)/pc(5))*(flow(5)/flow(igas(i)));
        else
       pumprel(igas(i))=(pc(i)/pc(4))*(flow(5)/flow(igas(i))) ;   
        end
    end
    
% Now compute normalization factor:

	sumfs=0.0;
	for i=1:ngases
	   sumfs=sumfs+flow(igas(i))*pumprel(igas(i));
    end
	an=totalflow/sumfs;

	for i=1:ngases
	   Pump(igas(i))=an*pumprel(igas(i));
    end   
    
    
end