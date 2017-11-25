%% ============================================================================
%% ======= Cette fonction permet le calcul initial des valeurs de beta ========
%% ============== Pour l'integrale utilisé dans le calcul de la ===============
%% ================ convolution du profil Voight d'élargissement ==============
%% ===================== dans TeEscapeFactor.m ================================
%% ============================================================================

Beta = logspace(-5,3,2000);

for j=1:length(Beta)
    
        fun=@(y,w) exp(-y.^2)./((Beta(j).^2)+(w-y).^2);
        integrale(j)=integral2(fun,-inf,inf,-inf,inf);
    
    
    
end

AllIntegral = [Beta',integrale'];
save('AllIntegral_new.mat','AllIntegral')