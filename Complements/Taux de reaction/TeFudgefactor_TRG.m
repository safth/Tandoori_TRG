function [FUDGEground_2P,FUDGE1s_2P] = TeFudgefactor_TRG(gaz)

%Matrice avec les fudgeFactor de Donnelly dans Teoes 2004

%Tous les fudge sont 1 pour les excitation vers les 1s et le quenching
%Donc seulement pour les excitations vers les 2p (Gs -> 2p et 1s -> 2p)

    if gaz ==3          %% Argon                       
        FUDGEground_2P=   [1.1  1  1  1.3  1.25  1  1  1  1  1.1]'; % 2p1 ... 2p10
        

        FUDGE1s_2P(3,:) = [1 1 1 1.2 1 1 1 1 1 1];      %1s3
        FUDGE1s_2P(5,:) = [1 1 1 1.2 1 1 1 0.9 1 1];    %1s5
        FUDGE1s_2P(1,:) = [1 1 1 1 1 1 1 1 1 1];        %RIEN
        FUDGE1s_2P(2,:) = [1 1 1 1 1 1 1 1 1 1];        %1s2
        FUDGE1s_2P(4,:) = [1 1 1 1 1 1 1 1 1 1];        %1s4

    elseif gaz ==2      %% Neon
        FUDGEground_2P =  [1 1 1 1 1 1 1 1 1 1]';                % 2p1 ... 2p10

        FUDGE1s_2P(3,:) = [1 1 1 1 1 1 1 1 1 1]; %1s3
        FUDGE1s_2P(5,:) = [1 1 1 1 1 1 1 1 1 1]; %1s3
        FUDGE1s_2P(1,:) = [1 1 1 1 1 1 1 1 1 1]; %RIEN
        FUDGE1s_2P(2,:) = [1 1 1 1 1 1 1 1 1 1]; %1s2
        FUDGE1s_2P(4,:) = [1 1 1 1 1 1 1 1 1 1]; %1s4

    elseif gaz ==4      %% Krypton
        FUDGEground_2P =  [0.92 0.7 0.5 0.8 0.94 0.7 0.6 0.8 0.7 0.7]';  % 2p1 ... 2p10

        FUDGE1s_2P(3,:) = [1.0 1.0 0.3 1.0 1.0 1.0 0.6 1.0 1.0 1.0];      %1s3
        FUDGE1s_2P(5,:) = [1.0 1.0 1.0 1.0 1.0 1.5 0.6 1.8 1.8 1.0];%1s5
        FUDGE1s_2P(1,:) = [1 1 1 1 1 1 1 1 1 1];        %RIEN
        FUDGE1s_2P(2,:) = [1 1 1 1 1 1 1 1 1 1];        %1s2
        FUDGE1s_2P(4,:) = [1 1 1 1 1 1 1 1 1 1];        %1s4

    elseif gaz ==5      %% Xenon
        FUDGEground_2P =  [1.03 1.1 0.90 1.1 1.26 1 1.1 1.3 1 1.1]';           % 2p1 ... 2p10
                     
        FUDGE1s_2P(3,:) = [0.05 6.5 1 10 0.3 1 1 1 1 1];         %1s3
        FUDGE1s_2P(5,:) = [0.05 1 1 1  0.3 1 2 0.8 1.1 1];%1s5
        FUDGE1s_2P(1,:) = [1 1 1 1 1 1 1 1 1 1];         %RIEN
        FUDGE1s_2P(2,:) = [1 1 1 1 1 1 1 1 1 1];         %1s2
        FUDGE1s_2P(4,:) = [1 1 1 1 1 1 1 1 1 1];         %1s4
    end
end

