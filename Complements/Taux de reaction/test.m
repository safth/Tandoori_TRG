
exposant=1;
Te=[0.1:0.1:10];
for i=1:length(Te)
Rate_exposant(i)=ratetrapz_TRG(Te(i),'sections_efficaces/Argon/gs-2p1.txt',   exposant, 0, 0)



fun= @(E) (exp(-E./Te(i))).*E.^(1/2);
norm=integral(fun,0,inf);

Rate_Maxwell(i)=ratetrapz_old(Te(i),'sections_efficaces/Argon/gs-2p1.txt', norm)
%Rate_Maxwell_norm=ratetrapz_test(Te,'sections_efficaces/Argon/gs-2p1.txt', 1)
end

plot(Te,Rate_exposant,'r')
hold
plot(Te,Rate_Maxwell,'b')