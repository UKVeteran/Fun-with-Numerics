
clear all
siMIN=0;
siMAX=10.0;
dsi=0.1;
srMIN=0;
srMAX=10.0;
dsr=0.1;
SI=[siMIN:dsi:siMAX];
nsi=round(1+((siMAX-siMIN)/dsi));
SR=[srMIN:dsr:srMAX];
nsr=round(1+((srMAX-srMIN)/dsr));
    for ISI=1:nsi;
        for ISR=1:nsr;
            
            S=SR(ISR)+1i*SI(ISI);
           
            G(ISR,ISI)=((pi.^(S/2))/(gammai((S/2))))+((pi.^((1-S)/2))/(gammai((1-S)/2)));
        end;
    end;
    contour(SI, SR, G, 50);
   