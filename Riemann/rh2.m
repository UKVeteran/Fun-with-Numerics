clear all
siMIN=13.000;
siMAX=15.000;
dsi=0.01;
srMIN=0.000;
srMAX=2.000;
dsr=0.05;
SI=siMIN:dsi:siMAX;
nsi=round(1+((siMAX-siMIN)/dsi));
SR=srMIN:dsr:srMAX;
nsr=round(1+((srMAX-srMIN)/dsr));
    for ISI=1:nsi;
        for ISR=1:nsr;
            
            S=SR(ISR)+1i*SI(ISI);
           
            G1=pi^(S/2)/(gammai(S/2));
            S2=SR(ISR)-1i*SI(ISI);
            G2=pi^((1-S)/2)/(gammai(S2/2));
            G=G2+G1;
            F(ISI,ISR)=abs(G);
        end;
    end;
    contour(SR, SI, F, 200)