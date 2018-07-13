clear all
siMIN=0;
siMAX=1;
dsi=0.01;
srMIN=0.0;
srMAX=1.0;
dsr=0.05;
SI=siMIN:dsi:siMAX;
nsi=round(1+((siMAX-siMIN)/dsi));
SR=srMIN:dsr:srMAX;
nsr=round(1+((srMAX-srMIN)/dsr));
    for ISI=1:nsi;
        for ISR=1:nsr;
            
            S=SR(ISR)+1i*SI(ISI);
           
            G1=pi^(S/2)/(gammai(S/2));
            G2=pi^((conj(S))/2)/(gammai((conj(S))/2));
            G=G1-G2;
            F(ISI,ISR)=abs(G);
        end;
    end;
    contour(SR, SI, F, 100)