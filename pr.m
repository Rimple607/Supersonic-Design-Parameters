p_r = [50,100,150,200,250,300,350,400,450,500];
disp(['p_r : ',num2str(p_r)]);
Pxo = 0.875;        ga = 1.4;      en = 0.025;     Tr = 1;     px1 = 4;

Pix = p_r.*(1/Pxo);
%disp(['Pix : ',num2str(Pix)]);

Mx1 = (((Pix.^((ga-1)/ga)).-1).*(2/(ga-1))).^(1/2);
Mx2 = (((1/Pxo)^((ga-1)/ga) - 1)*(2/(ga-1)))^(1/2);
%disp(['Mx1 : ',num2str(Mx1)]);
%disp(['Mx2 : ',num2str(Mx2)]);

Mx1_c = ((((ga+1)/2).*(Mx1.^2))./(1.+(((ga-1)/2).*(Mx1.^2)))).^(1/2);
Mx2_c = ((((ga+1)/2)*(Mx2^2))/(1 + (ga-1)*(Mx2^2)*(1/2)))^(1/2);
%disp(['Mx1_c : ',num2str(Mx1_c)]);
%disp(['Mx2_c : ',num2str(Mx2_c)]);

M1_c = (Mx1_c.+(en*Mx2_c*(Tr^(1/2))))./((1 + en)*(1 + en*Tr))^(1/2);
%disp(['M1_c : ',num2str(M1_c)]);

M1 =  ((2.*(M1_c.^2))./((ga+1).-((ga-1).*(M1_c.^2)))).^(1/2);
%disp(['M1 : ',num2str(M1)]);

pb = px1.*((((2*ga)/(ga+1)).*(M1.^2)).-((ga-1)/(ga+1)));
%disp(['pb : ',num2str(pb)]);

R = pb./px1;
disp(['R : ',num2str(R)]);
plot (p_r, R);
xlabel (' PR ');
ylabel (' RR ');








fprintf('\n');

