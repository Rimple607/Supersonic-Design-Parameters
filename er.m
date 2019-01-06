w = [0.2,0.15,0.1,0.08,0.07,0.06,0.05,0.04,0.03,0.025];
disp(w);

Mx1 = 2.22;     Mx2 = 0.52;     Tr = 1;    ga = 1.4;      px = 4;

m1c = (Mx1_c.+((Mx2_c*(Tr^(1/2))).*w))./(((1.+w).*(1.+(Tr.*w))).^(1/2));
%disp(['m1c : ',num2str(m1c)]);

m1 =  ((2.*(m1c.^2))./((ga+1).-((ga-1).*(m1c.^2)))).^(1/2);
%disp(['m1 : ',num2str(m1)]);

m1b = (((2/(ga-1)).+(m1.^2))./((((2*ga)/(ga-1)).*(m1.^2)).-1)).^(1/2);
%disp(['m1b : ',num2str(m1b)]);

p1b = px.*((((2*ga)/(ga+1)).*(m1.^2)).-((ga-1)/(ga+1)));
%disp(['p1b : ',num2str(p1b)]);

R = p1b./px;
disp(['R : ',num2str(R)]);

plot (w, R)
xlabel (' ER ');
ylabel (' RR ');

fprintf('\n');





