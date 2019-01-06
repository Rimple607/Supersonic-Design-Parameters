%% study of pressure recovery
disp('For first stage : ');      fprintf('\n');

%% Writing Cavity conditions and given values
pc = 3;         ga = 1.4;       mc = 1.5;        At = 10;          Pp = 16;       L1 = 6;   
re = 18;        Cd = 0.9;       mf = 0.12;       To = 400;         Ds = 75;       L2 = 6;
en = 0.025;     Pr = 380;       Tr = 1;          Pxo = 0.875;

%% Stage 1
%% Calculating Mach number                             
Ma = fzero(@(M) flowisentropic(M), [1; 4]);         
disp(['mach number is Ma :', num2str(Ma)]);
  
%% Calculating stagnation pressure
Po = pc*(1 + (ga-1)*(mc^2)*(1/2))^(ga/(ga-1));
disp(['Stagnation pressure is : ',num2str(Po)]);

%% Calculating Pseudo throat area
Apt = (At*Pp)/Po;
disp(['Pseudo throat height is : ',num2str(Apt)]);

%% Calculating Mach number downstream of shock
Mb = sqrt((2/(ga-1) + Ma^2)/( ((2*ga)/(ga-1))*(Ma^2) - 1));
disp(['Mach number downstream is : ',num2str(Mb)]);

%% Stagnation pressure downstream of shock
Pob = Po*(1 + ((2*ga)/(ga+1))*(Ma^2 - 1))*((1 + (1/2)*(ga-1)*(Mb^2))/(1 + (1/2)*(ga-1)*(Ma^2)))^(ga/(ga-1));
disp(['Stagnation pressure downstream is : ',num2str(Pob)]);

%% Calculating static pressure downstream
Pb = Pob/((1 + (1/2)*(ga-1)*(Mb^2))^(ga/(ga-1)));
disp(['Static pressure downstream is : ',num2str(Pb)]);     

Px =  (Pb*95)/100;
Pix = Pr/Pxo;
Pxi = 1/Pix;
P3o = re*(Pxo);

%% Calculating Mach number of Suction and Motive gas at x
Mx1 = sqrt(((Pix)^((ga-1)/ga) - 1)*(2/(ga-1))); 
disp(['Motive Gas Mach number at x is : ',num2str(Mx1)]);     
          
Mx2 = sqrt(((1/Pxo)^((ga-1)/ga) - 1)*(2/(ga-1)));
disp(['Suction Gas Mach number at x is : ',num2str(Mx2)]);
    
%% Calculating Characteristic Mach number for both gases
Mx1_c = sqrt((((ga+1)/2)*(Mx1^2))/(1 + (ga-1)*(Mx1^2)*(1/2)));
disp(['Motive Gas characteristic Mach number at x is : ',num2str(Mx1_c)]);     

Mx2_c = sqrt((((ga+1)/2)*(Mx2^2))/(1 + (ga-1)*(Mx2^2)*(1/2)));
disp(['Suction Gas characteristic Mach number at x is : ',num2str(Mx2_c)]);

%% Calculating characteristic Mach number of mixed stream from individual mach numbers
M1_c = (Mx1_c + en*Mx2_c*sqrt(Tr))/sqrt((1 + en)*(1 + en*Tr));
disp(['Mixed Stream characteristic Mach number : ',num2str(M1_c)]);

%% Calculating Mixed Stream Mach number
M1 =  sqrt((2*(M1_c^2))/((ga+1) - (ga-1)*(M1_c^2)));       
disp(['Mixed Stream Mach number : ',num2str(M1)]);

%% Calculating Mach number downstream of Transverse shock 
M1b = sqrt((2/(ga-1) + M1^2)/(((2*ga)/(ga-1))*(M1^2) - 1));
disp(['Mach number downstream of Transverse shock : ',num2str(M1b)]);

%% Calculating Static pressure at the exit
P1b = Px*(((2*ga)/(ga+1))*(M1^2) - (ga-1)/(ga+1));
disp(['Static pressure downstream of Transverse shock : ',num2str(P1b)]);    

P1b_l = (P1b*80)/100;   % After considering Frictional loss
disp(['After considering frictional losses, Static pressure at exit is :', num2str(P1b_l)]);

%% Calculating Stagnation pressure downstream of Transverse shock
Po1b = P1b_l*((1 + (ga-1)*(1/2)*(M1b^2))^(ga/(ga-1)));
disp(['Stagnation pressure downstream of Transverse shock : ',num2str(Po1b)]);

%%Calculating suction to motive gas area ratio
Asm = (Pr)*((((Pxi)^(1/ga))*sqrt(1 - (Pxi)^((ga-1)/ga)))/((Pxo)^(1/ga)*sqrt(1 - (Pxo)^((ga-1)/ga))))*en*sqrt(Tr);
disp(['Suction to motive gas area ratio :',num2str(Asm)]); 

%% Calculating Throat to motive gas area ratio
Atm = (((1/Pix)^(1/ga))*sqrt(1 - (1/Pix)^((ga-1)/ga)))/((2/(ga+1))*sqrt(1 - (2/(ga+1))));
disp(['Throat to motive gas area ratio :',num2str(Atm)]);

%% Calculating throat to constant area duct ratio 
P23 = P1b_l/Po1b;
At2 = (P3o/Pr)*(1/((1+en)*(1 + en*Tr)))*((((P23)^(1/ga))*sqrt(1 - (P23)^((ga-1)/ga)))/((2/(ga+1))^(1/(ga-1))*sqrt(1 - 2/(ga+1))));
disp(['Throat to constant duct area ratio :',num2str(At2)]);

%% Calulating throat area at specified conditions
Pi = (Pr*Po)*133.3;

At = (mf*sqrt(To))/(0.0404*Pi*Cd);
At = At*(10)^6
disp(['Throat area at specified conditions :',num2str(At)]);

%% Calculating motive gas area
As = (3.14*((Ds)^2))/4; 
Am = As/Asm;
disp(['motive gas area :',num2str(Am)]); 

%% Calculating area of constant duct
At1 = Atm*Am;
disp(['Throat area :',num2str(At1)]);

A2 = At1/At2;      
disp(['constant duct area :',num2str(A2)]);

%% Calculating motive gas diameter
Dm = sqrt((4*Am)/3.14);
disp(['motive gas diameter :',num2str(Dm)]);

%% Calculating constant duct diameter
Dc = sqrt((4*A2)/3.14);
disp(['constant duct diameter :',num2str(Dc)]);

%% Calculating throat diameter
Dt = sqrt((4*At1)/3.14);
disp(['throat diameter :',num2str(Dt)]);

%% Calculating diameter of outer tube
Aw = As + Am;
Dw = sqrt((Aw*4)/3.14);
disp(['outer tube diameter :',num2str(Dw)]);

%% Calculating convergence angle
thd = 2*(atand((sqrt(Aw/A2) - 1)/(L1 + L2)));
disp(['convergence angle :',num2str(thd)]);                fprintf('\n');

disp('Stage 2: ');        fprintf('\n');

%% Writing given conditions 
re = 18;        ga = 1.4;    mf = 4;                       Ti = 300;          L1 = 6;
Cd = 0.9;       Ds = 140;    Po = (Po1b*75)/100;           Mb = M1b;          L2 = 8;
en = 0.030;     Pr = 392;    Tr = 1;                   

disp(['Mach number at the inlet Mb is : ',num2str(Mb)]);
disp(['Practically achievable Stagnation pressure at inlet P0 is : ',num2str(Po)]);

%% Calculating static pressure downstream
Pb = Po/((1 + (1/2)*(ga-1)*(Mb^2))^(ga/(ga-1)));
disp(['Static pressure downstream is : ',num2str(Pb)]);

Px =  (Pb*94.42)/100;
pxo = Px/Po;
Pix = Pr/pxo;          
P3o = re*(pxo);

%% Calculating Mach number of Suction and Motive gas at x
Mx1 = sqrt(((Pix)^((ga-1)/ga) - 1)*(2/(ga-1)));
disp(['Motive Gas Mach number at x : ',num2str(Mx1)]);     

Mx2 = sqrt(((1/pxo)^((ga-1)/ga) - 1)*(2/(ga-1)));     
disp(['Suction Gas Mach number at x is : ',num2str(Mx2)]);

%% Calculating Characteristic Mach number for both gases
Mx1_c = sqrt((((ga+1)/2)*(Mx1^2))/(1 + (ga-1)*(Mx1^2)*(1/2)));
disp(['Motive Gas characteristic Mach number at x : ',num2str(Mx1_c)]);     

Mx2_c = sqrt((((ga+1)/2)*(Mx2^2))/(1 + (ga-1)*(Mx2^2)*(1/2)));
disp(['Suction Gas characteristic Mach number at x : ',num2str(Mx2_c)]);
     
%% Calculating characteristic Mach number of mixed stream from individual mach numbers
M1_c = (Mx1_c + en*Mx2_c*sqrt(Tr))/sqrt((1 + en)*(1 + en*Tr));
disp(['Mixed Stream characteristic Mach number : ',num2str(M1_c)]);

%% Calculating Mixed Stream Mach number
M1 = sqrt((2*(M1_c)^2)/((ga+1) - (ga-1)*(M1_c^2)));
disp(['Mixed Stream Mach number : ',num2str(M1)]);

%% Calculating Mach number downstream of Transverse shock 
M1b = sqrt((2/(ga-1) + M1^2)/(((2*ga)/(ga-1))*(M1^2) - 1));
disp(['Mach number downstream of Transverse shock : ',num2str(M1b)]);

%% Calculating Static pressure at the exit
P1b = Px*(((2*ga)/(ga+1))*(M1^2) - (ga-1)/(ga+1));
disp(['Static pressure downstream of Transverse shock  : ',num2str(P1b)]);        

P1b1 = (P1b*84)/100;
disp(['Static pressure after considering frictional losses  : ',num2str(P1b1)]);        

%% Calculating Stagnation pressure downstream of Transverse shock
Po1b = P1b*((1 + (ga-1)*(1/2)*(M1b^2))^(ga/(ga-1)));
disp(['Stagnation pressure downstream of Transverse shock  : ',num2str(Po1b)]);    

%%Calculating suction to motive gas area ratio    K
Asm = (Pr)*((((1/Pix)^(1/ga))*sqrt(1 - (1/Pix)^((ga-1)/ga)))/((pxo)^(1/ga)*sqrt(1 - (pxo)^((ga-1)/ga))))*en*sqrt(Tr);
disp(['Suction to motive gas area ratio :',num2str(Asm)]);     

%% Calculating Throat to motive gas area ratio
Atm = (((1/Pix)^(1/ga))*sqrt(1 - (1/Pix)^((ga-1)/ga)))/((2/(ga+1))*sqrt(1 - (2/(ga+1))));
disp(['Throat to motive gas area ratio :',num2str(Atm)]);     

%% Calculating throat to constant area duct ratio 
P23 = P1b/Po1b;
At2 = (P3o/Pr)*(1/((1+en)*(1 + en*Tr)))*((((P23)^(1/ga))*(1 - (P23)^((ga-1)/ga))^(1/2))/((2/(ga+1))^(1/(ga-1))*(1 - 2/(ga+1))^(1/2)));
disp(['Throat to constant duct area ratio :',num2str(At2)]);     

%% Calulating throat area at specified conditions
Pi = (Pr*Po)*133.3;

At1 = (mf*sqrt(Ti))/(0.0404*Pi*Cd);
At1 = At1*(10)^6;
disp(['Throat area at specified conditions :',num2str(At1)]);     

%% Calculating motive gas area
As = (3.14*((Ds)^2))/4;
Am = As/Asm;
disp(['motive gas area :',num2str(Am)]);     

%% Calculating area of constant duct
At = Atm*Am;
disp(['throat area :',num2str(At)]);     

A2 = At1/At2;
disp(['constant duct area :',num2str(A2)]);

%% Calculating motive gas diameter
Dm = sqrt((4*Am)/3.14);
disp(['motive gas diameter :',num2str(Dm)]);

%% Calculating constant duct diameter
Dc = sqrt((4*A2)/3.14);
disp(['constant duct diameter :',num2str(Dc)]);

%% Calculating throat diameter
Dt = sqrt((4*At1)/3.14);
disp(['throat diameter :',num2str(Dt)]);

%%Calculating diameter of outer tube
Aw = Am + As;
Dw = sqrt((Aw*4)/3.14);
disp(['outer tube diameter :',num2str(Dw)]);

%%Calculating convergence angle
thd = 2*(atand((sqrt(Aw/A2) - 1)/(L1 + L2)));
disp(['convergence angle theta : ',num2str(thd)]);

