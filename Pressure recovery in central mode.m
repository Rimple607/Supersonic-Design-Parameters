% to study the central injection mode of pressure recovery

%% Stage 1 
%% Input parameters
ms = 3;       Px = 4;       Tr = 1.33;      G = 1.4;        Cd = 0.9;    mc = 1.5;
w = 0.025;    Pr = 300;     To_s = 400;     re = 18;        Ds = 75;     Pc = 3;

% Calculating parameters
%% Motive gas flow rate
mm = ms/w;
disp(['Motive gas flow rate mm : ',num2str(mm)]);

%% Calculating Mach number
M_c = fzero(@(M) flowisentropic(M), [1; 4]);
disp(['Mach number of flow from cavity M_c : ',num2str(M_c)]);
  
%% Calculating stagnation pressure
Po = Pc*(1 + (G-1)*(mc^2)*(1/2))^(G/(G-1));
disp(['Stagnation pressure of flow Po : ',num2str(Po)]);

%% Calculating Mach number downstream of shock
Mb = sqrt((2/(G-1) + M_c^2)/( ((2*G)/(G-1))*(M_c^2) - 1));
disp(['Mach number downstream of shock Mb : ',num2str(Mb)]);

%% Stagnation pressure downstream of shock
Pob = Po*(1 + ((2*G)/(G+1))*(M_c^2 - 1))*((1 + (1/2)*(G-1)*(Mb^2))/(1 + (1/2)*(G-1)*(M_c^2)))^(G/(G-1));
disp(['stagnation pressure downstream of shock Pob : ',num2str(Pob)]);

%% Calculating static pressure downstream
Psb = Pob/((1 + (1/2)*(G-1)*(Mb^2))^(G/(G-1)));
disp(['Static pressure downstream of shock Psb : ',num2str(Psb)]);

%% Motive gas stagnation pressure
Po_m = Pr*Pob;
disp(['Motive gas stagnation pressure Po_m : ',num2str(Po_m)]);

%% Motive gas Mach number
Pix = Po_m/Px;         disp(['Pix : ',num2str(Pix)]);
Pxi = 1/Pix;           disp(['Pxi : ',num2str(Pxi)]);

Mx_m = sqrt((Pix^((G-1)/G) - 1 )*(2/(G-1)));
disp(['Motive gas Mach number Mx_m : ',num2str(Mx_m)]);

Mxm_c = sqrt((((G+1)/2)*(Mx_m^2))/(1 + ((G-1)/2)*(Mx_m^2)));
disp(['Motive gas characteristic Mach number Mxm_c : ',num2str(Mxm_c)]);

%% Suction gas Mach number
Pox = Pob/Px;         disp(['Pox : ',num2str(Pox)]);
Pxo = 1/Pox;         disp(['Pxo : ',num2str(Pxo)]);

Mx_s = sqrt((Pox^((G-1)/G) - 1 )*(2/(G-1)));
disp(['Suction gas Mach number Mx_s : ',num2str(Mx_s)]);

Mxs_c = sqrt((((G+1)/2)*(Mx_s^2))/(1 + ((G-1)/2)*(Mx_s^2)));
disp(['Suction gas characteristic Mach number Mxs_c : ',num2str(Mxs_c)]);

%% Mixed stream Mach number
M1_c = (Mxm_c + w*Mxs_c*sqrt(Tr))/sqrt((1+w)*(1+(w*Tr)));
disp(['Mixed stream characteristic Mach number : ',num2str(M1_c)]);

M1 =  sqrt((2*(M1_c^2))/((G+1) - (G-1)*(M1_c^2)));
disp(['Mixed stream Mach number  M1 : ',num2str(M1)]);

%% Mach number downstream shock 
M1b = sqrt((2/(G-1) + M1^2)/(((2*G)/(G-1))*(M1^2) - 1));
disp(['Mach number after the shock M1b : ',num2str(M1b)]);

%% Static pressure after the shock
Pb = Px*(((2*G)/(G+1))*(M1^2) - (G-1)/(G+1));
disp(['Static pressure after the shock Pb : ',num2str(Pb)]);
Pb = (Pb*84)/100;
disp(['after considering frictional losses static pressure Pb : ',num2str(Pb)]);

%% Stagnation pressure downstream 
Pbo = Pb*((1 + (G-1)*(1/2)*(M1b^2))^(G/(G-1)));
disp(['Stagnation pressure after the shock Pbo : ',num2str(Pbo)]);

%% Suction to motive gas area ratio
Asm = Pr*w*sqrt(Tr)*((Pxi^(1/G))*sqrt(1 - Pxi^((G-1)/G)))/((Pxo^(1/G))*sqrt(1 - Pxo^((G-1)/G)));
disp(['Suction to Motive gas area ratio Asm : ',num2str(Asm)]);

%% Motive throat to ejector entry area ratio
Atm = (((Pxi)^(1/G))*sqrt(1 - (1/Pix)^((G-1)/G)))/((2/(G+1))*sqrt(1 - (2/(G+1))));
disp(['Motive throat to Ejector entry area ratio Atm : ',num2str(Atm)]);

%% Motive throat to constant area duct ratio
P23 = Pb/Pbo;
P3o = re/Pox;

At2 = (P3o/Pr)*(1/((1+w)*(1 + w*Tr)))*((((P23)^(1/G))*sqrt(1 - (P23)^((G-1)/G)))/((2/(G+1))^(1/(G-1))*sqrt(1 - 2/(G+1))));
disp(['Motive throat to constant area duct ratio At2 : ',num2str(At2)]);

%% Suction gas inlet area
As = 3.14*(Ds^2)/4;
disp(['Suction gas inlet area As : ',num2str(As)]);

%% Motive gas inlet area 
Am = As/Asm;
disp(['Motive gas inlet area Am : ',num2str(Am)]);

Dm = sqrt(Am*4/3.14);
disp(['Motive gas diameter Dm : ',num2str(Dm)]);

%% Calulating throat area 
Pi = Po_m*133.3;
To_m = To_s/Tr;         

mm = mm*10^(-3);
At_c = (mm*sqrt(To_m))/(0.0404*Pi*Cd);
At_c = At_c*(10)^6;
disp(['Throat area at specified conditions At_c : ',num2str(At_c)]);

At = Atm*Am;
disp(['Throat area At : ',num2str(At)]);

Dt = sqrt(At*4/3.14);
disp(['Throat diameter Dt : ',num2str(Dt)]);

%% Calculating constant duct area
A2 = At/At2;
disp(['Constant duct area A2 : ',num2str(A2)]); 

D2 = sqrt(A2*4/3.14);
disp(['constant duct diameter : ',num2str(D2)]);

%% Convergence angle
Ax = As + Am;
Dx = sqrt(Ax*4/3.14);         disp(['Outer tube diameter : ',num2str(Dx)]);
L1 = 6*D2;                    disp(['Length of constant pressure mixing chamber L1 : ',num2str(L1)]);
L2 = 6*D2;                    disp(['Length of constant area duct L2 : ',num2str(L2)]);
theta = 2*atand(((Dx/D2) - 1)/(L1 + L2));           disp(['Convergence angle theta : ',num2str(theta)]);


%% Stage 2
%% Input parametres
ms = 3;       Pob = (Pbo*75.6)/100;      Tr = 1.33;      G = 1.4;        Cd = 0.9;      mm = 4;
w = 0.030;    Pr = 392;                  To_s = 400;     re = 18;         Ds = 140;

%% Calculating static pressure downstream
Psb = Pob/((1 + (1/2)*(G-1)*(Mb^2))^(G/(G-1)));
disp(['Static pressure downstream of shock Psb : ',num2str(Psb)]);

Px = (Psb*99)/100
Pob

%% Motive gas stagnation pressure
Po_m = Pr*Pob;
disp(['Motive gas stagnation pressure Po_m : ',num2str(Po_m)]);

%% Motive gas Mach number
Pix = Po_m/Px;         disp(['Pix : ',num2str(Pix)]);
Pxi = 1/Pix;           disp(['Pxi : ',num2str(Pxi)]);

Mx_m = sqrt((Pix^((G-1)/G) - 1 )*(2/(G-1)));
disp(['Motive gas Mach number Mx_m : ',num2str(Mx_m)]);

Mxm_c = sqrt((((G+1)/2)*(Mx_m^2))/(1 + ((G-1)/2)*(Mx_m^2)));
disp(['Motive gas characteristic Mach number Mxm_c : ',num2str(Mxm_c)]);     

%% Suction gas Mach number
Pox = Pob/Px;         disp(['Pox : ',num2str(Pox)]);
Pxo = 1/Pox;         disp(['Pxo : ',num2str(Pxo)]);

Mx_s = sqrt((Pox^((G-1)/G) - 1 )*(2/(G-1)));
disp(['Suction gas Mach number Mx_s : ',num2str(Mx_s)]);

Mxs_c = sqrt((((G+1)/2)*(Mx_s^2))/(1 + ((G-1)/2)*(Mx_s^2)));
disp(['Suction gas characteristic Mach number Mxs_c : ',num2str(Mxs_c)]);

%% Mixed stream Mach number
M1_c = (Mxm_c + w*Mxs_c*sqrt(Tr))/sqrt((1+w)*(1+(w*Tr)));
disp(['Mixed stream characteristic Mach number : ',num2str(M1_c)]);

M1 =  sqrt((2*(M1_c^2))/((G+1) - (G-1)*(M1_c^2)));
disp(['Mixed stream Mach number  M1 : ',num2str(M1)]);

%% Mach number downstream shock 
M1b = sqrt((2/(G-1) + M1^2)/(((2*G)/(G-1))*(M1^2) - 1));
disp(['Mach number after the shock M1b : ',num2str(M1b)]);

%% Static pressure after the shock
Pb = Px*(((2*G)/(G+1))*(M1^2) - (G-1)/(G+1));
disp(['Static pressure after the shock Pb : ',num2str(Pb)]);
Pb = (Pb*84)/100;
disp(['after considering frictional losses static pressure Pb : ',num2str(Pb)]);

%% Stagnation pressure downstream 
Pbo = Pb*((1 + (G-1)*(1/2)*(M1b^2))^(G/(G-1)));
disp(['Stagnation pressure after the shock Pbo : ',num2str(Pbo)]);

%% Suction to motive gas area ratio
Asm = Pr*w*sqrt(Tr)*((Pxi^(1/G))*sqrt(1 - Pxi^((G-1)/G)))/((Pxo^(1/G))*sqrt(1 - Pxo^((G-1)/G)));
disp(['Suction to Motive gas area ratio Asm : ',num2str(Asm)]);

%% Motive throat to ejector entry area ratio
Atm = (((Pxi)^(1/G))*sqrt(1 - (1/Pix)^((G-1)/G)))/((2/(G+1))*sqrt(1 - (2/(G+1))));
disp(['Motive throat to Ejector entry area ratio Atm : ',num2str(Atm)]);

%% Motive throat to constant area duct ratio
P23 = Pb/Pbo;
P3o = re/Pox;

At2 = (P3o/Pr)*(1/((1+w)*(1 + w*Tr)))*((((P23)^(1/G))*sqrt(1 - (P23)^((G-1)/G)))/((2/(G+1))^(1/(G-1))*sqrt(1 - 2/(G+1))));
disp(['Motive throat to constant area duct ratio At2 : ',num2str(At2)]);

%% Suction gas inlet area
As = 3.14*(Ds^2)/4;
disp(['Suction gas inlet area As : ',num2str(As)]);

%% Motive gas inlet area 
Am = As/Asm;
disp(['Motive gas inlet area Am : ',num2str(Am)]);

Dm = sqrt(Am*4/3.14);
disp(['Motive gas diameter Dm : ',num2str(Dm)]);

%% Calulating throat area 
Pi = Po_m*133.3;
To_m = To_s/Tr;         

At_c = (mm*sqrt(To_m))/(0.0404*Pi*Cd);
At_c = At_c*(10^6);
disp(['Throat area at specified conditions At_c : ',num2str(At_c)]);

At = Atm*Am;
disp(['Throat area At : ',num2str(At)]);

Dt = sqrt(At*4/3.14);
disp(['Throat diameter Dt : ',num2str(Dt)]);

%% Calculating constant duct area
A2 = At/At2;
disp(['Constant duct area A2 : ',num2str(A2)]); 

D2 = sqrt(A2*4/3.14);
disp(['constant duct diameter : ',num2str(D2)]);

%% Convergence angle
Ax = As + Am;
Dx = sqrt(Ax*4/3.14);         disp(['Outer tube diameter : ',num2str(Dx)]);
L1 = 6*D2;                    disp(['Length of constant pressure mixing chamber L1 : ',num2str(L1)]);
L2 = 8*D2;                    disp(['Length of constant area duct L2 : ',num2str(L2)]);
theta = 2*atand(((Dx/D2) - 1)/(L1 + L2));           disp(['Convergence angle theta : ',num2str(theta)]);










