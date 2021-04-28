function [output] = hierarchy_calibration_model(p,s)

%% Define fixed parameters %%

q(1)  = 9;          % n_M - Conversion factor from mmol MCPA to mmol C  [mmol C mmol MCPA^-1]
q(2)  = 0.1;        % n_S  - Switch function parameter
q(3)  = 0.35;       % th_V - Average volumetric soil water content (cm^3 cm^-3)
q(4)  = 1.1;        % rho_B - Bulk density of soils (g cm^-3)
q(7)  = -0.006;     % Field capacity matrix potential 
q(9)  = 1.2e-5;
q(11) = 0.09;       % Freundlich coeff of MCPA sorption isotherm (mmol MCPA g^-1 soil/(mmol MCPA cm^-3)^nF_MCPA)
q(12) = 0.8576;     % Freundlich exponent of MCPA sorption isotherm (1)

if s == 1
    q(5)  = 10;         % Temperature in celsius of the experiment
    q(6)  = -0.006;     % Current matrix potential
    q(8)  = 0.015;
    C_tot = 9e-04/20;   % Total initial MCPA conc. mmol g-1 soil
    f_1   = 10^p(27);   % Conversion factor [carbon/gene]
    n_H   = p(21);      % Hill coeficient  [1]
    K_G   = 10^p(22);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(20);   % *  Conversion factor [transcripts/gene]
elseif s == 2
    q(5)  = 10;
    q(6)  = -0.006;
    q(8)  = 0.015/20;
    C_tot = 9e-04;
    f_1   = 10^p(10);  % Conversion factor [carbon/gene]
    n_H   = p(2);      % Hill coeficient  [1]
    K_G   = 10^p(3);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(1);   % *  Conversion factor [transcripts/gene]
elseif s == 3
    q(5)  = 10;
    q(6)  = -0.3;
    q(8)  = 0.015;
    C_tot = 9e-04/20;
    f_1   = 10^p(33);  % Conversion factor [carbon/gene]
    n_H   = p(21);     % Hill coeficient  [1]
    K_G   = 10^p(29);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(20);  % *  Conversion factor [transcripts/gene]
elseif s == 4
    q(5)  = 10;
    q(6)  = -0.3;
    q(8)  = 0.015/20;
    C_tot = 9e-04;
    f_1   = 10^p(18);  % Conversion factor [carbon/gene]
    n_H   = p(2);      % Hill coeficient  [1]
    K_G   = 10^p(12);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(1);   % *  Conversion factor [transcripts/gene]
elseif s == 5
    q(5)  = 20;
    q(6)  = -0.006;
    q(8)  = 0.015;
    C_tot = 9e-04/20;
    f_1   = 10^p(27);  % Conversion factor [carbon/gene]
    n_H   = p(21);     % Hill coeficient  [1]
    K_G   = 10^p(22);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(20);  % *  Conversion factor [transcripts/gene]
elseif s == 6
    q(5)  = 20;
    q(6)  = -0.006;
    q(8)  = 0.015/20;
    C_tot = 9e-04;
    f_1   = 10^p(10);  % Conversion factor [carbon/gene]
    n_H   = p(2);      % Hill coeficient  [1]
    K_G   = 10^p(3);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(1);   % *  Conversion factor [transcripts/gene]
elseif s == 7
    q(5)  = 20;
    q(6)  = -0.3;
    q(8)  = 0.015;
    C_tot = 9e-04/20;
    f_1   = 10^p(33);  % Conversion factor [carbon/gene]
    n_H   = p(21);     % Hill coeficient  [1]
    K_G   = 10^p(29);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(20);  % *  Conversion factor [transcripts/gene]
else 
    q(5)  = 20;
    q(6)  = -0.3;
    q(8)  = 0.015/20;
    C_tot = 9e-04;
    f_1   = 10^p(18);  % Conversion factor [carbon/gene]
    n_H   = p(2);      % Hill coeficient  [1]
    K_G   = 10^p(12);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T   = 10^p(1);   % *  Conversion factor [transcripts/gene]
end

tspan   = linspace(0,30,31);

%% Options for the solver %%

abstol = 1e-9;
reltol = 1e-7;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:5); % 'Events',@stopevent);%,'NonNegative',1:7);
  
K_F   = q(11);
n_F   = q(12);
C_L0  = 0;
th_V  = q(3);      % Average volumetric soil water content [cm^3 cm^-3]
rho_B = q(4);      % Bulk density of soils [g cm^-3]

fun=@(C_L0) MCPA_init(C_L0,C_tot,K_F,n_F,th_V,rho_B);
        
%  switch off solver progress information on solver progress of fsolve
f_opts = optimoptions('fsolve','Display','none');
        
% calculate initial MCPA in solution phase
C_L0 = fsolve(fun,C_L0,f_opts); % MCPA in solution
        
% DNA %
DNA_T     = 10^p(35)*f_1;    % total initíal gene abundance [gene copies cm-3]
        
%% Set INITIAL CONDITIONS
        
c(1) = DNA_T;                   % Active DNA _tfdA  [mmolC g-1]: DNA^ac
c(2) = 0;                       % Active 14C_DNA _tfdA  [mmolC g-1]: DNA^ac
c(3) = 0;                       % NER pool [mmol C g-1]
c(4) = C_L0;                    % MCPA in solution (mass per volume of water) (mmol C cm^(-3))
c(5) = 0;                       % Total CO2 [mmolC g-1]   

%% Running ode %%
       
 try
     warning off
     tic
     if s == 1
         [ty,cu] = ode15s(@model_new_three,tspan,c,o_opts,p',q);
     elseif s == 2
         [ty,cu] = ode15s(@model_new_one,tspan,c,o_opts,p',q);
     elseif s == 3
         [ty,cu] = ode15s(@model_new_four,tspan,c,o_opts,p',q);
     elseif s == 4
         [ty,cu] = ode15s(@model_new_two,tspan,c,o_opts,p',q);
     elseif s == 5
         [ty,cu] = ode15s(@model_new_three,tspan,c,o_opts,p',q);
     elseif s == 6
         [ty,cu] = ode15s(@model_new_one,tspan,c,o_opts,p',q);
     elseif s == 7
         [ty,cu] = ode15s(@model_new_four,tspan,c,o_opts,p',q);
     else
         [ty,cu] = ode15s(@model_new_two,tspan,c,o_opts,p',q);
     end

 catch ME
     warning off
     ty = length(tspan);
     cu = ones(length(tspan),length(c))*1e+99;     
 end
 
 if length(cu) < length(tspan)
    cu = ones(length(tspan),length(c))*1e+99;
end

if isreal(cu)==0
    cu = ones(length(tspan),length(c))*1e+99;    
end

%% Post-proccessing %%

    DNA_a   = cu(:,1);  % [mmol C g-1]
    DNA_14  = cu(:,2);  % [mmol C g-1]
    P       = cu(:,4);  % [mmol C cm-3]
    CO2     = cu(:,5);  % [mmol C g-1]
    unit    = 1000*200/9; % Unit correction from mmol to mg MCPA
        
        Initial     = C_tot*q(8);    % [mmol C g-1]
        alpha       = q(9);
        RNA         = f_T.*(P.^(n_H)).*((K_G^(n_H)+P.^(n_H)).^(-1)) + alpha; % [transcripts/gene]
        C_14_MCPA   = 6.6e-7;   % [mmol C g-1]
        
%% Compute full time series of measured parameters from model output %%

    mineralA (:,1)       = CO2.*100/Initial;               % [%]
    transcriptsA (:,1)   = RNA.*(DNA_a*f_1^(-1));          % [transcripts g-1]
    DNA_totA (:,1)       = (DNA_a)*f_1^(-1);   % [gene g-1]
    P_total (:,1)        = (P.*th_V/rho_B+K_F.*P.^n_F)*unit;      % [mg C g-1]
    C_14_min(:,1)        = (DNA_14)/C_14_MCPA;
    
%% Final output file %%

output = vertcat(mineralA([3 6 8 11 14 16 18 20 22 24 27 29]),transcriptsA([7 11 16 27]),...
                DNA_totA([1 7 11 16 27]),...
                P_total([1 7 11 17 27]),...
                C_14_min([6 16 29]));
end