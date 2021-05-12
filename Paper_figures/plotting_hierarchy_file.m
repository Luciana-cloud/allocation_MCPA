function [mineralA,transcriptsA,DNA_totA,P_total,C_14_min,ty,C_cum,CUE_E,CUE_M,CUE_T,CUE_ET,CUE_EM,CUE_EG] = plotting_hierarchy_file(p,s)

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
    q(5)   = 10;         % Temperature in celsius of the experiment
    q(6)   = -0.006;     % Current matrix potential
    q(8)   = 0.015;
    C_tot  = 9e-04/20;   % Total initial MCPA conc. mmol g-1 soil
    f_1    = 10^p(27);   % Conversion factor [carbon/gene]
    n_H    = p(21);      % Hill coeficient  [1]
    K_G    = 10^p(22);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(20);   % *  Conversion factor [transcripts/gene]    
    mu_max = 10^p(23);   % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(5);    % ** Decay rate of active?DNA?_tfdA [d-1]
    Y_P    = p(24);      %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(25);   %    Maintenance term 
    K_M    = 10^p(26);   %    Michaelis menten constant (mmol C cm-3)
    a_CO2  = 10^p(9);    % ** Rate of NER decay
    Q10    = p(28);      %    Temperature function constant
elseif s == 2
    q(5)   = 10;
    q(6)   = -0.006;
    q(8)   = 0.015/20;
    C_tot  = 9e-04;
    f_1    = 10^p(10);  % Conversion factor [carbon/gene]
    n_H    = p(2);      % Hill coeficient  [1]
    K_G    = 10^p(3);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(1);   % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(4);   % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(5);   % ** Decay rate of active?DNA?_tfdA [d-1]
    Y_P    = p(6);      %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(7);   %    Maintenance term  
    K_M    = 10^p(8);   %    Michaelis menten constant (mmol C cm-3) 
    a_CO2  = 10^p(9);   % ** Rate of NER decay 
    Q10    = p(11);     %    Temperature function constant
elseif s == 3
    q(5)   = 10;
    q(6)   = -0.3;
    q(8)   = 0.015;
    C_tot  = 9e-04/20;
    f_1    = 10^p(33);  % Conversion factor [carbon/gene]
    n_H    = p(21);     % Hill coeficient  [1]
    K_G    = 10^p(29);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(20);  % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(23);  % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(13);  % ** Decay rate of active?DNA?_tfdA [d-1]
    Y_P    = p(30);     %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(31);  %    Maintenance term 
    K_M    = 10^p(32);  %    Michaelis menten constant (mmol C cm-3)
    a_CO2  = 10^p(17);  % ** Rate of NER decay
    Q10    = p(34);     %    Temperature function constant
elseif s == 4
    q(5)   = 10;
    q(6)   = -0.3;
    q(8)   = 0.015/20;
    C_tot  = 9e-04;
    f_1    = 10^p(18);  % Conversion factor [carbon/gene]
    n_H    = p(2);      % Hill coeficient  [1]
    K_G    = 10^p(12);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(1);   % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(4);   % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(13);  % ** Decay rate of active?DNA?_tfdA [d-1] 
    Y_P    = p(14);     %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(15);  %    Maintenance term 
    K_M    = 10^p(16);  %    Michaelis menten constant (mmol C cm-3)
    a_CO2  = 10^p(17);  % ** Rate of NER decay
    Q10    = p(19);     %    Temperature function constant
elseif s == 5
    q(5)   = 20;
    q(6)   = -0.006;
    q(8)   = 0.015;
    C_tot  = 9e-04/20;
    f_1    = 10^p(27);   % Conversion factor [carbon/gene]
    n_H    = p(21);      % Hill coeficient  [1]
    K_G    = 10^p(22);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(20);   % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(23);   % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(5);    % ** Decay rate of active?DNA?_tfdA [d-1]
    Y_P    = p(24);      %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(25);   %    Maintenance term 
    K_M    = 10^p(26);   %    Michaelis menten constant (mmol C cm-3)
    a_CO2  = 10^p(9);    % ** Rate of NER decay
    Q10    = p(28);      %    Temperature function constant
elseif s == 6
    q(5)   = 20;
    q(6)   = -0.006;
    q(8)   = 0.015/20;
    C_tot  = 9e-04;
    f_1    = 10^p(10);  % Conversion factor [carbon/gene]
    n_H    = p(2);      % Hill coeficient  [1]
    K_G    = 10^p(3);   % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(1);   % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(4);   % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(5);   % ** Decay rate of active?DNA?_tfdA [d-1]
    Y_P    = p(6);      %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(7);   %    Maintenance term  
    K_M    = 10^p(8);   %    Michaelis menten constant (mmol C cm-3) 
    a_CO2  = 10^p(9);   % ** Rate of NER decay
    Q10    = p(11);     %    Temperature function constant
elseif s == 7
    q(5)   = 20;
    q(6)   = -0.3;
    q(8)   = 0.015;
    C_tot  = 9e-04/20;
    f_1    = 10^p(33);  % Conversion factor [carbon/gene]
    n_H    = p(21);     % Hill coeficient  [1]
    K_G    = 10^p(29);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(20);  % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(23);  % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(13);  % ** Decay rate of active?DNA?_tfdA [d-1]
    Y_P    = p(30);     %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(31);  %    Maintenance term 
    K_M    = 10^p(32);  %    Michaelis menten constant (mmol C cm-3)
    a_CO2  = 10^p(17);  % ** Rate of NER decay
    Q10    = p(34);     %    Temperature function constant
else 
    q(5)   = 20;
    q(6)   = -0.3;
    q(8)   = 0.015/20;
    C_tot  = 9e-04;
    f_1    = 10^p(18);  % Conversion factor [carbon/gene]
    n_H    = p(2);      % Hill coeficient  [1]
    K_G    = 10^p(12);  % Activation coefficient of MCPA for gene expression [mmol C cm-3]
    f_T    = 10^p(1);   % *  Conversion factor [transcripts/gene]
    mu_max = 10^p(4);   % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
    a_a    = 10^p(13);  % ** Decay rate of active?DNA?_tfdA [d-1] 
    Y_P    = p(14);     %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
    m      = 10^p(15);  %    Maintenance term 
    K_M    = 10^p(16);  %    Michaelis menten constant (mmol C cm-3)
    a_CO2  = 10^p(17);  % ** Rate of NER decay
    Q10    = p(19);     %    Temperature function constant
end

tspan   = [0 30];

%% Options for the solver %%

abstol = 1e-9;
reltol = 1e-7;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:6); % 'Events',@stopevent);%,'NonNegative',1:7);
  
K_F   = q(11);
n_F   = q(12);
C_L0  = 0;
th_V  = q(3);      % Average volumetric soil water content (cm^3 cm^-3)
rho_B = q(4);      % Bulk density of soils (g cm^-3)

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
c(6) = 0;                       % Initial cumulative pesticide-derived C in the biomass

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
    NER     = cu(:,3);  % [mmol C g-1]
    P       = cu(:,4);  % [mmol C cm-3]
    CO2     = cu(:,5);  % [mmol C g-1]
    R       = cu(:,6);  % [mmol C g-1]
    
        Initial     = C_tot*q(8);    % [mmol C g-1]
        alpha       = q(9);
        RNA         = f_T.*(P.^(n_H)).*((K_G^(n_H)+P.^(n_H)).^(-1)) + alpha; % [transcripts/gene]
        C_14_MCPA   = 6.6e-7;   % [mmol C g-1]
        
%% Compute full time series of measured parameters from model output %%

    mineralA (:,1)       = CO2.*100/Initial;                % [%]
    transcriptsA (:,1)   = RNA.*(DNA_a*f_1^(-1));           % [transcripts g-1]
    DNA_totA (:,1)       = (DNA_a)*f_1^(-1);                % [gene g-1]
    P_total (:,1)        = P.*th_V/rho_B+K_F.*P.^n_F;       % [mg C g-1]
    C_14_min(:,1)        = DNA_14.*100/Initial;             % [%]
    C_cum(:,1)           = R.*100/Initial;                  % Cumulative pesticide-derived C in the biomass [%]
    
%% CUE %%

% Constant %

n_M     = q(1);    % n_M - Conversion factor from mmol MCPA to mmol C  [mmol C mmol MCPA^-1]
th_V    = q(3);    % th_V - Average volumetric soil water content (cm^3 cm^-3)
rho_B   = q(4);    % rho_B - Bulk density of soils (g cm^-3)
T       = q(5);    % Temperature in celsius of the experiment
psi     = q(6);    % Current matrix potential
psi_fc  = q(7);    % Field capacity matrix potential 
a_14    = q(8);    % the fraction of 14C in the applied MCPA
alphaA  = q(9);
K_F_P   = q(11);   % Freundlich coeff of MCPA sorption isotherm (mmol MCPA g^-1 soil/(mmol MCPA cm^-3)^nF_MCPA)
n_F_P   = q(12);   % Freundlich exponent of MCPA sorption isotherm (1)

% Rates %

AR              = Q10^((T-10)/10);
mu_P            = mu_max.*AR.*(P.^(n_H).*(K_G^(n_H)+P.^(n_H)).^(-1)+alphaA/f_T)...
                 .*(P.*(P+K_M).^(-1));
r_T             = m;                    % Exogeneous maintenance
r_Ex            = m*(P.*(P+K_M).^(-1)); % Endogenuos maintenance
r_En            = r_T - r_Ex;             
r_respiration   = a_14.*DNA_a.*mu_P.*(1-Y_P).*Y_P^(-1);
r_mainte_Pd     = r_Ex.*DNA_a*a_14*AR;
r_mainte_Bd     = r_En.*DNA_14*AR;
r_growth        = DNA_a.*(mu_P).*a_14;
r_decay         = DNA_14*a_a*AR;
r_NER           = NER.*a_CO2.*AR; 

% CUE %

CUE_E  = DNA_14./(DNA_14+CO2); % Identical to CUE calculation based on experimental data

CUE_M  = 1 - ((r_respiration+r_mainte_Pd+r_mainte_Bd)./(r_growth+r_respiration+...
               r_mainte_Pd)); % Considering Maintenance

CUE_T  = 1 - ((r_respiration+r_mainte_Pd+r_mainte_Bd+r_NER)./(r_growth+...
               r_respiration+r_mainte_Pd)); % Considering microbial turnover

CUE_ET = (r_growth-r_decay)./(r_growth+r_respiration+r_mainte_Pd+r_mainte_Bd+r_NER); % Similar to CUE calculation based on experimental data, but based on rates, considering maintenance and turnover

CUE_EM = (r_growth-r_decay)./(r_growth+r_respiration+r_mainte_Pd+r_mainte_Bd); % Similar to CUE calculation based on experimental data, without microbial turnover, but with maintenance

CUE_EG = (r_growth-r_decay)./(r_growth+r_respiration); % Similar to CUE calculation based on experimental data, only growth

end