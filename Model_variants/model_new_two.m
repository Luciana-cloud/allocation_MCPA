function vf_ = model_new_two(t,x,p1,q)

% Model versions are in core the same, but we performed a hierarchical model 
% calibration assuming different bacterial subpopulation under the two
% different initial concentrations of MCPA applied (* parameters), and possible 
% physiological and morphological bacterial changes under different humidity 
% levels (** parameters). 
% The conditions simulated by this model version are: 
% initial concentration = 20 mg/Kg
% humidity level        = pF 3.5
% temperature level     = both (10 and 20°C)

%% STOCKS

DNA_a      = x(1); % Active DNA [mmol C g-1]  
DNA_14     = x(2); % Active 14C DNA [mmol C g-1]  
NER        = x(3); % NER pool [mmol C g-1]
P          = x(4); % MPCA in solution [mmol C cm-3]

%% PARAMETERS %% * between populations parameters (1 and 20 mg/Kg); 
%               ** between humidity levels parameters (pF 1.8 and 3.5)

% BACTERIAL PARAMETERS
f_T    = 10^p1(1);      % *  Conversion factor [transcripts/gene]
n_H    = p1(2);         % *  Hill coeficient  [1]
K_G    = 10^p1(12);     %    Activation coefficient of MCPA for gene expression [mmol C cm-3]
mu_max = 10^p1(4);      % *  Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
a_a    = 10^p1(13);     % ** Decay rate of active?DNA?_tfdA [d-1] 
Y_P    = p1(14);        %    Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
m      = 10^p1(15);     %    Maintenance term 
K_M    = 10^p1(16);     %    Michaelis menten constant (mmol C cm-3)
a_CO2  = 10^p1(17);     % ** Rate of NER decay
f_1    = 10^p1(18);     %    Conversion factor [carbon/gene]

% FUNCTION PARAMETERS 
Q10      = p1(19);      %    Temperature function constant

%% VALUES OF CONSTANT %%

th_V    = q(3);    % th_V - Average volumetric soil water content (cm^3 cm^-3)
rho_B   = q(4);    % rho_B - Bulk density of soils (g cm^-3)
T       = q(5);    % Temperature in celsius of the experiment
a_14    = q(8);    % Fraction of 14C in the applied MCPA
alphaA  = q(9);    % Basal gene expression
K_F_P   = q(11);   % Freundlich coeff of MCPA sorption isotherm (mmol MCPA g^-1 soil/(mmol MCPA cm^-3)^nF_MCPA)
n_F_P   = q(12);   % Freundlich exponent of MCPA sorption isotherm (1)

%% BIOKINETIC FUNCTIONS %%

% Temperature function
AR = Q10^((T-10)/10);
% Substrate (MCPA)  dependent specific growth rate of bacterial pesticide
% degraders bacteria  [d-1]
mu_P = mu_max*AR*(P^(n_H)*(K_G^(n_H)+P^(n_H))^(-1)+alphaA/f_T)*(P*(P+K_M)^(-1));
% Total Maintenance 
r_T  = m;
% Exogeneous maintenance
r_Ex = m*(P*(P+K_M)^(-1));
% Endogenuos maintenance
r_En = r_T - r_Ex;

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(6,1);

% Active DNA_tfdA  - [mmol C g-1]
vf_(1) = DNA_a*(mu_P) - r_En*DNA_a*AR - (DNA_a)*a_a*AR;

% Active 14C_DNA_tfdA  - [mmol C g-1]
vf_(2) = DNA_a*(mu_P)*a_14 - r_En*DNA_14*AR - DNA_14*a_a*AR;

% NER - [mmol C g-1]
vf_(3) = DNA_14*a_a*AR - NER*a_CO2*AR;

% Solution phase MCPA concentration - mmmol C cm^(-3) soil solution)
vf_(4) = -(mu_P*DNA_a*(th_V^(-1)*rho_B)*(Y_P^(-1)) + r_Ex*DNA_a*(th_V^(-1)*rho_B)*AR)*...
          (1+th_V^(-1)*rho_B*P^(n_F_P-1)*K_F_P*n_F_P)^(-1);

% Total CO2 [mmol C g-1]
vf_(5)= a_14*DNA_a*mu_P*(1-Y_P)*Y_P^(-1) + NER*a_CO2*AR + ...
        r_En*DNA_14*AR + r_Ex*DNA_a*a_14*AR;

% Cumulative pesticide-derived C in the biomass [mmol C g-1]
vf_(6)= DNA_a*mu_P*a_14;    
end