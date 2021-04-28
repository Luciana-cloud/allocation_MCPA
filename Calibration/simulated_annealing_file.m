
%% Simulated Annealing %%

% Example file for simulated annealing calibration %

ObjectiveFunction = @hierarchy_calibration_file;

load('model_parameters_hierarchy.mat')  % example file


ub = [3 10 6 log10(5) log10(0.1) 0.9 5 5 log10(0.1) -9 3 ...
      6 log10(0.1) 0.9 5 5 log10(0.1) -9 3 ...
      3 10 6 log10(5) 0.9 5 5 -9 3 ...
      6 0.9 5 5 -9 3 5.5];              % upper boundary
        
lb = [-3 1 -10 log10(0.15) log10(1e-3) 0.1 -5 -5 log10(1e-5) -13 1.1 ...
      -10 log10(1e-3) 0.1 -5 -5 log10(1e-5) -13 1.1 ...
      -3 1 -10 log10(0.15) 0.1 -5 -5 -13 1.1 ...
      -10 0.1 -5 -5 -13 1.1 3] ;        % lower boundary

p0 = p;

options = optimoptions(@simulannealbnd,'MaxFunEvals',75000);
[p,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,p0,lb,ub,options);%,options);

save('final_calibration.mat','p')