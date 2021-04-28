
function [SSE] = hierarchy_calibration_file(p)

% calibration file including the sum of squared errors (SSE) as objective
% function to be minimized by the simulated annealing algorithm. Due to the
% different order of magnitudes of the data (mineralization, transcripts
% and genes), I took the genes and transcript data in log-form. 

    %% Calling Data %%

    load('mineralization.txt')
    load('total_dna.txt')
    load('total_rna.txt')
    load('residual.txt')
    load('C14_Cmin.txt')

    %% Calibration %%
    
    s = 6; % (20, 1.8, 20°)
    
    [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035            = (mineralization(61:72,5));
    mineral_2035std         = (mineralization(61:72,6));
    dna_2035                = log(total_dna(27:29,5));
    dna_2035std             = log(total_dna(27:29,6));
    rna_2035                = log(total_rna(27:29,5));
    rna_2035std             = log(total_rna(27:29,6));
    residual_2035           = (residual(26:29,5));
    residual_2035std        = (residual(26:29,6));
    C14_CminI               = C14_Cmin(16:18,5)/150;
    C14_CminI_std           = C14_Cmin(16:18,6)/150;
%     
    one                     = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*(log(output(13:15))-rna_2035)./rna_2035std,...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std)); 
                                  
%%
    s = 8;  % (20, 3.5, 20°)

    [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035         = mineralization(85:end,5);
    mineral_2035std      = mineralization(85:end,6);
    dna_2035             = log(total_dna(37:39,5));
    dna_2035std          = log(total_dna(37:39,6));
    rna_2035             = log(total_rna(37:39,5));
    rna_2035std          = log(total_rna(37:39,6));
    residual_2035        = residual(36:39,5);
    residual_2035std     = residual(36:39,6);
    C14_CminI            = C14_Cmin(22:24,5)/150;
    C14_CminI_std        = C14_Cmin(22:24,6)/150;

    two                     = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:15))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std));

%%

    s = 4;  % (20, 3.5, 10°)

    [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035         = mineralization(37:48,5);
    mineral_2035std      = mineralization(37:48,6);
    dna_2035             = log(total_dna(17:19,5));
    dna_2035std          = log(total_dna(17:19,6));
    rna_2035             = log(total_rna(17:19,5));
    rna_2035std          = log(total_rna(17:19,6));
    residual_2035        = residual(16:19,5);
    residual_2035std     = residual(16:19,6);
    C14_CminI            = C14_Cmin(22:24,5)/150;
    C14_CminI_std        = C14_Cmin(22:24,6)/150;

    three                     = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:15))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std));

%%
    s = 2;  % (20, 1.8, 10°)

    [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035       = mineralization(13:24,5);
    mineral_2035std    = mineralization(13:24,6);
    dna_2035           = log(total_dna(7:9,5));
    dna_2035std        = log(total_dna(7:9,6));
    rna_2035           = log(total_rna(7:10,5));
    rna_2035std        = log(total_rna(7:10,6));
    residual_2035      = residual(6:9,5);
    residual_2035std   = residual(6:9,6);
    C14_CminI          = C14_Cmin(4:6,5)/150;
    C14_CminI_std      = C14_Cmin(4:6,6)/150;
% %     
    four           = vertcat(1*((output(1:12)-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:16))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*((output(22:25)-residual_2035)./residual_2035std)); 

%%

    s = 5;  % (1, 1.8, 20°)
    [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035        = (mineralization(49:60,5));
    mineral_2035std     = (mineralization(49:60,6));
    dna_2035            = log(total_dna(22:24,5));
    dna_2035std         = log(total_dna(22:24,6));
    rna_2035            = log(total_rna(22:24,5));
    rna_2035std         = log(total_rna(22:24,6));
    residual_2035       = vertcat(residual(21:24,5));
    residual_2035std    = vertcat(residual(21:24,6));
    C14_CminI           = C14_Cmin(13:15,5)/150;
    C14_CminI_std       = C14_Cmin(13:15,6)/150;
% % % % % % 
    five                = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:15))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std)); 
                                  
%%
    s = 7;  % (1, 3.5, 20°)
   [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035        = (mineralization(73:84,5));
    mineral_2035std     = (mineralization(73:84,6));
    dna_2035            = log(total_dna(32:34,5));
    dna_2035std         = log(total_dna(32:34,6));
    rna_2035            = log(total_rna(32:34,5));
    rna_2035std         = log(total_rna(32:34,6));
    residual_2035       = (residual(31:34,5));
    residual_2035std    = (residual(31:34,6));
    C14_CminI           = C14_Cmin(19:21,5)/150;
    C14_CminI_std       = C14_Cmin(19:21,6)/150;

    six                = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:15))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std)); 
                                  
%%
    s = 3;  % (1, 3.5, 10°)
   [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035        = (mineralization(25:36,5));
    mineral_2035std     = (mineralization(25:36,6));
    dna_2035            = log(total_dna(12:14,5));
    dna_2035std         = log(total_dna(12:14,6));
    rna_2035            = log(total_rna(12:15,5));
    rna_2035std         = log(total_rna(12:15,6));
    residual_2035       = (residual(11:14,5));
    residual_2035std    = (residual(11:14,6));
    C14_CminI           = C14_Cmin(19:21,5)/150;
    C14_CminI_std       = C14_Cmin(19:21,6)/150;

    seven                = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:16))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std));
                                  
%%
% 
    s = 1;  % (1, 1.8, 10°)

    [output] = hierarchy_calibration_model(p,s);
    
    mineral_2035        = mineralization(1:12,5);
    mineral_2035std     = mineralization(1:12,6);
    dna_2035            = log(total_dna(2:4,5));
    dna_2035std         = log(total_dna(2:4,6));
    rna_2035            = log(total_rna(2:5,5));
    rna_2035std         = log(total_rna(2:5,6));
    residual_2035       = residual(1:4,5);
    residual_2035std    = residual(1:4,6);
    C14_CminI           = C14_Cmin(1:3,5)/150;
    C14_CminI_std       = C14_Cmin(1:3,6)/150;
% % 
    eight                = vertcat(1*(((output(1:12))-mineral_2035)./mineral_2035std),...
                                      1*((log(output(13:16))-rna_2035)./rna_2035std),...
                                      1*((log(output(18:20))-dna_2035)./dna_2035std),...
                                      1*(((output(22:25))-residual_2035)./residual_2035std)); 
                                  
%%

SSES = (vertcat(one,two,three,four,five,six,seven,eight)).^2;

subplot(2,4,1);
plot(one.^2)
title('20, 1.8, 20°')
subplot(2,4,2);
plot(two.^2)
title('20, 3.5, 20°')
subplot(2,4,3);
plot(three.^2)
title('20, 3.5, 10°')
subplot(2,4,4);
plot(four.^2)
title('20, 1.8, 10°')
subplot(2,4,5);
plot(five.^2)
title('1, 1.8, 20°')
subplot(2,4,6);
plot(six.^2)
title('1, 3.5, 20°')
subplot(2,4,7);
plot(seven.^2)
title('1, 3.5, 10°')
subplot(2,4,8);
plot(eight.^2)
title('1, 1.8, 10°')
SSE  = sum(SSES)
end