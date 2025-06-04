function [fit_func, C_fit, gof, beta_fit] = createFit_Cvals(z_fit, delta_fit, C1_forced, C2_forced, C3_forced)
% Input
% z             - [m] Depth vector
% delta         - [m] Wall deflection
% C1_forced     - [-] Forced cantilever coefficient; can be set to a value like 0.5 or nan
% C2_forced     - [-] Forced parabolic coefficient; can be set to a value like 0.5 or nan
% C3_forced     - [-] Forced kick-in coefficient; can be set to a value like 0.5 or nan

% Output
% fit_func      - [m] Fitting function based on C1, C2, C3 called by: fit_func(x)
% C_fit         - [-] Wall deflection coefficients C1, C2, C3

% Fitted wall height according to given data
Hw_fit  = max(z_fit) - min(z_fit)   ; % [m] 
    
% Volume loss along fitted wall
VLW_fit     = trapz(z_fit, delta_fit)    ; % [m^2] VLW_fit

% Average wall deflection along fitted wall
delta_x_avg_fit     = VLW_fit/Hw_fit    ; % [m] 

% Average wall deflection for the fitted wall
% beta_fit    = mean(delta_fit)/Hw_fit        ; % [-]
beta_fit    = delta_x_avg_fit/Hw_fit        ; % [-]

count = 0;
if isnan(C1_forced)     ; C1_txt = 'C1';                
else                    ; C1_txt = num2str(C1_forced);  count = count + 1;
end

if isnan(C2_forced)     ; C2_txt = 'C2';                
else                    ; C2_txt = num2str(C2_forced);  count = count + 1; 
end

if isnan(C3_forced)     ; C3_txt = ['(1-' C1_txt '-' C2_txt ')']; 
else                    ; C3_txt = num2str(C3_forced);  count = count + 1;
end

if count == 2
    idC = find(isnan([C1_forced, C2_forced, C3_forced]));
    if      idC == 1; C1_txt = ['(1-' C2_txt '-' C3_txt ')'];
    elseif  idC == 2; C2_txt = ['(1-' C1_txt '-' C3_txt ')'];
    elseif  idC == 3; C3_txt = ['(1-' C1_txt '-' C2_txt ')'];
    end
end

% Define fittype and options according to forced C1 or not
if count < 2
    ft = fittype([ C1_txt '*' num2str(beta_fit) '* ( 2 * ' num2str(Hw_fit) '                                  - 2 * (x -' num2str(z_fit(1)) '))' ...
               '+' C2_txt '*' num2str(beta_fit) '* (-6 / ' num2str(Hw_fit) ' * ((x -' num2str(z_fit(1)) '))^2 + 6 * (x -' num2str(z_fit(1)) '))' ...
               '+' C3_txt '*' num2str(beta_fit) '* ( 2 * ((x -' num2str(z_fit(1)) ') ))'                                                      ], ...
               'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    
    if count == 1;  opts.StartPoint = 1     ;
    else        ;   opts.StartPoint = [0 1] ;
    end

    % Fit model to data and create fit function
    [fitresult, gof]    = fit(z_fit, delta_fit, ft, opts)   ; %
    fit_func            = @(x) feval(fitresult,x)           ; %

    % Extract wall deflection coefficients C1, C2, C3
    if isnan(C1_forced) ; C1_fit  = fitresult.C1    ; % [-] C1, cantilever (calculated)
    else                ; C1_fit  = C1_forced       ; % [-] C1, cantilever (forced)
    end
    if isnan(C2_forced) ; C2_fit  = fitresult.C2    ; % [-] C2, parabolic (calculated)
    else                ; C2_fit  = C2_forced       ; % [-] C2, parabolic (forced)
    end
    if isnan(C3_forced) ; C3_fit  = 1 - C1_fit - C2_fit ; % [-] C3, kick-in type (calculated)
    else                ; C3_fit  = C3_forced           ; % [-] C3, kick-in type (forced)
    end
    
    C_fit   = [C1_fit, C2_fit, C3_fit]      ; % [-]

else
    % Calculate the unknown coefficient
    if isnan(C1_forced)
        C1_fit = 1 - C2_forced - C3_forced;
        C2_fit = C2_forced;
        C3_fit = C3_forced;
    elseif isnan(C2_forced)
        C1_fit = C1_forced;
        C2_fit = 1 - C1_forced - C3_forced;
        C3_fit = C3_forced;
    elseif isnan(C3_forced)
        C1_fit = C1_forced;
        C2_fit = C2_forced;
        C3_fit = 1 - C1_forced - C2_forced;
    end

    % Build direct function without fitting
    fit_func = @(x) ...
        C1_fit * beta_fit * ( 2 * Hw_fit                     - 2 * (x - z_fit(1))) + ...
        C2_fit * beta_fit * (-6 / Hw_fit * (x - z_fit(1)).^2 + 6 * (x - z_fit(1))) + ...
        C3_fit * beta_fit * ( 2 * (x - z_fit(1)) );

    % Output results
    C_fit = [C1_fit, C2_fit, C3_fit];
    gof = struct('rsquare', NaN, 'rmse', NaN); % No fit performed
    
end
