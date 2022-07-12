function Bxyz = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, varargin)
% ======
% CON2020_MODEL_XYZ (Cartesian)
% ======
% Code to calculate the perturbation magnetic field produced by the Connerney et al. 1981 (CAN) current sheet, which is
% represented by a finite disk of current.
%  This disk has variable parameters including (among others) the current density, and current sheet inner edge, outer
%   edge and thickness.
%  The disk is centered on the magnetic equator (shifted in longitude and tilted as specified by model parameters xp__cs_rhs_azimuthal_angle_of_tilt_degs and xt__cs_tilt_degs)
%   parameters of an internal field model like VIP4 or JRM09)
%  This 2020 version includes a radial current per Connerney et al. (2020),
%   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JA028138
%  For more details about the model and the development of this code please see the PDF at 
%   https://github.com/marissav06/con2020_idl/blob/main/con2020_final_code_documentation_sept13_2021.pdf
%
% Use in one of the following ways:
%  Use default current sheet model parameter structure:  B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj) 
%  Use your own current sheet model parameter structure: B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, use_these_params)
%  Obtain the default model parameters:             params = con2020_model_xyz('default_values')
%  Then you can edit the structure and use, e.g. params.r1__outer_rj = 51.400001,
%    then B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, params)
%
% Required inputs (System III Cartesian, right handed):
%  eq_type - equation type: 'integral', 'analytic' or 'hybrid',
%   or set to 'default_values' to return a structure of all default values.
%  x_rj      - SYSIII x position, in Rj, Values must be -200 < x_rj < 200.
%  y_rj      - SYSIII y position, in Rj, Values must be -200 < x_rj < 200.
%  z_rj      - SYSIII z position, in Rj, Values must be -200 < x_rj < 200.
% x_rj, y_rj and z_rj can be scalars or 1D arrays (nx1), but only one eq_type.
%
% Unless an option structure is provided it will default to parameters from Connerney et al., 2020.
% Optional input of a structure: use_these_params
% with the structure fields:
%  use_these_params.mu_i_div2__current_density_nT           - mu0i0/2 term (current sheet current density), in nT
%  use_these_params.i_rho__radial_current_intensity_MA      - radial current term from Connerney et al., 2020 (set this to zero to turn radial currents off as in Connerney et al. 1981)
%  use_these_params.r0__inner_rj                            - inner edge of current disk in Rj
%  use_these_params.r1__outer_rj                            - outer edge of current disk in Rj
%  use_these_params.d__cs_half_thickness_rj                 - current sheet half thickness in Rj
%  use_these_params.xt__cs_tilt_degs                        - current sheet tilt in degrees
%  use_these_params.xp__cs_rhs_azimuthal_angle_of_tilt_degs - current sheet longitude (right handed) in degrees
%  use_these_params.error_check                             - 1 to check that inputs are valid (Default),
%                                                             or set to 0 to skip input checks (faster).
%
% To retrieve the default value parameters, run 
% params = con2020_model_xyz('default_values')
% then you may edit the values in the structure, then use as B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, params)
%
% Outputs:
%  B - Cartesian Magnetic field vector from current sheet model, [Bx, By, Bz], units of nT.
%
% This code can take a hybrid approach to calculating the current sheet field, using the integral equations in some regions
% and the analytic equations in others.
% Following Connerney et al. 1981, figure A1, and Edwards et al. (2001), figure 2, the choice of integral vs. analytic
% equations is most important near rho = r0 and z = 0.
% By default, this hybrid method uses the analytic equations everywhere except |Z| < D*1.5 and |Rho-R0| < 2
%    Analytic equations:
%        For the analytic equations, we use the equations provided by Edwards et al. 2001:
%         https://doi.org/10.1016/S0032-0633(00)00164-1
%        Other analytic approximations to the CAN sheet equations are provided in Connerney et al., 1981
%         https://doi.org/10.1029/JA086iA10p08370
%    Integral equations:
%        For the integral equations we use the Bessel functions from Connerney et al. 1981, eqs. 14, 15, 17, 18
%        We do not integrate lambda from zero to infinity, but vary the integration limit depending on the value of the
%        Bessel functions.
%
% Updates:
% by Marissa Vogt (mvogt@bu.edu), March 2021,
% RJ Wilson (Rob.Wilson@lasp.colorado.edu) did some speedups and re-formatting of lines, also March 2021
%
% Converted to MATLAB by Marty Brennan (martin.brennan@jpl.nasa.gov), June 2021
% RJ Wilson did some reformatting, June 2021, and added
% int_tabulated_rjw2_sub as a subfunction, rather than separate file int_tabulated_rjw2.m
% which was then replaced by some in-line code instead of calling the subfunction
% RJ Wilson split initial Matlab and IDL code in to Cartesian and a Spherical wrapper code and updated this help text,
% in August 2021, to make con2020_model_xyz and con2020_model_rtp.
% RJW Wilson renamed i_rho__radial_current_density_nT to i_rho__radial_current_intensity_MA in June 2022.

switch numel(varargin) % faster than if exist('use_these_params','var')
    case 1
        use_these_params = varargin{1};    	
        if ~isstruct(use_these_params)
            error('Must be a structure of terms to use in code')
        end
        mu_i_div2__current_density_nT           = double(use_these_params.mu_i_div2__current_density_nT          );
        r0__inner_rj                            = double(use_these_params.r0__inner_rj                           );
        r1__outer_rj                            = double(use_these_params.r1__outer_rj                           );
        d__cs_half_thickness_rj                 = double(use_these_params.d__cs_half_thickness_rj                );
        xt__cs_tilt_degs                        = double(use_these_params.xt__cs_tilt_degs                       );
        xp__cs_rhs_azimuthal_angle_of_tilt_degs = double(use_these_params.xp__cs_rhs_azimuthal_angle_of_tilt_degs);
        i_rho__radial_current_intensity_MA      = double(use_these_params.i_rho__radial_current_intensity_MA     );
        error_check = double(use_these_params.error_check);
        
        if strcmpi(eq_type,'default_values') % case insensitive % Not yet checked if eq_type is a character string
            error('ERROR: Can not request default_values and also provide a set of values to use. To get default values, just pass the first argument.')
        end
    case 0
        % Make sure all of these numbers are doubles!
        Default_values.mu_i_div2__current_density_nT           = 139.6; % current density (nT)
        Default_values.r0__inner_rj                            =   7.8; % inner radius (Rj)
        Default_values.r1__outer_rj                            =  51.4; % outer radius (Rj)
        Default_values.d__cs_half_thickness_rj                 =   3.6; % half-height  (Rj)
        Default_values.xt__cs_tilt_degs                        =   9.3; % dipole tilt (Deg.)
        Default_values.xp__cs_rhs_azimuthal_angle_of_tilt_degs = 155.8; % dipole longitude (right handed) (Deg.), Table 1 xp = 204.2 but that value is in left handed SIII
        Default_values.i_rho__radial_current_intensity_MA      =  16.7; % radial current term from Connerney et al., 2020, in mega-Amps.
        % NOTE: The default value (16.7 MA) is the average value from Connerney et al 2020. This value was shown to vary from one
        % pass to the next, where Table 2 (units of MA) provides radial current intensity values for 23 of the first 24 perijoves
        % (units mistakenly given as nT in the article text).
        Default_values.error_check = 1; % input error check: 1 = yes, 0 = no
        
        if strcmpi(eq_type,'default_values') % case insensitive. % Not yet checked if eq_type is a character string, Matlab doesn't care
            Bxyz = Default_values;
            fprintf('Returning structure of Default terms used in code.\n');
            return
        end
        
        mu_i_div2__current_density_nT           = Default_values.mu_i_div2__current_density_nT          ;
        r0__inner_rj                            = Default_values.r0__inner_rj                           ;
        r1__outer_rj                            = Default_values.r1__outer_rj                           ;
        d__cs_half_thickness_rj                 = Default_values.d__cs_half_thickness_rj                ;
        xt__cs_tilt_degs                        = Default_values.xt__cs_tilt_degs                       ;
        xp__cs_rhs_azimuthal_angle_of_tilt_degs = Default_values.xp__cs_rhs_azimuthal_angle_of_tilt_degs;
        i_rho__radial_current_intensity_MA      = Default_values.i_rho__radial_current_intensity_MA     ;
        error_check = Default_values.error_check ;
        
    otherwise
        error('ERROR: code accepts either 4 or 5 input arguments, not %d.',numel(varargin)+4)
end

Deg2Rad = pi/180;

N_input = length(x_rj);
scalar_input = (N_input == 1); % scalar or not

eq_type = lower(eq_type);

if error_check(1)
    % Check optional input values to be numeric, and doubles!
    if  (~isnumeric(error_check)) || (numel(error_check)~=1), error('ERROR: Field  error_check must be a scalar double number');end
    if  (~isnumeric(mu_i_div2__current_density_nT          )) || (numel(mu_i_div2__current_density_nT          )~=1), error('ERROR: Field       mu_i_div2__current_density_nT     must be a scalar double number');end
    if  (~isnumeric(r0__inner_rj                           )) || (numel(r0__inner_rj                           )~=1), error('ERROR: Field               r0__inner_rj              must be a scalar double number');end
    if  (~isnumeric(r1__outer_rj                           )) || (numel(r1__outer_rj                           )~=1), error('ERROR: Field               r1__outer_rj              must be a scalar double number');end
    if  (~isnumeric(d__cs_half_thickness_rj                )) || (numel(d__cs_half_thickness_rj                )~=1), error('ERROR: Field          d__cs_half_thickness_rj        must be a scalar double number');end
    if  (~isnumeric(xt__cs_tilt_degs                       )) || (numel(xt__cs_tilt_degs                       )~=1), error('ERROR: Field             xt__cs_tilt_degs            must be a scalar double number');end
    if  (~isnumeric(xp__cs_rhs_azimuthal_angle_of_tilt_degs)) || (numel(xp__cs_rhs_azimuthal_angle_of_tilt_degs)~=1), error('ERROR: Field xp__cs_rhs_azimuthal_angle_of_tilt_degs must be a scalar double number');end
    if  (~isnumeric(i_rho__radial_current_intensity_MA     )) || (numel(i_rho__radial_current_intensity_MA     )~=1), error('ERROR: Field    i_rho__radial_current_intensity_MA   must be a scalar double number');end
    
    %if  (error_check~=0) && (error_check~=1), error('ERROR: Field  error_check must be 0 or 1');end
    % if here error_check must be scalar and not 0
    if error_check~=1, error('ERROR: Field  error_check must be 0 or 1');end
    error_check = true; 
    if mu_i_div2__current_density_nT <= 0           , error('mu_i_div2__current_density_nT must be > 0'           ),end
    if                  r0__inner_rj <= 0           , error(                 'r0__inner_rj must be > 0'           ),end
    if                  r1__outer_rj <= r0__inner_rj, error(                 'r1__outer_rj must be > r0__inner_rj'),end
    if       d__cs_half_thickness_rj <= 0           , error(      'd__cs_half_thickness_rj must be > 0'           ),end


    % Check if equation input is type char and matches a selection
    if (isa(eq_type,'char') == 0) || (size(eq_type,1)~=1)
        error('ERROR: First argument equation_type must be a character string and can only be ''hybrid'' (default), ''analytic'' or ''integral''');
    end
    if sum(strcmp(eq_type,{'integral','analytic','hybrid'}))==0
        error('ERROR: First argument equation_type must be a string and can only be ''hybrid'' (default), ''analytic'' or ''integral''');
    end
    
    % Check inputs x_rj, y_rj and z_rj are all numbers, and same size (also scalar or 1D only)
    if (~isnumeric(x_rj)) || (size(x_rj,2) ~= 1), error('ERROR: Second argument x_rj must be a scalar number or 1D column array of numbers'); end
    if (~isnumeric(y_rj)) || (size(y_rj,2) ~= 1), error('ERROR: Third  argument y_rj must be a scalar number or 1D column array of numbers'); end
    if (~isnumeric(z_rj)) || (size(z_rj,2) ~= 1), error('ERROR: Fourth argument z_rj must be a scalar number or 1D column array of numbers'); end
    if (N_input ~= length(y_rj)), error('ERROR: Second argument x_rj must be the same size as 3rd argument y_rj'); end
    if (N_input ~= length(z_rj)), error('ERROR: Second argument x_rj must be the same size as 4th argument z_rj'); end
    % do this check to be sure that user hasn't got position in km.
    if scalar_input
        if (x_rj <= -200) || (x_rj >= 200) || (y_rj <= -200) || (y_rj >= 200) || (z_rj <= -200) || (z_rj >= 200)
            error('ERROR: Positions x_rj, y_rj and z_rj must all be in units of Rj, and each >-200 and <200 only, and not outside that range (did you use km instead?)');
        end
    else
        if (min(x_rj) <= -200) || (max(x_rj) >= 200), error('ERROR: Second argument, Position x_rj, must be in units of Rj and >-200 and <200 only, and not outside that range (did you use km instead?)'); end
        if (min(y_rj) <= -200) || (max(y_rj) >= 200), error('ERROR: Third  argument, Position y_rj, must be in units of Rj and >-200 and <200 only, and not outside that range (did you use km instead?)'); end
        if (min(z_rj) <= -200) || (max(z_rj) >= 200), error('ERROR: Fourth argument, Position z_rj, must be in units of Rj and >-200 and <200 only, and not outside that range (did you use km instead?)'); end
    end
end

xp__cs_rhs_azimuthal_angle_of_tilt_degs = xp__cs_rhs_azimuthal_angle_of_tilt_degs - 180.0; % shift needed because of the way we define the rotation angle 
dipole_shift = xp__cs_rhs_azimuthal_angle_of_tilt_degs * Deg2Rad; % xp__cs_rhs_azimuthal_angle_of_tilt_degs is longitude of the dipole. dipole_shift used here and at end of code
theta_cs     = xt__cs_tilt_degs                        * Deg2Rad; % dipole tilt is xt__cs_tilt_degs
cos_dipole_shift = cos(dipole_shift);
sin_dipole_shift = sin(dipole_shift);
cos_theta_cs     = cos(theta_cs);
sin_theta_cs     = sin(theta_cs);

% Covert to Doubles, and rename input variables so as not to over-write them (only an IDL issue)
xSYSIII = double(x_rj);
ySYSIII = double(y_rj);
zSYSIII = double(z_rj);

% Two rotations are needed to align the IAU_JUPITER frame with the MAG frame
% first by (if using defaults) +155.2 degrees about Z, second by -9.3 degrees about Y.
% so first rotate around Z
x0 =  xSYSIII.*cos_dipole_shift + ySYSIII.*sin_dipole_shift;
y0 = -xSYSIII.*sin_dipole_shift + ySYSIII.*cos_dipole_shift;
z0 =  zSYSIII; % no change as we rotate about Z
% Second rotate around Y
x1 =  x0.*cos_theta_cs + z0.*sin_theta_cs;
y1 =  y0; % no chance as we rotate about Y %RJW - NOT NEEDED??? - But used in ATAN equivalent code near end of function (could just use y)
z1 = -x0.*sin_theta_cs + z0.*cos_theta_cs;
% Calculate rho
rho1_sq = x1.*x1 + y1.*y1;
rho1 = sqrt(rho1_sq); %cylindrical radial distance


abs_z1 = abs(z1);

% Decide whether to use integral equations or analytic equations
switch eq_type(1) % instead of whole char string, just take first letter
    case 'a' % 'analytic'
        if scalar_input
            do_integral = false;
        else
            do_integral = zeros(N_input,1,'logical');
        end
    case 'h' %'hybrid'
        % scalar form is IF ((abs_z1 LE d__cs_half_thickness_rj*1.5d) and (abs(rho1-r0__inner_rj) LE 2.0d)) THEN do_integral = 1 ELSE do_integral = 0;
        do_integral = (abs_z1 <= d__cs_half_thickness_rj*1.5) & (abs(rho1-r0__inner_rj) <= 2);
    case 'i' %'integral'
        if scalar_input
            do_integral = true;
        else
            do_integral = ones( N_input,1,'logical');
        end
    otherwise
        error('ERROR: case statement has unrecognized string - was your equation_type lower case?') % if no_error_check, this will still catch bad inputs for equation_type
end

if scalar_input
    brho1 = 0;
    bz1   = 0;
    if do_integral == 0
        n_ind_analytic = 1;
        n_ind_integral = 0;
        ind_analytic   = 1;
        ind_integral   = [];
    else
        n_ind_analytic = 0;
        n_ind_integral = 1;
        ind_analytic   = [];
        ind_integral   = 1;
    end
else
    ind_analytic   = find(do_integral == 0);
    n_ind_analytic = length(ind_analytic);
    ind_integral   = find(do_integral == 1); % ==1 is really ~=0, but as logical, this is the same.
    n_ind_integral = length(ind_integral);
    % preallocate arrays with NaNs (For safety)
    brho1 = NaN(N_input,1); % set all to NaN
    bz1   = brho1;            % copy of Brho of array of NaNs
end



if (n_ind_integral ~= 0)
    %     ;% Integral equations - Brho and Bz eqs. vary depending on region with respect to the disk
    dlambda_brho    = 1E-4  ;% default step size for Brho function
    dlambda_bz      = 5E-5  ;% default step size for Bz function
    
    check1 = abs(abs_z1(ind_integral) -  d__cs_half_thickness_rj);
    check2 =     abs_z1(ind_integral) <= d__cs_half_thickness_rj*1.1; % 1 or 0 (true or false)
    for zcase = 1:6
        switch zcase
            case 1
                check3 = (check2 == 1) & (check1 >= 0.7);
                lambda_max_brho =   4 ;% default integration limit for Brho function;
                lambda_max_bz   = 100 ;
            case 2
                check3 = (check2 == 0) & (check1 >= 0.7);
                lambda_max_brho =   4 ;% default integration limit for Brho function;
                lambda_max_bz   =  20 ;% Small Z or default in else
            case 3
                check3 = (check2 == 1) & (check1 >= 0.1) & (check1 < 0.7);
                lambda_max_brho =  40 ;% Z > D             ;
                lambda_max_bz   = 100 ;
            case 4
                check3 = (check2 == 0) & (check1 >= 0.1) & (check1 < 0.7);
                lambda_max_brho =  40 ;% Z > D             ;
                lambda_max_bz   =  20 ;% Small Z or default in else
            case 5
                check3 = (check2 == 1) & (check1 < 0.1);
                lambda_max_brho = 100 ;% Z very close to D
                lambda_max_bz   = 100 ;
            case 6
                check3 = (check2 == 0) & (check1 < 0.1);
                lambda_max_brho = 100 ;% Z very close to D;
                lambda_max_bz   =  20 ;% Small Z or default in else
            otherwise
                error('Should never get to this else case statement')
        end
        
        if scalar_input
            if ~check3, continue; end
            ind_case   = 1;
            n_ind_case = 1;
            %         n_ind_case_minus_1 = 0; % DO NOT NEED THIS FOR MATLAB INDEXING
        else
            % ind_case = find(check3 == 1);
            ind_case = find(check3);
            n_ind_case = length(ind_case);
            if ~n_ind_case
                continue
            end
            %         n_ind_case_minus_1 = n_ind_case ; % DO NOT NEED THIS FOR MATLAB INDEXING
        end
        
        lambda_int_brho = dlambda_brho:dlambda_brho:(lambda_max_brho-dlambda_brho);
        lambda_int_bz   = dlambda_bz  :dlambda_bz  :(lambda_max_bz  -dlambda_bz  );
        
        beselj_rho_r0_0   = besselj(0,lambda_int_brho*r0__inner_rj); % Only 6 sets of values
        beselj_z_r0_0     = besselj(0,lambda_int_bz  *r0__inner_rj); % Only 6 sets of values
        
        for zi = 1:n_ind_case
            ind_for_integral = ind_integral(ind_case(zi));% sub-indices of sub-indices!
            
            besselj_rho_rho1_1 = besselj(1,lambda_int_brho*rho1(ind_for_integral));
            besselj_z_rho1_0   = besselj(0,lambda_int_bz  *rho1(ind_for_integral));
            
            if (abs_z1(ind_for_integral) > d__cs_half_thickness_rj) % Connerney et al. 1981 eqs. 14 and 15
                brho_int_funct = besselj_rho_rho1_1.*beselj_rho_r0_0.*sinh(d__cs_half_thickness_rj*lambda_int_brho).*exp(-abs_z1(ind_for_integral).*lambda_int_brho)./lambda_int_brho;
                bz_int_funct   = besselj_z_rho1_0  .*beselj_z_r0_0  .*sinh(d__cs_half_thickness_rj*lambda_int_bz  ).*exp(-abs_z1(ind_for_integral).*lambda_int_bz  )./lambda_int_bz  ;
                
                %brho1(ind_for_integral) = mu_i_div2__current_density_nT*2*int_tabulated_rjw2_sub(lambda_int_brho,brho_int_funct);  %#ok<AGROW>
                % Can use Trapezoidal rule approx
                brho1(ind_for_integral) = mu_i_div2__current_density_nT*2*dlambda_brho * (sum(brho_int_funct) - (brho_int_funct(1) + brho_int_funct(end))*0.5 ); %#ok<AGROW>

                if z1(ind_for_integral) < 0, brho1(ind_for_integral) = -brho1(ind_for_integral); end %#ok<AGROW>
                
            else % Connerney et al. 1981 eqs. 17 and 18
                brho_int_funct = besselj_rho_rho1_1.*beselj_rho_r0_0.*(    sinh(z1(ind_for_integral).*lambda_int_brho).*exp(-d__cs_half_thickness_rj*lambda_int_brho))./lambda_int_brho;
                bz_int_funct   = besselj_z_rho1_0  .*beselj_z_r0_0  .*(1 - cosh(z1(ind_for_integral).*lambda_int_bz  ).*exp(-d__cs_half_thickness_rj*lambda_int_bz  ))./lambda_int_bz  ;
                
                %brho1(ind_for_integral) = mu_i_div2__current_density_nT*2*int_tabulated_rjw2_sub(lambda_int_brho,brho_int_funct);
                % Can use Trapezoidal rule approx
                brho1(ind_for_integral) = mu_i_div2__current_density_nT*2* dlambda_brho * (sum(brho_int_funct) - (brho_int_funct(1) + brho_int_funct(end))*0.5 ); %#ok<AGROW>
            end
            
            %bz1(ind_for_integral)   = mu_i_div2__current_density_nT*2*int_tabulated_rjw2_sub(lambda_int_bz  ,bz_int_funct  );
            % Can use Trapezoidal rule approx
            bz1(ind_for_integral)   = mu_i_div2__current_density_nT*2* dlambda_bz * (sum(bz_int_funct) - (bz_int_funct(1) + bz_int_funct(end))*0.5 ); %#ok<AGROW>
            
        end
        if scalar_input, break; end % If only scalar, stop doing for loops once we found the one we want
    end
end

% Work out the finite sheet here - re-use some bits for the final analysis
[brho_finite, bz_finite] =                     con2020_model_xyz_analytic( rho1, z1, rho1_sq, d__cs_half_thickness_rj, r1__outer_rj, mu_i_div2__current_density_nT, scalar_input);

% Do now near Analytic bit
if (n_ind_analytic ~= 0)
    [brho1(ind_analytic), bz1(ind_analytic)] = con2020_model_xyz_analytic( rho1(ind_analytic), z1(ind_analytic), rho1_sq(ind_analytic), d__cs_half_thickness_rj, r0__inner_rj, mu_i_div2__current_density_nT, scalar_input);
end

% Check that brho1 and bz1 do not contain NaNs or Infs
if error_check
    if sum(isfinite(brho1)) ~= N_input, error('ERROR: some brho1 values are NaN or Infs'); end
    if sum(isfinite( bz1 )) ~= N_input, error('ERROR: some  bz1  values are NaN or Infs'); end
end
% RJW - DO SCALAR VERSION HERE FOR SPEED


% New to CAN2020 (not included in CAN1981): radial current produces an azimuthal field, so Bphi is nonzero
bphi1 = 2.7975*i_rho__radial_current_intensity_MA./rho1;
if scalar_input
    if abs_z1 < d__cs_half_thickness_rj, bphi1 =  bphi1 * abs_z1 / d__cs_half_thickness_rj; end
    
    if z1     > 0                      , bphi1 = -bphi1                                   ; end
else
    ind = find(abs_z1 < d__cs_half_thickness_rj);
    if ~isempty(ind), bphi1(ind) =  bphi1(ind) .* abs_z1(ind) / d__cs_half_thickness_rj; end
    
    ind = find(    z1 > 0);
    if ~isempty(ind), bphi1(ind) = -bphi1(ind)                                         ; end
end


% Account for finite nature of current sheet by subtracting the field values
% calculated brho_finite and bz_finite earlier
brho1       = brho1 - brho_finite;
bz1         = bz1   - bz_finite  ;


% brho1, bphi1, and bz1 here are the ultimately calculated brho and bz values from the CAN model
% the remaining calculations just rotate the field back into SIII

%Calculate 'magnetic longitude' and convert the field into cartesian coordinates
cos_phi1 = x1./rho1 ;
sin_phi1 = y1./rho1 ;

bx1 = brho1.*cos_phi1 - bphi1.*sin_phi1;
by1 = brho1.*sin_phi1 + bphi1.*cos_phi1;

% Now convert back to SYSIII

% Rotate back by dipole tilt amount, into coordinate system that is aligned with Jupiter's spin axis
bx = bx1.*cos_theta_cs - bz1.*sin_theta_cs;
%by = by1; % just using by1 below
bz = bx1.*sin_theta_cs + bz1.*cos_theta_cs;

% Finally, shift back to SIII longitude
bx2 = bx .*cos_dipole_shift - by1.*sin_dipole_shift; % used sin(-a) = asin(a) & cos(-a) = cos(a)
by2 = by1.*cos_dipole_shift + bx .*sin_dipole_shift;
% bz2 = bz so just using bz below 

Bxyz = [bx2, by2, bz]; % size n x 3

end
%%
% function trap = int_tabulated_rjw2_sub(x,f)
% % Numerical Integration based on Trapezoidal rule only https://www.math24.net/trapezoidal-rule
% % Assumes equally spaced x-data
% nseg = length(x) - 1;
% h    = (x(end) - x(1)) / nseg;
% trap = h * (sum(f) - (f(1) + f(end))*0.5 );
% end
%%
function [brho1, bz1] = con2020_model_xyz_analytic(rho1, z1, rho1_sq, d__cs_half_thickness_rj,r, mu_i_div2__current_density_nT, scalar_input)

    % Analytic equations
    % Connerney et al. 1981's equations for the field produced by a semi-infinite disk of thickness D, inner edge R0, outer edge R1 -
    %  see their equations A1 through A9
    % the analytic equations for Brho and Bz vary depending on the region with respect to the current disk
    
    % Doing these 3 equations on the whole array to save getting confused by indices,  Will swap to just required indices later
    z1md = z1-d__cs_half_thickness_rj; % Matt's zmd
    z1pd = z1+d__cs_half_thickness_rj; % Matt's zpd
    
    r_sq = r*r; % Matt's a2
    
    if scalar_input
        if rho1 < r
            ind_LT   = 1;
            n_ind_LT = 1;
            %         ;% ind_GE = [];
            n_ind_GE = 0;
        else
            %         ;% ind_LT = [];
            n_ind_LT = 0;
            ind_GE   = 1;
            n_ind_GE = 1;
        end
        brho1 = NaN;
        bz1   = NaN;
    else
        ind_LT = find(rho1 <  r); n_ind_LT = length(ind_LT);  %Matt's small approx
        ind_GE = find(rho1 >= r); n_ind_GE = length(ind_GE);  %Matt's large approx
        brho1 = NaN(numel(rho1),1);
        bz1   = brho1;
    end
    
    if n_ind_LT %(n_ind_LT ~= 0) % Matt's small approx
        zmd2 = z1md(ind_LT).*z1md(ind_LT);
        zpd2 = z1pd(ind_LT).*z1pd(ind_LT);
        f1_sq = zmd2 + r_sq;
        f2_sq = zpd2 + r_sq;
        f1 = sqrt(f1_sq);
        f2 = sqrt(f2_sq);
        f1_cubed = f1_sq.*f1;
        f2_cubed = f2_sq.*f2;
        
        % calculate the terms which make equations 9a and 9b
     	rhoov2   = rho1(ind_LT)/2;
     	rho2ov4  = rhoov2  .* rhoov2;
     	rho3ov16 = rho2ov4 .* rhoov2/2;

        % these bits are used to form 9a
     	f3 = (r_sq - 2*zmd2)./(f1_sq.*f1_sq.*f1);
        f4 = (r_sq - 2*zpd2)./(f2_sq.*f2_sq.*f2);

        terma0 = rhoov2  .*(1./f1 - 1./f2);
        terma1 = rho3ov16.*(f3    - f4   );
     	brho1(ind_LT) = mu_i_div2__current_density_nT*(terma0 + terma1);
 	
        % now equation 9b
        termb0 = log((z1pd(ind_LT) + f2)./(z1md(ind_LT) + f1));
    	termb1 = rho2ov4.*(z1pd(ind_LT)./f2_cubed - z1md(ind_LT)./f1_cubed);
        bz1(  ind_LT) = mu_i_div2__current_density_nT*(termb0 + termb1);
    end
    if n_ind_GE % i.e. (n_ind_GE ~= 0)
        zmd2 = z1md(ind_GE).*z1md(ind_GE);
        zpd2 = z1pd(ind_GE).*z1pd(ind_GE);
        f1_sq = zmd2 + rho1_sq(ind_GE);
        f2_sq = zpd2 + rho1_sq(ind_GE);
        f1 = sqrt(f1_sq);
        f2 = sqrt(f2_sq);
        f1_cubed = f1_sq.*f1;
        f2_cubed = f2_sq.*f2;

      
    	%equation 13a
        terma0 = (1./rho1(ind_GE)).*(f1 - f2);
        terma1 = (rho1(ind_GE).*r_sq/4).*(1./f2_cubed - 1./f1_cubed);
        % from Python z1_clip = z1.clip(max=D,min=-D)
        if scalar_input
            if     z1 >  d__cs_half_thickness_rj
                z1_clip =  d__cs_half_thickness_rj;
            elseif z1 < -d__cs_half_thickness_rj
                z1_clip = -d__cs_half_thickness_rj;
            else
                z1_clip = z1;
            end
        else
            z1_clip = z1(ind_GE);
            z1_clip(z1_clip >  d__cs_half_thickness_rj) =  d__cs_half_thickness_rj;
            z1_clip(z1_clip < -d__cs_half_thickness_rj) = -d__cs_half_thickness_rj;
        end
        terma2 = (2./rho1(ind_GE)).*z1_clip;
        brho1(ind_GE) = mu_i_div2__current_density_nT*(terma0 + terma1 + terma2);
	
        %equation 13b - same as before
   	    bz1(ind_GE) = mu_i_div2__current_density_nT*(log((z1pd(ind_GE)+f2)./(z1md(ind_GE)+f1)) + (r_sq/4).*((z1pd(ind_GE)./f2_cubed) - (z1md(ind_GE)./f1_cubed)));
    end
end 
