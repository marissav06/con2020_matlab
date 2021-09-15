function brtp = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, varargin) % varagin = use_these_params in this context
% ======
% CON2020_MODEL_RTP (Spherical)
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
%  Use default current sheet model parameter structure:  B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads)
%  Use your own current sheet model parameter structure: B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, use_these_params)
%  Obtain the default model parameters:             params = con2020_model_rtp('default_values')
%  Then you can edit the structure and use, e.g. params.r1__outer_rj = 51.400001,
%    then B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, params)
%
% Required inputs (System III Spherical, right handed):
%  eq_type - equation type: 'integral', 'analytic' or 'hybrid',
%   or set to 'default_values' to return a structure of all default values.
%  r_rj       - radial distance, in Rj.                    Value(s) must be 0 <      r_rj   <  200.
%  colat_rads - colatitude, in radians.                    Value(s) must be 0 <= colat_rads <=  pi.
%  elong_rads - East longitude, right handed, in radians.  Value(s) must be 0 <= elong_rads <= 2pi.
% r_rj, colat_rads and elong_rads can be scalars or 1D arrays (nx1), but only one eq_type.
%
% Unless an option structure is provided it will default to parameters from Connerney et al., 2020.
% Optional input of a structure: use_these_params
% with the structure fields:
%  use_these_params.mu_i_div2__current_density_nT           - mu0i0/2 term (current sheet current density), in nT
%  use_these_params.i_rho__radial_current_density_nT        - radial current term from Connerney et al., 2020 (set this to zero to turn radial currents off as in Connerney et al. 1981)
%  use_these_params.r0__inner_rj                            - inner edge of current disk in Rj
%  use_these_params.r1__outer_rj                            - outer edge of current disk in Rj
%  use_these_params.d__cs_half_thickness_rj                 - current sheet half thickness in Rj
%  use_these_params.xt__cs_tilt_degs                        - current sheet tilt in degrees
%  use_these_params.xp__cs_rhs_azimuthal_angle_of_tilt_degs - current sheet longitude (right handed) in degrees
%  use_these_params.error_check                             - 1 to check that inputs are valid (Default),
%                                                             or set to 0 to skip input checks (faster).
%
% To retrieve the default value parameters, run
% params = con2020_model_rtp('default_values')
% then you may edit the values in the structure, then use as B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, params)
%
% Outputs:
%  B - Spherical Magnetic field vector from current sheet model, [Br, Btheta, Bphi], units of nT.
%
% This code takes a hybrid approach to calculating the current sheet field, using the integral equations in some regions
% and the analytic equations in others.
% Following Connerney et al. 1981, figure A1, and Edwards et al. (2001), figure 2, the choice of integral vs. analytic
% equations is most important near rho = r0 and z = 0.
% By default, this code uses the analytic equations everywhere except |Z| < D*1.5 and |Rho-R0| < 2
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

if strcmpi(eq_type,'default_values') % case insensitive. % Not yet checked if eq_type is a character string, Matlab doesn't care
    brtp = con2020_model_xyz(eq_type);
    return
end

% Covert to Doubles, and rename input variables so as not to over-write them (only an IDL issue)
r_in  = double(    r_rj  ); % units Rj
theta = double(colat_rads); % SYSIII Colat (radians)
phi   = double(elong_rads); % SYSIII ELong or Azimuth (radians)

% Very basis error check that theta and phi are probably in radians
while 1
    if numel(varargin) == 1
        if varargin{1}.error_check == 1
            % skip error check that values are probably radians.
            break
        end
    end
    % Not checking value of r_rj here.  Cartesian Code will check x, y and z are all within +/-200.
    if numel(theta) == 1
        if  (theta < 0) || (    theta  >   pi); error('Input colat_rads must be in radians and from 0 to  pi only!'); end
        if  (phi   < 0) || (    phi    > 2*pi); error('Input elong_rads must be in radians and from 0 to 2pi only!'); end
        break
    end
    % assume vector, not scalar
    if (min(theta) < 0) || (max(theta) >   pi); error('Input colat_rads must be in radians and from 0 to  pi only!'); end
    if (min(phi  ) < 0) || (max(phi  ) > 2*pi); error('Input elong_rads must be in radians and from 0 to 2pi only!'); end
    break % escape WHILE loop!
end

% Convert to cartesian coordinates and rotate into magnetic longitude
% (x,y,z) are the shifted (phi) coordinates
sin_theta = sin(theta);
cos_theta = cos(theta);
sin_phi   = sin(phi);
cos_phi   = cos(phi);

x = r_in.*sin_theta.*cos_phi;
y = r_in.*sin_theta.*sin_phi;
z = r_in.*cos_theta;

% Do calculation in cartesian
switch numel(varargin)
    case 0
        bxyz = con2020_model_xyz(eq_type, x, y, z);
    case 1
        bxyz = con2020_model_xyz(eq_type, x, y, z, varargin{1});
    otherwise
        error('ERROR: code accepts either 4 or 5 input arguments, not %d.',numel(varargin)+4)
end

% Convert to spherical coordinates
brtp = [...
    bxyz(:,1).*cos_phi.*sin_theta + bxyz(:,2).*sin_phi.*sin_theta + bxyz(:,3).*cos_theta, ... %br
    bxyz(:,1).*cos_phi.*cos_theta + bxyz(:,2).*sin_phi.*cos_theta - bxyz(:,3).*sin_theta, ... %bt
    -bxyz(:,1).*sin_phi           + bxyz(:,2).*cos_phi                                    ... %bp
    ]; % size n x 3
end
