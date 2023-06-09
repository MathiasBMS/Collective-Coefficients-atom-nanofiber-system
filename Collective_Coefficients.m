
% This code aims to calculate the coefficients in the master equation for a
% system of atoms residing above the surface of a single-mode nanofiber as 
% a function of the interatomic distance a. This is done by splitting the
% coefficients into guided and radiation contribution. The guided
% contribution is calculated analytically, while the radiation contribution
% is calculated by performing the numerical scheme presented in the paper. 

% Written by Mathias B.M. Svendsen, May 2023

clear

% -- Define global variables -- 

global n_1 %Refractive index of the fiber
global n_2 %Refractive index in the medium outside the fiber
global d_alpha %Dipole moment of atom alpha
global d_beta %Dipole moment of atom beta
global r_f %Radius of the fiber
global thetas 
global c
global hbar
global eps_0
global mu_0

% -- Define physical constants --

n_2 = 1; % Index of refraction in vacuum
l_a = 852.347e-9; % Transition wavelength
c = 299792458; % Speed of light
k_a = 2*pi/l_a; % Transition wavenumber
w_a = k_a*c; % Transition frequency
x_a = 100e-9; % Distance from fiber
phi = 0; % Angular coordinate of the atoms
% Planck constant (and reduced planck constant)
h_pl = 6.62607004e-34;
hbar = h_pl ./ (2 .* pi);
% Permittivity of free space
eps_0 = 8.854187817e-12;
% Permeability of free space
mu_0 = 1.25663706e-6;
r_f = 250e-9; % Radius of fiber

conv = 0.0001; % Convergence factor at zero to avoid numerical instabilities
w_c_min=2; % Minimum cut-off for a=1
w_c_max=5; % Maximum cut-off for a=1

a_vec = linspace(1,1,1).*l_a; % Array of atom separations from 1 to 1/10 transition wavelength
w_c_max_vec = w_c_max./(a_vec./l_a); % Max cut-offs for each atom separation a
w_c_min_vec = w_c_min./(a_vec./l_a); % Min cut-offs for each atom separation a
w_c_vec = sort([w_c_min_vec,w_c_max_vec]); % Array of all cut-offs for all values of a
w_vec_norm_L = linspace(conv,w_c_min,floor(30*w_c_min));
w_vec_L = w_vec_norm_L.*c*k_a; % w_vec up to w_c_min
w_vec_norm = w_vec_norm_L;
for j=1:1:(length(w_c_vec)-1)
    w_vec_norm = [w_vec_norm,linspace(w_c_vec(j),w_c_vec(j+1),ceil(30*(w_c_vec(j+1)-w_c_vec(j))))];
end

w_vec = w_vec_norm.*w_a; % w_vec including all cut-offs

% Dipole moment

d_dir = [0,0,1];
dNorm = norm(d_dir);
dMag = 1;

d_alpha = dMag/dNorm.*d_dir;
d_beta = d_alpha;

% Refractive index of fiber by Sellmaier equation
l_am = l_a.* 1e6;
n_1 = sqrt(1 + ((0.6961663 .* (l_am.^2)) ./ (l_am.^2 - 0.0684043.^2)) ...
    + ((0.4079246 .* (l_am.^2)) ./ (l_am.^2 - 0.1162414.^2)) ...
    + ((0.8974794 .* (l_am.^2)) ./ (l_am.^2 - 9.896161.^2)));


gamma = norm(d_alpha)^2*k_a^3/(3*pi*eps_0*hbar); % Spontaneous decay rate

% -- GUIDED PART --

% -- Calculate the guided contribution to the dipole-dipole interaction and
% decay rates --

for t=1:1:length(a_vec)
    r_alpha=[r_f+x_a,phi,0]; % Position of atom alpha
    r_beta = [r_f+x_a,phi,a_vec(t)]; % Position of atom beta
    [V_gd,Gamma_gd] = GuidedVandGamma(r_alpha,r_beta,k_a,r_f);
    V_gd_vec(t) = V_gd;
    Gamma_gd_vec(t) = Gamma_gd;
end

% -- RADIATION PART --

% -- Calculate the longitudinal part of the radiation contribution to the 
% dipole-dipole interaction and decay rates --

for t=1:1:length(a_vec)
    r_alpha=[r_f+x_a,phi,0]; % Position of atom alpha
    r_beta = [r_f+x_a,phi,a_vec(t)]; % Position of atom beta
    r_ab = (r_alpha-r_beta);
    r_ab_hat = r_ab./norm(r_ab);
    dr_prod = dot(d_alpha,r_ab_hat)*dot(conj(d_beta),r_ab_hat);
    dd_prod = dot(d_alpha,d_beta);
    V_rd_long(t) = ((k_a*c)^2)/(hbar*eps_0*c^2)*(3*dr_prod-dd_prod)/(4*pi)*(1./(a_vec(t).^3.*k_a^2)); %Longitudinal contribution to the dipole-dipole interaction
end

% -- Calculate the transverse part of the radiation contribution to the 
% dipole-dipole interaction and decay rates

num_theta = 400; % Number of theta points
sumProf = zeros(length(w_vec),num_theta); % Matrix of the sum of the profile functions for all omega and theta

epsilon = 0.000001; % Convergence factor
thetas = linspace(0+epsilon,pi/2,num_theta); % Array of thetas (symmetric around pi/2)

thresh = 10^(-6); % Threshhold for truncation of mode sum
sumProf_w_a = Sum_Radiation_Profile_Function(thresh, w_a, r_alpha, r_beta); % Calculate the sum of profile functions at w_a
for t=1:1:length(a_vec)
    a = a_vec(t);
    prod_w_a = sumProf_w_a.*exp(1i*w_a/c*cos(thetas)*a); % Product of the sum of the profile functions and the osc function for all thetas and for a specific omega and a
    dImGd = 2*pi/2*trapz(thetas,prod_w_a); % Calculate ImG as the integral over theta, multiply by 2 due to symmetry in theta
    Gamma_rd_vec(t) = 2*w_a^2/(hbar*eps_0*c^2)*dImGd;
end

for j=1:1:length(w_vec)
    w = w_vec(j);
    sumProf(j,:) = Sum_Radiation_Profile_Function(thresh, w, r_alpha, r_beta); % Calculate the sum of profile functions for a specific freq and all theta
end

% -- Calculate the dipole-dipole interaction by using the presented
% numerical solution scheme --

for t=1:1:length(a_vec)
    a=a_vec(t);
    dImGd_vec = zeros(1,length(w_vec));
    for j=1:1:length(w_vec)
        w = w_vec(j);
        prod = sumProf(j,:).*exp(1i*w/c*cos(thetas)*a); % Product of the sum of the profile functions and the osc function for all thetas and for a specific omega and a
        dImGd_vec(j) = 2*pi/2*trapz(thetas,prod); % Calculate ImG as the integral over theta, multiply by 2 due to symmetry in theta
    end
    pos_min = find(w_vec_norm==w_c_min_vec(t)); % Find the position in the omega array to truncate the integral (where w=w_c)
    pos_max = find(w_vec_norm==w_c_max_vec(t)); % Find the position in the omega array to truncate the integral (where w=w_c)

    w_vec_L_a = w_vec(1:pos_min(1)); % Define an array of omega up to the position of the min truncation frequency
    w_vec_R_a = w_vec(pos_min(1):pos_max(1));
    dImGd_L_a = dImGd_vec(1:pos_min(1)); % Define ImG up to the position of truncation
    dImGd_R_a = dImGd_vec(pos_min(1):pos_max(1));
    tau_vec = linspace(0,20.*(a./l_a),1000)./(k_a*c); % Define array of tau-values for integral
    int_tau = zeros(1,length(tau_vec));
    for b=1:1:length(tau_vec)
        tau = tau_vec(b);
        s=0;
        % Do integral over all w_c and avg
        for j=1:1:length(w_vec_R_a)
            w_vec_a = [w_vec_L_a,w_vec_R_a(1:j)]; % Array of frequency up to the truncation frequency
            dImGd_a = [dImGd_L_a,dImGd_R_a(1:j)];
            w_a_vec = linspace(k_a*c,k_a*c,length(w_vec_a)); % Array of transition frequency w_a

            integrand = dImGd_a.*sin(w_vec_a.*tau); % Frequency integrand
            integ = trapz(w_vec_a,integrand); % Integral over frequency
            s=s+integ; % Sum for all cut-off frequencies
        end
        int_tau(b) = s/(length(w_vec_R_a)); % Do average over all cut-offs
    end
    V_rd_avg(t) = 2*(k_a*c)^2/(pi*hbar*eps_0*c^2)*trapz(tau_vec,cos(k_a*c.*tau_vec).*int_tau)./gamma+V_rd_long(t)./gamma; %Do integral over tau
end


% -- FUNCTIONS --

function sum_prof_vec = Sum_Radiation_Profile_Function(thresh, w, r_alpha,r_beta)

%   Sum_Radiation_Profile_Function: Calculating the sum over mode order m 
%   and polarization l of the radiation profile function of the product
%   between the radiation profile functions and the dipole moments. for 
%   num_theta values of theta between 0 and pi/2 for a specific value of 
%   the mode frequency w for atoms at position r_alpha and r_beta. The sum 
%   over the mode order is truncated when the contribution is smaller than 
%   the threshold value thresh.
%
%   Input:
%   thresh - Threshold value determining where to truncate the sum
%   w - Frequency of the radiation mode
%   r_alpha - Coordinates of atom alpha
%   r_alpha - Coordinates of atom beta
%
%   Output:
%   sum_prof_vec - An array of the sum of the product between the radiation
%   profile functions and the dipole moments for all values of theta

global n_1 %Refractive index of the fiber
global n_2 %Refractive index in the medium outside the fiber
global d_alpha %Dipole moment of atom alpha
global d_beta %Dipole moment of atom beta
global r_f %Radius of the fiber
global thetas

c = 299792458; % Speed of light
eps_0 = 8.854187817e-12; % Permeability of free space
% permeability of free space
mu_0 = 1.25663706e-6;

k = w/c;

% Define the Hankel function and derivative
hankelj = @(m,j,x) besselj(m,x)+(-1)^(j-1)*1i*bessely(m,x);
hankeljp = @(m,j,x) 0.5*(hankelj(m-1,j,x)-hankelj(m+1,j,x));

% Bessel function derivative
besseljp = @(m,x) 0.5*(besselj(m-1,x)-besselj(m+1,x));

% Define the coefficients V,M,L as function handles
V1_fh = @(m,beta,kappa,q) m*k*beta/(r_f*kappa^2*q^2)*(n_2^2-n_1^2)*besselj(m,kappa*r_f)*conj(hankelj(m,1,q*r_f));
V2_fh = @(m,beta,kappa,q) m*k*beta/(r_f*kappa^2*q^2)*(n_2^2-n_1^2)*besselj(m,kappa*r_f)*conj(hankelj(m,2,q*r_f));

M1_fh = @(m,beta,kappa,q) 1/kappa*besseljp(m,kappa*r_f)*conj(hankelj(m,1,q*r_f))-1/q*besselj(m,kappa*r_f)*conj(hankeljp(m,1,q*r_f));
M2_fh = @(m,beta,kappa,q) 1/kappa*besseljp(m,kappa*r_f)*conj(hankelj(m,2,q*r_f))-1/q*besselj(m,kappa*r_f)*conj(hankeljp(m,2,q*r_f));

L1_fh = @(m,beta,kappa,q) n_1^2/kappa*besseljp(m,kappa*r_f)*conj(hankelj(m,1,q*r_f))-n_2^2/q*besselj(m,kappa*r_f)*conj(hankeljp(m,1,q*r_f));
L2_fh = @(m,beta,kappa,q) n_1^2/kappa*besseljp(m,kappa*r_f)*conj(hankelj(m,2,q*r_f))-n_2^2/q*besselj(m,kappa*r_f)*conj(hankeljp(m,2,q*r_f));

% Define the components of the radiation profile function outside the fiber
e_r_1_out = @(m,beta,kappa,q,r,C1,D1) 1i/q^2*(beta*q*C1*hankeljp(m,1,q*r)+1i*m*w*mu_0/r*D1*hankelj(m,1,q*r));
e_r_2_out = @(m,beta,kappa,q,r,C2,D2) 1i/q^2*(beta*q*C2*hankeljp(m,2,q*r)+1i*m*w*mu_0/r*D2*hankelj(m,2,q*r));
e_r_out = @(m,beta,kappa,q,r,C1,D1,C2,D2) e_r_1_out(m,beta,kappa,q,r,C1,D1)+e_r_2_out(m,beta,kappa,q,r,C2,D2);

e_phi_1_out = @(m,beta,kappa,q,r,C1,D1) 1i/q^2*(1i*m*beta/r*C1*hankelj(m,1,q*r)-q*w*mu_0*D1*hankeljp(m,1,q*r));
e_phi_2_out = @(m,beta,kappa,q,r,C2,D2) 1i/q^2*(1i*m*beta/r*C2*hankelj(m,2,q*r)-q*w*mu_0*D2*hankeljp(m,2,q*r));
e_phi_out = @(m,beta,kappa,q,r,C1,D1,C2,D2) e_phi_1_out(m,beta,kappa,q,r,C1,D1)+e_phi_2_out(m,beta,kappa,q,r,C2,D2);

e_z_out = @(m,beta,kappa,q,r,C1,C2) C1*hankelj(m,1,q*r)+C2*hankelj(m,2,q*r);

rho_alpha= r_alpha(1); % Radial coordinate of atom alpha
rho_beta = r_beta(1); % Radial coordinate of atom beta
phi = 0; % Angular coordinate equal for all atoms

tic
sum_prof_vec = zeros(1,length(thetas));
l=1; % Polarization

for t=1:1:length(thetas)
    sum_prof = 0; % Initialize to zero
    theta = thetas(t);
    %Parameters associated with the radiation field
    beta_rd = k*cos(theta);
    kappa = sqrt(k^2*n_1^2-beta_rd^2);
    q = sqrt(k^2*n_2^2-beta_rd^2);
    ed_prod=1;
    ed_prod_prev=1;
    m = 0;
    while (ed_prod>thresh || ed_prod_prev>thresh) %Truncate at m_c when ed_prod(m) and ed_prod(m-1) < thresh
        V1 = V1_fh(m,beta_rd,kappa,q);
        V2 = V2_fh(m,beta_rd,kappa,q);

        M1 = M1_fh(m,beta_rd,kappa,q);
        M2 = M2_fh(m,beta_rd,kappa,q);

        L1 = L1_fh(m,beta_rd,kappa,q);
        L2 = L2_fh(m,beta_rd,kappa,q);
        eta = eps_0*c*sqrt((n_2^2*V1*conj(V1)+L1*conj(L1))/(V1*conj(V1)+n_2^2*M1*conj(M1)));

        C1A = (-1)^1*1i*pi*q^2*r_f/(4*n_2^2)*(L1-mu_0*c*l*eta*V1);

        C2A = (-1)^2*1i*pi*q^2*r_f/(4*n_2^2)*(L2-mu_0*c*l*eta*V2);

        D1A = (-1)^1*pi*q^2*r_f/4*(eps_0*c*V1-l*eta*M1);
        D2A = (-1)^2*pi*q^2*r_f/4*(eps_0*c*V2-l*eta*M2);
        A = (16*pi^2*k^2/q^3*(n_2^2*C1A*conj(C1A)+mu_0/eps_0*D1A*conj(D1A)))^(-1/2);

        C1 = C1A*A;
        C2 = C2A*A;
        D1 = D1A*A;
        D2 = D2A*A;

        e_x_out_alpha = cos(phi)*e_r_out(m,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2)-sin(phi)*e_phi_out(m,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2);
        e_y_out_alpha = sin(phi)*e_r_out(m,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2)+cos(phi)*e_phi_out(m,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2);
        e_z_out_alpha = e_z_out(m,beta_rd,kappa,q,rho_alpha,C1,C2);
        e_out_alpha = [e_x_out_alpha,e_y_out_alpha,e_z_out_alpha];

        e_x_out_beta = cos(phi)*e_r_out(m,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2)-sin(phi)*e_phi_out(m,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2);
        e_y_out_beta = sin(phi)*e_r_out(m,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2)+cos(phi)*e_phi_out(m,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2);
        e_z_out_beta = e_z_out(m,beta_rd,kappa,q,rho_beta,C1,C2);
        e_out_beta = [e_x_out_beta,e_y_out_beta,e_z_out_beta];

        ed_prod_prev = ed_prod;
        ed_prod = 2*dot(d_alpha,e_out_alpha)*dot(e_out_beta,d_beta); % Multiply by 2 due to symmetry in l

        sum_prof = sum_prof + ed_prod; %Sum over m
        sum_prof_m(m+1) = ed_prod;
        m=m+1;
    end
    for j=(-(m-1)):1:(-1) % Do the sum from m=-1 to m=-m_c
        V1 = V1_fh(j,beta_rd,kappa,q);
        V2 = V2_fh(j,beta_rd,kappa,q);

        M1 = M1_fh(j,beta_rd,kappa,q);
        M2 = M2_fh(j,beta_rd,kappa,q);

        L1 = L1_fh(j,beta_rd,kappa,q);
        L2 = L2_fh(j,beta_rd,kappa,q);
        eta = eps_0*c*sqrt((n_2^2*V1*conj(V1)+L1*conj(L1))/(V1*conj(V1)+n_2^2*M1*conj(M1)));

        C1A = (-1)^1*1i*pi*q^2*r_f/(4*n_2^2)*(L1-mu_0*c*l*eta*V1);

        C2A = (-1)^2*1i*pi*q^2*r_f/(4*n_2^2)*(L2-mu_0*c*l*eta*V2);

        D1A = (-1)^1*pi*q^2*r_f/4*(eps_0*c*V1-l*eta*M1);
        D2A = (-1)^2*pi*q^2*r_f/4*(eps_0*c*V2-l*eta*M2);
        A = (16*pi^2*k^2/q^3*(n_2^2*C1A*conj(C1A)+mu_0/eps_0*D1A*conj(D1A)))^(-1/2);

        C1 = C1A*A;
        C2 = C2A*A;
        D1 = D1A*A;
        D2 = D2A*A;

        e_x_out_alpha = cos(phi)*e_r_out(j,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2)-sin(phi)*e_phi_out(j,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2);
        e_y_out_alpha = sin(phi)*e_r_out(j,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2)+cos(phi)*e_phi_out(j,beta_rd,kappa,q,rho_alpha,C1,D1,C2,D2);
        e_z_out_alpha = e_z_out(j,beta_rd,kappa,q,rho_alpha,C1,C2);
        e_out_alpha = [e_x_out_alpha,e_y_out_alpha,e_z_out_alpha];

        e_x_out_beta = cos(phi)*e_r_out(j,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2)-sin(phi)*e_phi_out(j,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2);
        e_y_out_beta = sin(phi)*e_r_out(j,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2)+cos(phi)*e_phi_out(j,beta_rd,kappa,q,rho_beta,C1,D1,C2,D2);
        e_z_out_beta = e_z_out(j,beta_rd,kappa,q,rho_beta,C1,C2);
        e_out_beta = [e_x_out_beta,e_y_out_beta,e_z_out_beta];

        sum_prof = sum_prof + 2*dot(d_alpha,e_out_alpha)*dot(e_out_beta,d_beta); %Sum over m, multiply by two due to symmetry in l
    end
    sum_prof_vec(t) = sum_prof;
end
toc
end

function [V_gd,Gamma_gd] = GuidedVandGamma(r_alpha,r_beta,k_a,r_f)

%   GuidedVandGamma: Calculating the guided contribution to the
%   dipole-dipole interaction and the decay rate
%
%   Input:
%   r_alpha - Coordinates of atom alpha
%   r_alpha - Coordinates of atom beta
%   k_a - Transition wavenumber
%   r_f - fiber radius (in m)
%
%   Output:
%   V_gd - Guided contribution to the dipole-dipole interaction strength
%   between two atoms located at position r_alpha and r_beta in vecinity to
%   a fiber of radius r_f
%   Gamma_gd - Guided contribution to the decay rate between two atoms 
%   located at position r_alpha and r_beta in vecinity to a fiber of radius 
%   r_f

global n_1 % Refractive index of the fiber
global n_2 % Refractive index in the medium outside the fiber
global d_alpha % Dipole moment of atom alpha
global d_beta % Dipole moment of atom beta
global hbar % Reduced Planck constant
global c % Speed of light
global eps_0 % Permittivity in vacuum

w_a = k_a*c;
l_a = 2*pi/k_a;

% Hankel function of first and second kind
hankelj = @(j,m,x) besselj(m,x)+((-1)^(j-1))*1i*bessely(m,x);
% Derivative of hankel function
hankeljD = @(j,m,x) 0.5*(hankelj(j,m-1,x)-hankelj(j,m+1,x));
% Derivative of bessel function
besseljD = @(m,x) 0.5*(besselj(m-1,x)-besselj(m+1,x));

% Propagation constant
b_g = guidedModeBeta(n_1,n_2,l_a,r_f);
% Propagation constant derivative
bp = guidedModeBetaD(n_1,n_2,l_a,r_f);

% Parameters associated with the guided field
kappa = ((k_a^2*n_1^2)-b_g^2)^0.5;
q = (b_g^2-(k_a^2*n_2^2))^0.5;

% Normalization constant
C = guidedModeProfileNorm(n_1,n_2,l_a,r_f,b_g);
% Bessel function derivatives
Jp = @(x) 0.5*(besselj(0,x)-besselj(2,x));
Kp = @(x) -0.5*(besselk(0,x)+besselk(2,x));

% Constant
s = ((1/(kappa^2*r_f^2))+(1/(q^2*r_f^2)))/ ...
    ((Jp(kappa*r_f)/((kappa*r_f)*besselj(1,kappa*r_f)))+...
    (Kp(q*r_f)/((q*r_f)*besselk(1,q*r_f))));

% Mode profile function
e_r_g = @(r) 1i*C*(((1-s)*besselk(0,q*r))+((1+s)*besselk(2,q*r)));
e_phi_g = @(r)-C*(((1-s)*besselk(0,q*r))-((1+s)*besselk(2,q*r)));
e_z_g = @(r) C*((2*q)/b_g)*besselk(1,q*r);

Gamma_gd = 0;
V_gd = 0;
% For loop performing the sum over f and p
for f = [-1,1]
    for p = [-1,1]
        z_alpha = r_alpha(3);
        phi_alpha = r_alpha(2);
        rho_alpha =  r_alpha(1);
        e_mu_alpha = (([1;0;0] .* cos(phi_alpha) + [0;1;0] .* sin(phi_alpha)) .* e_r_g(rho_alpha)) ...
                + (([-1;0;0] .* sin(phi_alpha) + [0;1;0] .* cos(phi_alpha)) .*  p.* e_phi_g(rho_alpha)) ...
                + f*([0;0;1] .* e_z_g(rho_alpha));
        z_beta = r_beta(3);
        phi_beta = r_beta(2);
        rho_beta =  r_beta(1);
        e_mu_beta = (([1;0;0] .* cos(phi_beta) + [0;1;0] .* sin(phi_beta)) .* e_r_g(rho_beta)) ...
                + (([-1;0;0] .* sin(phi_beta) + [0;1;0] .* cos(phi_beta)) .*  p.* e_phi_g(rho_beta)) ...
                + f*([0;0;1] .* e_z_g(rho_beta));
        L_mualpha = sqrt(w_a*bp/(4*pi*hbar*eps_0))*dot(conj(d_alpha),e_mu_alpha)*exp(1i*p*phi_alpha+1i*f*b_g*z_alpha); 
        L_mubeta = sqrt(w_a*bp/(4*pi*hbar*eps_0))*dot(conj(d_beta),e_mu_beta)*exp(1i*p*phi_beta+1i*f*b_g*z_beta);
        Gamma_gd = Gamma_gd+2*pi*L_mualpha*conj(L_mubeta);
        if f*(z_alpha-z_beta)>0
            V_gd = V_gd+pi*1i*L_mualpha*conj(L_mubeta);
        else
            V_gd = V_gd-pi*1i*L_mualpha*conj(L_mubeta);
        end
    end
end

end

function beta = guidedModeBeta(n_1, n_2, l_a, r_f)
%   guidedModeBeta: Solution of eigenvalue equation for the propagation
%   constant of a step-profile fiber waveguide
%
%   Input:
%   n_1 - refractive index of fiber core
%   n_2 - refractive index of surrounding medium
%   l_a - free space wavelength (in m)
%   r_f - fiber radius (in m)
%
%   Output:
%   beta - propagation constant (in m^-1)

% Free space transition wavenumber
k_a = (2*pi)/l_a;
% Parameter associated with the fiber
U = @(B) r_f*((k_a^2*n_1^2)-B^2)^0.5;
% Parameter associated with the surroundings
W = @(B) r_f*(B^2-(k_a^2*n_2^2))^0.5;
% Fiber parameter
V = k_a*r_f*(n_1^2-n_2^2)^0.5;

% Bessel function derivatives
Jp = @(x) 0.5*(besselj(0,x)-besselj(2,x));
Kp = @(x) -0.5*(besselk(0,x)+besselk(2,x));

% Terms in the eigenvalue equation
f1 = @(B) (Jp(U(B))/(U(B)*besselj(1,U(B))))+(Kp(W(B))/(W(B)*besselk(1,W(B))));
f2 = @(B) (Jp(U(B))/(U(B)*besselj(1,U(B))))+(n_2^2/n_1^2)*(Kp(W(B))/(W(B)*besselk(1,W(B))));
f3 = @(B) (B/(k_a*n_1))^2*(V/(U(B)*W(B)))^4;

f = @(B) (f1(B)*f2(B))-f3(B);

% The solution is bound by the inequality k_a * n_cl < beta < k_a * n_co
% The interval passed to fzero is truncated slightly using the parameter x,
% which must be large enough such that we avoid the 1/x singularities in f,
% but small enough that we do not lose any solutions
x = 0.0001;
guess = [(k_a*n_2)+x*((k_a*n_1)-(k_a*n_2)),(k_a*n_2)+(1-x)*((k_a*n_1)-(k_a*n_2))];

% Find solution using inbuilt function fzero
beta = fzero(f,guess);

clearvars -except beta

end

function bp = guidedModeBetaD(n_1,n_2,l_a,r_f)
%   guidedModeBetaD: Calculates the derivative of the propagation constant 
%   for a step-profile fiber waveguide with respect to atomic frequency
%
%   Input:
%   n_1 - refractive index of fiber core
%   n_2 - refractive index of surrounding medium
%   l_a - free space wavelength (in m)
%   rho - fiber radius (in m)
%
%   Output:
%   bp - propagation constant derivative

% Free space wavelength values
l_aa = [l_a-(10*1e-9),l_a+(10*1e-9)];
% Corresponding frequencies
w_aa = (2*pi*3e8)./l_aa;

% Bessel function derivatives
Jp = @(x) 0.5*(besselj(0,x)-besselj(2,x));
Kp = @(x) -0.5*(besselk(0,x)+besselk(2,x));

Ba = zeros(1,length(l_aa));

for loop = 1:length(l_aa)
    
    l_a = l_aa(loop);
    
    % Free space transition wavenumber
    k_a = (2*pi)/l_a;
    
    % Parameter associated with the fiber
    U = @(B) r_f*((k_a^2*n_1^2)-B^2)^0.5;
    % Parameter associated with the surroundings
    W = @(B) r_f*(B^2-(k_a^2*n_2^2))^0.5;
    % Fiber parameter
    V = k_a*r_f*(n_1^2-n_2^2)^0.5;

    % Terms in the eigenvalue equation
    f1 = @(B) (Jp(U(B))/(U(B)*besselj(1,U(B))))+(Kp(W(B))/(W(B)*besselk(1,W(B))));
    f2 = @(B) (Jp(U(B))/(U(B)*besselj(1,U(B))))+(n_2^2/n_1^2)*(Kp(W(B))/(W(B)*besselk(1,W(B))));
    f3 = @(B) (B/(k_a*n_1))^2*(V/(U(B)*W(B)))^4;
   
    f = @(B) (f1(B)*f2(B))-f3(B);
    
    % The solution is bound by the inequality k_a * n_cl < beta < k_a * n_co
    % The interval passed to fzero is truncated slightly using the parameter x,
    % which must be large enough such that we avoid the 1/x singularities in f,
    % but small enough that we do not lose any solutions
    x = 0.0001;
    guess = [(k_a*n_2)+x*((k_a*n_1)-(k_a*n_2)),(k_a*n_2)+(1-x)*((k_a*n_1)-(k_a*n_2))];
    
    % Find solution using inbuilt function fzero
    beta = fzero(f,guess);
    
    Ba(loop) = beta;
    
end

% Calculate gradient
bp = (Ba(1)-Ba(2))/(w_aa(1)-w_aa(2));

clearvars -except bp

end


