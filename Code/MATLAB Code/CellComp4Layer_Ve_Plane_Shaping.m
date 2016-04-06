function [Ve] = CellComp4Layer_Ve_Plane_Shaping(Xi, Yi, Zi, I_M, I_D, I_G,x_max, z_max, t_max, d_x, d_z, d_t, h_F, Ya, rot, Ri)
%%CELLCOMP4LAYER_VE_PLANE_SHAPING
% CELLCOMP4LAYER_VE_PLANE_SHAPING returns the electric field (longitudinal)
% due to subretinal stimulation with point or disk electrodes in a 
% plane in the x-z dimensions (where y is in the direction normal to the
% surface of the retina, toward the electrode). The retina is modelled
% using a 4-layer description, including Insulator, 'Other' Layers, Nerve Fibre
% Layer, and Vitreous (see diagram below for geometry). The integration 
% coefficients are calculated in an accompanying file, 
% CellComp4Layer_Shaping_Coefficients.m.
%
% All units are S.I.
%       _____________________________________________________
%       1.Insulator
%       _____________________________________________________
%                   [|||]        [|||]        [|||]
%                    Electrodes (at layer boundary)
%       
%       2.'Other' cell layer (incl. GCL)
%       _____________________________________________________
%       3.Nerve fibre layer
%       _____________________________________________________
%       4.Vitreous
%       _____________________________________________________
%
% The extracellular potential for the NFL is calculated using a modified version
% of the self-consistent, linear, sub-threshold model presented in:
%
%   B. Tahayori, H. Meffin, E.N. Sergeev, I.M.Y. Mareels, A.N. Burkitt, and
%   D.N. Grayden (2014), "Modelling extracellular electrical stimulation:
%   IV. Effect of the cellular composition of neural tissue on its
%   spatio-temporal filtering properties", J. Neural Eng. 11.
%
% INPUTS:
%
% Xi, Yi and Zi         the x-, y-, and z- coordinates of the point source
%                       electrodes. (m)
% I_M, I_D and I_G      the stimulation amplitude, duration and interphase gap
%                       for each electrode. (A, s, s)
% z_max, x_max, t_max   z-, x- and time-extent of the simulation. (m, s)
% d_z, d_x, d_t         z , x and time step sizes. (m, s)
% h_F                   Nerve fibre layer thickness (m)
% Ya                    Depth of analysis plane below the surface of the
%                       retina (m)
% rot                   coordinate rotation (only set to a nonzero value if the
%                       output will be used to calculate the membrane
%                       potential (rad)
%                       for a fibre which is rotated w.r.t. the fibre bindle.
% Ri                    Radius of disk for each electrode (set to 0 for
%                       point source. (m)
%
% OUTPUTS:
%
% Ve                    The calculated extracellular potential for a full plane.
%
%
% EXAMPLE USAGE:
%
% Xi=0; Yi=30e-6; Zi=0;
% I_M=-1e-6; I_D=100e-6; I_G = 0;
% z_max=1000e-6; x_max=1000e-6; t_max=1000e-6;
% d_z=4e-6; d_x=4e-6; d_t=2e-6;
% h_F = 60e-6; Ya = -30e-6; rot = 0; Ri = 100e-6;
%
% [Ve] = CellComp4Layer_Ve_Plane_Shaping(Xi, Yi, Zi, I_M, I_D, I_G,x_max, ...
% z_max, t_max, d_x, d_z, d_t, h_F, Ya, rot, Ri);
%
%%%%%%%%%%%%%%%%%%%%%%%%% Created by: Tim Esler, 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define data folder for saving and reading from hdf5 files

dataFolder = '..\..\NEURON Code and Data\';

%% Error checking

if (length(Xi)-length(Yi)+1)*(length(Xi)-length(Zi)+1)*...
        (length(Xi)-length(I_M)+1)*(length(Xi)-length(I_D)+1)*...
        (length(Xi)-length(I_G)+1) ~= 1
    error('Inconsistent number of electrodes in electrode specifications (Xi, Yi, Zi, I_M, I_D, I_G)');
end
if length(unique(Yi)) > 1
    error('All electrodes must lie at a single height');
end


%% Define parameters

p = NTESparams('double');         % Call parameter function

% Unpack parameters
b = p.b;                % NTES radius (m)
d = p.d;                % Width of extracellular sheath (m)

C_m = p.C_m;            % Membrane capacitance (F/m^2)
R_m = p.R_m;            % Membrane unit area resistance (ohm.m^2)

rho_i = p.rho_i;        % Intracellular resistivity (ohm.m)
rho_e = p.rho_e;        % Extracellular resistivity (ohm.m)
r_m = p.r_m;            % Membrane unit length resistance (ohm.m)
r_i = p.r_i;            % Intracellular resistance (ohm/m)
r_e = p.r_e;            % Extracellular resistance (ohm/m)

%% Define sampling in time-space and Fourier domains

% Sampling space Fourier domains
kz_max = pi/d_z;
nz = length(d_z:d_z:z_max);
d_kz = kz_max/nz;
kzp = double(-kz_max:d_kz:kz_max);
kzp(fix(length(kzp)/2+1)) = kzp(fix(length(kzp)/2+1))*1e-12;

kx_max = pi/d_x;
nx = length(d_x:d_x:x_max);
d_kx = kx_max/nx;
kxp = double(-kx_max:d_kx:kx_max);
kxp(fix(length(kxp)/2+1)) = kxp(fix(length(kxp)/2+2))*1e-12;

% Sampling time Fourier domain
w_max = pi/d_t;
nt = length(d_t:d_t:t_max);
d_w = w_max/nt;
w = double(-w_max:d_w:w_max);
w(fix(length(w)/2+1)) = 1e-20;

% Sampling space domain
Zp = -z_max:d_z:z_max;

Xp = -x_max:d_x:x_max;

% Sampling time domain
T = -t_max:d_t:t_max;

% Create sample mesh
[kxp_m,kzp_m] = ndgrid(kxp,kzp);

% Apply rotation in Fourier space
kx_m = kxp_m*cos(rot) + kzp_m*sin(rot);
kz_m = -kxp_m*sin(rot) + kzp_m*cos(rot);
clear kxp_m kzp_m

%% Define electrotonic length constants, time constants
% in the Fourier domain

tau_m = R_m*C_m;                % Membrane time constant (s)

% Electrotonic length constants (static and frequency-dependent, for both
% current density (J) and voltage (V) boundary conditions)
L_0J = sqrt(r_m/(r_e+r_i));
L_0V = sqrt(r_m/r_i);

%% Initialise the data file
% This file will hold the calculated extracellular potential, and its
% Fourier domain representation in preliminary steps.

% IMPORTANT: for this simulation, chunk size should be set to [1 length(Zp) 1]
% to speed up IO. This is because the data is read and written twice, once
% taking slices of size length(Xp)*length(Zp) and once taking slices of
% length(Zp)*length(w). The chunk size is the largest common dimension
% between these steps.
if exist([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'], 'file')==2
  delete([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5']);
end
h5create([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_re',[length(Xp) length(Zp) inf], ...
    'FillValue',double(0),'Datatype','double', ...
    'ChunkSize',[1 length(Zp) 1]);
h5create([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_im',[length(Xp) length(Zp) inf], ...
    'FillValue',double(0),'Datatype','double', ...
    'ChunkSize',[1 length(Zp) 1]);

h2 = waitbar(0,'Calculating Extracellular Potential');

%% Iterate through temporal frequencies
% In each iteration of this loop, the extracellular potential is calculated
% for a single value of omega by calculating the 2D IFFT over the spatial
% dimensions.
for j = 1:length(w)
    w_m = w(j);
    L_J_m = L_0J./sqrt(1+1i*w_m*tau_m);
    L_V_m = L_0V./sqrt(1+1i*w_m*tau_m);
    
    %% Define admittivities and conductivities
    
    % NOTE: in this script the conductivities for the vitreous and "other" cell
    % layer are swapped to model subretinal stimulation
    
    % Insulating layer
    sigma_I = 1e-12;                        % Ohm.m
    
    % Lower layers
    sigma_L_xz = 0.1;                     % Ohm.m
    sigma_L_y  = 0.1;    
    
    % Nerve fibre layer (admittivity)
    xi_L_f_m = 1/rho_i * (1+(kz_m.^2).*L_J_m.^2)./(1+(kz_m.^2).*L_V_m.^2);
    xi_T_f = d/b/rho_e;
    
    % Vitreous layer
    sigma_V = 1.78;
    
    %% Define eta coefficients
    
    eta_I = sqrt((kx_m.^2*sigma_I    + kz_m.^2*sigma_I   )./sigma_I  );
    eta_L = sqrt((kx_m.^2*sigma_L_xz + kz_m.^2*sigma_L_xz)./sigma_L_y);
    eta_F = sqrt((kx_m.^2*xi_T_f     + kz_m.^2.*xi_L_f_m )./xi_T_f   );
    eta_V = sqrt((kx_m.^2*sigma_V    + kz_m.^2*sigma_V   )./sigma_V  );
    
    %% Iterate through the point sources
    
    Ve_tmp = zeros(length(Xp),length(Zp));
    
    for i = 1:length(I_M)
        %% Define point source stimulation in the time Fourier domain
        
        % Biphasic pulse - with phase gap, I_G
        % Apply Lanczos sigma factor to reduce Gibbs phenomenon
        sig = sinc(T/d_t*pi*0.5/length(T));
        I_hat = I_D(i)*I_M(i)/sqrt(2*pi).*sinc(I_D(i)*w_m/2/pi)...
            .*(exp(-1i*I_D(i)*w_m/2)-exp(-1i*w_m*(3*I_D(i)/2+I_G(i)))).*sig(j);
        
        %% Calculate extracellular voltage contributed by this source
        % in the Fourier domain
        
        % Define parameter m (which captures the stimulus waveform and
        % electrode geometry)
        if Ri(i) == 0
            m = I_hat./(2*pi*sigma_L_y);
        elseif Ri(i) > 0
            m = I_hat.*besselj(1,sqrt(kx_m.^2 + kz_m.^2).*Ri(i))./...
                (pi*sigma_L_y*Ri(i)*sqrt(kx_m.^2 + kz_m.^2));
        else
            error('Electrode radius should be >= 0');
        end
        
        % Spatial translation of electrode to specified location
        e_shift = exp(-1i*kx_m*Xi(i)-1i*kz_m*Zi(i));
        
        % Extracellular voltage in temporal Fourier domain for each retinal layer
        % GCL and other layers
        if Ya >= 0 && Ya < Yi(i)
            B_1 = CellComp4Layer_Shaping_Coefficients('B1',eta_L,eta_F,eta_V, ...
                eta_I,sigma_L_y,xi_T_f,sigma_V,sigma_I,Yi(i), h_F, m);
            B_2 = CellComp4Layer_Shaping_Coefficients('B2',eta_L,eta_F,eta_V, ...
                eta_I,sigma_L_y,xi_T_f,sigma_V,sigma_I,Yi(i), h_F, m);
            Ve_tmp = Ve_tmp + fftshift(ifft2(ifftshift(2*pi/d_z/d_x* ...
                (B_1.*exp(-eta_L*Ya) + B_2.*exp(eta_L*Ya) + m./(2*eta_L).* ...
                exp(eta_L*(Ya-Yi(i)))).*e_shift)));
        % Nerve fibre layer
        elseif Ya < 0 && Ya >= -h_F
            C_1 = CellComp4Layer_Shaping_Coefficients('C1',eta_L,eta_F,eta_V, ...
                eta_I,sigma_L_y,xi_T_f,sigma_V,sigma_I,Yi(i), h_F, m);
            C_2 = CellComp4Layer_Shaping_Coefficients('C2',eta_L,eta_F,eta_V, ...
                eta_I,sigma_L_y,xi_T_f,sigma_V,sigma_I,Yi(i), h_F, m);
            Ve_tmp = Ve_tmp + fftshift(ifft2(ifftshift(2*pi/d_z/d_x* ...
                (C_1.*exp(-eta_F*Ya) + C_2.*exp(eta_F*Ya)).*e_shift)));
        else
            error('Ya (depth of analysis plane) should be between -h_F and Yi');
        end
    end
    
    % Save current iteration to data file
    h5write([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_re', ...
        double(real(Ve_tmp)),[1 1 j],[length(Xp) length(Zp) 1]);
    h5write([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_im', ...
        double(imag(Ve_tmp)),[1 1 j],[length(Xp) length(Zp) 1]);
    
    waitbar(j/length(w),h2)
end
close(h2)

h2 = waitbar(0,'Calculating 3D FT');

%% Calculate voltage in the time-space domain

% Calculate the IFT along the remaining dimension (time)
for i = 1:length(Xp)
    Ve_tmp_re = h5read([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_re', ...
        [i 1 1],[1 length(Zp) length(w)]);
    Ve_tmp_im = h5read([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_im', ...
        [i 1 1],[1 length(Zp) length(w)]);
    Ve_tmp = complex(Ve_tmp_re,Ve_tmp_im);
    Ve_tmp = real(fftshift(ifft(ifftshift(sqrt(2*pi)/d_t*squeeze(Ve_tmp),2),[],2),2));
    Ve_tmp = permute(Ve_tmp,[3 1 2]);
    h5write([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_re', ...
        real(Ve_tmp),[i 1 1],[1 length(Zp) length(w)]);
    h5write([dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'],'/Ve_im', ...
        imag(Ve_tmp),[i 1 1],[1 length(Zp) length(w)]);
    waitbar(i/length(Xp),h2)
end
close(h2)

Ve = [dataFolder,'Ve_depth',num2str(Ya*1e6),'.h5'];
clearvars -except Ve

end