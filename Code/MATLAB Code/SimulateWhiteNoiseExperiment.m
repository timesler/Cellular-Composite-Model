%%SIMULATEWHITENOISEEXPERIMENT
% SIMULATEWHITENOISEEXPERIMENT prepares calculated extracellular potential input 
% for driving a cell in a NEURON simulation, by completing the following steps:
% 	1. Define simulation geometry and parameters:
% 		i.  Tissue geometry: sets soma location, reads in a neuron reconstruction,
% 		    define thickness of each of layers.
% 		ii. Electrode geometry: defines electrode location and stimulus amplitude,
% 			defines the simulation's spatial and temporal extent and sampling
% 			parameters.
%			Note that all coordinates are rotated after being defined so that the cell
% 			axon aligns with the orientation of the nerve fibre layer.
%	2. Calculate the time course of extracellular potential in a specify volume via
%      repeated calls to CellComp4Layer_Ve_Plane_Shaping (each call calculates
% 	   potential in a plane.
% 	3. Plot the simulation geometry, including the electrode array and the cell
% 	   reconstruction, with each domain (e.g. soma, axon, dendrite, SOCB) color-coded.
% 	4. From the 3D grid of extracellular potential (Ve), interpolate its value at the
% 	   centre of each cell compartment. This is done separately for each electrode, 
% 	   which is appropriate since the effect of each electrode can be combined using
% 	   superposition later (the CCM is linear).
% 	5. Calculate the "activating function" for each compartment. This is the equivalent
% 	   injected current induced by the extracellular stimulation.
% 	6. Save AF data to text file which can then be read by NEURON and used to drive an
% 	   "IClamp" mechanism in each cell.
%
% See "help CellComp4Layer_Ve_Plane_Shaping" for a description of the 4-layer model.
%
% The cellular composite model was first introduced in:
% 
%   B. Tahayori, H. Meffin, E.N. Sergeev, I.M.Y. Mareels, A.N. Burkitt, and
%   D.N. Grayden (2014), "Modelling extracellular electrical stimulation:
%   IV. Effect of the cellular composition of neural tissue on its
%   spatio-temporal filtering properties", J. Neural Eng. 11.
%
% Created by: Tim Esler, 2015
%%

addpath('..\Trees Package 1.15')
dataFolder = '..\..\NEURON Code and Data\';

%% Parameters

cellName = '2016Jan26_c2';
mkdir([dataFolder,cellName]) % folder for generated data and code

% Define location of soma
Ys = 10e-6; % this is equivalent to 10 um into the GCL (see "help CellComp4Layer_Ve_Plane_Shaping")
Xs = 0e-6;
Zs = 1000e-6*sin(pi/6)/sin(4*pi/6);
axonAction = 3; % 1 = artificially extend axon, 2 = leave axon as is, 3 = remove axon except initial segment

regenVe = false; % Do we need to recalculate extracellular potential or use the saved data?

% Define simulation size and step sizes (a single simulation produces Ve for a single y-plane)
t_max = 1500e-6;
x_max = 5000e-6;
z_max = 5000e-6;
d_z = 10e-6;
d_x = d_z;
d_y = 4e-6;
d_t = 5e-6;

% Path of cell reconstruction
SWCPath = '..\..\NEURON Code and Data\Cells\Reconstructions\2014Sept08_c2.swc';

% Nerve fibre layer thickness and rotation (0 rad = align NFL with z-axis)
h_F = 40e-6;
rot = 0;

% Load tree and process
start_trees;
load_tree(SWCPath);
tree = trees{1};
[tree] = PrepCell_Shaping(tree,Xs,Ys,Zs,axonAction);
Xc = tree.X;
Yc = tree.Y;
Zc = tree.Z;

% Define sampling vectors (Y is set to span the volume occupied by the reconstructed cell)
% Here we define the size of the volume that is simulated
X = -x_max:d_x:x_max;
Y = -12e-6:d_y:36e-6;
Z = -z_max:d_z:z_max;
T = -t_max:d_t:t_max;

% Define location of base electrode
% A wide field simulation is run with a single electrode.
% Since all electrodes are at the same y-height, the output of this can just be scaled
% and translated to get the potential from the remaining electrodes. The influence 
% of all electrodes can then be summed linearly at each point in space.
Xi = 0e-6;
Yi = 360e-6;
Zi = 0e-6;
Ri = 200e-6; % electrode disc radius

% Define electrode array geometry
X_MEA = [-2000; -1000; 0; 1000; ...
    -1500; -500; 500; 1500; ...
    -2000; -1000; 0; 1000; ...
    -1500; -500; 500; 1500; ...
    -2000; -1000; 0; 1000]*1e-6;
Y_MEA = ones(20,1)*Yi;
Z_MEA = [-2000*cos(pi/6); -2000*cos(pi/6); -2000*cos(pi/6); -2000*cos(pi/6); ...
    -1000*cos(pi/6); -1000*cos(pi/6); -1000*cos(pi/6); -1000*cos(pi/6); ...
    0; 0; 0; 0; ...
    1000*cos(pi/6); 1000*cos(pi/6); 1000*cos(pi/6); 1000*cos(pi/6); ...
    2000*cos(pi/6); 2000*cos(pi/6); 2000*cos(pi/6); 2000*cos(pi/6)]*1e-6;

% Define point source amplitudes and pulse durations
I_M = 100e-6;
I_D = 500e-6;
I_G = 50e-6;

SimSize = length(Y);

%% Rotate all coordinates such that the axon of the RGC cell aligns with the z-axis

% Clockwise rotation of 30 degrees
rotation = 15*pi/36;
rot_mat = [cos(rotation) 0 -sin(rotation);...
    0 1 0;...
    sin(rotation) 0 cos(rotation)];

% Soma coords
XYZs = rot_mat*[Xs; Ys; Zs];
Xs = XYZs(1);
Ys = XYZs(2);
Zs = XYZs(3);

% Cell coords
XYZc = rot_mat*[Xc'; Yc'; Zc'];
Xc = XYZc(1,:)';
Yc = XYZc(2,:)';
Zc = XYZc(3,:)';
tree.X = Xc;
tree.Y = Yc;
tree.Z = Zc;

% Electrode coords
XYZi = rot_mat*[Xi; Yi; Zi];
Xi = XYZi(1);
Yi = XYZi(2);
Zi = XYZi(3);

% Cell coords
XYZ_MEA = rot_mat*[X_MEA'; Y_MEA'; Z_MEA'];
X_MEA = XYZ_MEA(1,:)';
Y_MEA = XYZ_MEA(2,:)';
Z_MEA = XYZ_MEA(3,:)';

%% Calculate the extracellular voltage for a single electrode

if regenVe

    h = waitbar(0,'Calculating Rotated Membrane Potential');

	% iterate through each y-layer occupied by cell and calculate Ve in that plane
    for j = 1:SimSize

        disp(['Current plane: ', num2str(Y(j)*1e6)])
		
		% Calculate timecourse of Ve in current plane
        [Ve_L_plane] = CellComp4Layer_Ve_Plane_Shaping(...
            Xi, Yi, Zi, I_M, I_D, I_G, ...              % Electrode parameters
            x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
            h_F, Y(j), rot, Ri ...                      % Others
            );

        waitbar(j/SimSize,h)

    end

    close(h)

end

%% Save simulation parameters
save([dataFolder,cellName,'\SingleElectrodeWideField_Params_' cellName '.mat'],'t_max','x_max','z_max', ...
    'd_x','d_y','d_z','d_t','SWCPath','Ys','Xs','Zs','tree','X','Y', ...
    'Z','T','Xi','Yi','Zi','Ri','I_M','I_D','I_G','X_MEA','Y_MEA','Z_MEA', ...
    'dataFolder','cellName','Xc','Yc','Zc')

%% Plot simulation geometry

% scrsz = get(groot,'MonitorPositions');
% P4 = figure('outerposition',[scrsz(2,1) 31 scrsz(2,3)/2 scrsz(2,4)-31],'Color',[1 1 1]);

for i = 1:length(X_MEA)
    plot3(Z_MEA(i)*1e6,X_MEA(i)*1e6,Y_MEA(i)*1e6,'d','MarkerSize',6,'MarkerFaceColor','b')
    hold on
end
for i = 1:length(Xi)
    plot3(Zi(i)*1e6,Xi(i)*1e6,Yi(i)*1e6,'d','MarkerSize',6,'MarkerFaceColor','k')
    hold on
end
plot3(Zs*1e6,Xs*1e6,Ys*1e6,'d','MarkerSize',6,'MarkerFaceColor','r')
% axis equal
grid minor
xlabel('Z (\mum)')
ylabel('X (\mum)')
zlabel('Y (\mum)')
xlim([-2750 2750])
ylim([-2750 2750])
for i = 1:length(Xc)
    switch tree.R(i)
        case 1
            col = 'k';
        case 2
            col = 'b';
        case 3
            col = 'g';
        case 4
            col = 'r';
        case 5
            col = 'c';
        case 6
            col = 'm';
    end
    plot3(Zc(i)*1e6,Xc(i)*1e6,Yc(i)*1e6,'d','MarkerSize',2,'MarkerFaceColor',col,...
        'MarkerEdgeColor',col)
    hold on
end
hold off
view(0,90)

%% Interpolate Ve at each cell compartment in 3 dimensions

Ve = TrilinearInterpolation(Xc,Yc,Zc,X_MEA,Y_MEA,Z_MEA,dataFolder,cellName);

save([dataFolder,cellName,'\SingleElectrodeWideField_Ve_' cellName '.mat'],'Ve','Xc','Yc','Zc', ...
    'X_MEA','Y_MEA','Z_MEA')

%% Calculate activating function from the morphology and extracellular potential
% The activating function is the equivalent "injected" current in a particular 
% compartment given the extracellular potential (Ve) along the neurite.

AF = ActivatingFunction(tree,Ve);

save([dataFolder,cellName,'\SingleElectrodeWideField_AF_' cellName '.mat'],'AF','tree')

% Sort the activating function data by region (soma, axon, dendrite, initialseg,
% SOCB, narrowr). This is the order used when constructing the cell topology in
% generateNEURONTopology.m

AF = [AF(:,tree.R == 1,:) AF(:,tree.R == 2,:) AF(:,tree.R == 3,:) ...
    AF(:,tree.R == 4,:) AF(:,tree.R == 5,:) AF(:,tree.R == 6,:)];

% Downsample (lowpass filter) the AF matrix to a suitable NEURON
% sampling rate
AF_d = zeros(length(X_MEA),length(Xc),ceil(length(T)/5));
for i = 1:length(X_MEA)
    for j = 1:length(Xc)
        AF_d(i,j,:) = decimate(squeeze(AF(i,j,:)),5);
    end
end
T_d = T(1:5:end);

% Set final AF value to 0 for all segment (this value will be used as the
% injection for the rest of the NEURON simulation after 3 ms)
AF_d(:,:,end) = 0;

for i = 1:length(X_MEA)
    AF_E = squeeze(AF_d(i,:,:));
    dlmwrite([dataFolder,cellName,'\',cellName,'_AF_E',num2str(i),'.txt'],AF_E(:)*1e9,' ')
end

dlmwrite([dataFolder,cellName,'\',cellName,'_tvec.txt'],T_d*1e3-min(T_d*1e3),' ')

%% Create .nrn file for model topology and dynamics

p = NTESparams('double');

% Format parameters in the neuron tree structure using NEURON units
tree.X = tree.X*1e6;
tree.Y = tree.Y*1e6;
tree.Z = tree.Z*1e6;
tree.D = tree.D*1e6;
tree.Ra = p.rho_i*1e2;
tree.tstop = 2*t_max*1e3;
tree.dt = d_t*5*1e3;

generateNEURONTopology(tree,[dataFolder,cellName,'\',cellName,'_RGC.hoc']);

tree.X = tree.X*1e-6;
tree.Y = tree.Y*1e-6;
tree.Z = tree.Z*1e-6;
tree.D = tree.D*1e-6;

%% Prepare electrode amplitude file for reading by NEURON

wn_data = dlmread([dataFolder,cellName,'\AllStim_M2.txt'],' ',3,0);
wn_amps = wn_data(:,1:20);
dlmwrite([dataFolder,cellName,'\',cellName,'_ElectrodeAmplitudes.txt'],wn_amps(:),' ');