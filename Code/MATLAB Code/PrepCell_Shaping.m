function [tree] = PrepCell_Shaping(tree,Xs,Ys,Zs,axonAction)
% Perfoms a series of processing steps on the cell reconstruction, such as:
% 	1. rotations and translations
%	2. adding/extending/removing the axon
% 	3. defines compartment regions: soma, axon, dendrite, initial segment, SOCB, narrow region
% 	4. adds some electrical properties based on segment length and diameter - these are used
% 	   for the activating function calculation later.

% Resample tree to have equal length segments
comp_len = 10;
tree = resample_tree(tree,comp_len);

% Read in reconstruction data (swapping Y and Z)
Xc = tree.X;
Zc = tree.Y;
Yc = tree.Z;
Rc = tree.R;

idpar = idpar_tree(tree);
seglen = len_tree(tree);

% Translate soma to origin
Xc = Xc - Xc(1);
Yc = Yc - Yc(1);
Zc = Zc - Zc(1);

% Get axon
Xa = Xc(Rc==2);
Ya = Yc(Rc==2);
Za = Zc(Rc==2);
Ra = Rc(Rc==2);

% Get array of cumulative distance and find max distance
% In this loop, I also define the initial segment, SOCB and narrow region
cumDist = zeros(1,length(Xa));

tree.rnames = {'soma','axon','dend','initseg','SOCB','narrowr','unknown'};
Ra(1) = 4;

for j = 2:length(Xa)
    cumDist(j) = cumDist(j-1) + sqrt((Xa(j)-Xa(j-1))^2 + (Ya(j)-Ya(j-1))^2 + (Za(j)-Za(j-1))^2);
    if cumDist(j) <= 30
        Ra(j) = 4;
    elseif cumDist(j) > 30 && cumDist(j) <= 70
        Ra(j) = 5;
    elseif cumDist(j) > 70 && cumDist(j) <= 130
        Ra(j) = 6;
    end
end

Rc(Rc==2) = Ra;

maxDist = max(cumDist);

% Find vector of last 50 um of axon
[~, ind_last] = min(abs(maxDist-cumDist-50));
last_proj = [Xa(end)-Xa(ind_last); 0; Za(end)-Za(ind_last)];

% Align the last axon segment with the direction of the 9th electrode
angle_proj = atan2(last_proj(1),last_proj(3));

Xc_rot = Xc*cos(angle_proj+pi/2) - Zc*sin(angle_proj+pi/2);
Zc_rot = Xc*sin(angle_proj+pi/2) + Zc*cos(angle_proj+pi/2);

Xc = Xc_rot;
Zc = Zc_rot;

% Translate cell so that soma is at desired location
Xc = Xc*1e-6 + Xs;
Yc = Yc*1e-6 + Ys;
Zc = Zc*1e-6 + Zs;

if axonAction == 1 % extend
    % Extend axon until optic disk
    % Get distal axon
    Xa = Xc(Rc==2);
    Ya = Yc(Rc==2);
    Za = Zc(Rc==2);
    XE9 = -2000e-6;
    ZE9 = 0;
    NumSegs = round((sqrt((Xa(end)-XE9)^2+(Za(end)-ZE9)^2)+500e-6)/(comp_len*1e-6));
    Xext = Xa(end)-(1:NumSegs)*(comp_len*1e-6)*cos(pi/3);
    Zext = Za(end)-(1:NumSegs)*(comp_len*1e-6)*sin(pi/3);
    Yext = Ya(end)+(1:NumSegs)*comp_len*1e-8;

    Xc = [Xc; Xext'];
    Yc = [Yc; Yext'];
    Yc = 2*Ys - Yc;
    Zc = [Zc; Zext'];

    % Find last axon segment in data
    lastA = find(Rc == 2);
    lastA = lastA(end);

    idpar = [idpar; lastA; [1:(NumSegs-1)]'+length(idpar)];
    seglen = [seglen*1e-6; sqrt((Xext(2) - Xext(1))^2 + ...
        (Yext(2) - Yext(1))^2 + ...
        (Zext(2) - Zext(1))^2)*ones(NumSegs,1)]; 
    Rc = [Rc; 2*ones(NumSegs,1)];
    D = [tree.D; mean(tree.D)*ones(NumSegs,1)]*1e-6;

    tree.idpar = idpar;
    tree.dA = blkdiag(tree.dA,circshift(eye(NumSegs),1));
    tree.dA(length(Xc)-NumSegs+1,:) = 0;
    tree.dA(length(Xc)-NumSegs+1,lastA) = 1;
    
    tree.X = Xc;
    tree.Y = Yc;
    tree.Z = Zc;
    tree.D = D;
    tree.R = Rc;
    tree.seglen = seglen;
elseif axonAction == 2 % do nothing
    Yc = 2*Ys - Yc;
    seglen = seglen*1e-6;
    D = tree.D;
    tree.idpar = idpar;
    
    tree.X = Xc;
    tree.Y = Yc;
    tree.Z = Zc;
    tree.D = D;
    tree.R = Rc;
    tree.seglen = seglen;
else % remove axon except initial segment
    Yc = 2*Ys - Yc;
    seglen = seglen*1e-6;
    D = tree.D*1e-6;
    tree.idpar = idpar;
    tree.R = Rc;
    tree = delete_tree(tree,find(tree.R == 2),'-r');
    tree.X = Xc(Rc ~= 2);
    tree.Y = Yc(Rc ~= 2);
    tree.Z = Zc(Rc ~= 2);
    tree.D = D(Rc ~= 2);
    tree.R = Rc(Rc ~= 2);
    tree.seglen = seglen(Rc ~= 2);
end

tree.seglen(1) = tree.D(1);
tree.name = [tree.name '_processed'];

% Add electrical properties to cell
p = NTESparams('double');

tree.Ri = p.rho_i*tree.seglen./(pi*tree.D.^2/4);
tree.Cm = p.C_m*tree.seglen*pi.*tree.D;

