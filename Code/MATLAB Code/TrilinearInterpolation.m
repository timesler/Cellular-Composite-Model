function [Ve] = TrilinearInterpolation(Xc,Yc,Zc,X_MEA,Y_MEA,Z_MEA,dataFolder,cellName)

load([dataFolder,cellName,'\SingleElectrodeWideField_Params_',cellName,'.mat'])
Ve = zeros(length(X_MEA),length(Xc),length(T));
    
h = waitbar(0,'3D Interpolation');

for j = 1:length(X_MEA)
    %% Load simulation parameters
    load([dataFolder,cellName,'\SingleElectrodeWideField_Params_',cellName,'.mat'])
    
    Xc_e = Xc - X_MEA(j);
    Yc_e = Yc;
    Zc_e = Zc - Z_MEA(j);
    
    %% Perform 3D linear interpolation
    for i = 1:length(Xc_e)
        disp(i)
        % Get 8 nearest grid coordinates
        [x_idx,x_d] = knnsearch(X',Xc_e(i),'K',2);
        [y_idx,y_d] = knnsearch(Y',Yc_e(i),'K',2);
        [z_idx,z_d] = knnsearch(Z',Zc_e(i),'K',2);
        x_d = x_d/d_x;
        y_d = y_d/d_y;
        z_d = z_d/d_z;
        
        % Put lower index first
        if x_idx(1) > x_idx(2)
            x_idx = [x_idx(2) x_idx(1)];
            x_d = [x_d(2) x_d(1)];
        end
        if y_idx(1) > y_idx(2)
            y_idx = [y_idx(2) y_idx(1)];
            y_d = [y_d(2) y_d(1)];
        end
        if z_idx(1) > z_idx(2)
            z_idx = [z_idx(2) z_idx(1)];
            z_d = [z_d(2) z_d(1)];
        end
        
        currentFile1 = [dataFolder,'Ve_depth',num2str(Y(y_idx(1))*1e6),'.h5'];
        currentFile2 = [dataFolder,'Ve_depth',num2str(Y(y_idx(2))*1e6),'.h5'];
        
        Ve_Vol_Y1 = squeeze(h5read(currentFile1,'/Ve_re', ...
            [x_idx(1) z_idx(1) 1],[2 2 length(T)]));
        Ve_Vol_Y2 = squeeze(h5read(currentFile2,'/Ve_re', ...
            [x_idx(1) z_idx(1) 1],[2 2 length(T)]));
        
        % Get 8 corresponding extracellular potential values
        Ve_000 = squeeze(Ve_Vol_Y1(1,1,:));
        Ve_100 = squeeze(Ve_Vol_Y1(2,1,:));
        Ve_001 = squeeze(Ve_Vol_Y1(1,2,:));
        Ve_101 = squeeze(Ve_Vol_Y1(2,2,:));
        Ve_010 = squeeze(Ve_Vol_Y2(1,1,:));
        Ve_110 = squeeze(Ve_Vol_Y2(2,1,:));
        Ve_011 = squeeze(Ve_Vol_Y2(1,2,:));
        Ve_111 = squeeze(Ve_Vol_Y2(2,2,:));
        
        % Interpolate in the x-direction
        Ve_00 = Ve_000*x_d(2) + Ve_100*x_d(1);
        Ve_01 = Ve_001*x_d(2) + Ve_101*x_d(1);
        Ve_10 = Ve_010*x_d(2) + Ve_110*x_d(1);
        Ve_11 = Ve_011*x_d(2) + Ve_111*x_d(1);
        
        % Interpolate in the y-direction
        Ve_0 = Ve_00*y_d(2) + Ve_10*y_d(1);
        Ve_1 = Ve_01*y_d(2) + Ve_11*y_d(1);
        
        %Interpolate in the z-direction
        Ve(j,i,:) = Ve_0*z_d(2) + Ve_1*z_d(1);
        
        waitbar(i/length(Xc),h,['3D Interpolation: ' num2str(j)])
        
        
    end
    clearvars -EXCEPT Ve h j Xc Yc Zc X_MEA Y_MEA Z_MEA dataFolder cellName
end

close(h)