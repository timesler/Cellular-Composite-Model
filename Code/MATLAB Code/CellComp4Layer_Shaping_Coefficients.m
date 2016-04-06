function C = CellComp4Layer_Shaping_Coefficients(Coeff,eta_L,eta_F,eta_V,eta_I,sigma_L_y,xi_T_f,sigma_V,sigma_I,Yii, h_F, m)
%%CELLCOMP4LAYER_SHAPING_COEFFICIENTS
% CELLCOMP4LAYER_SHAPING_COEFFICIENTS defines the integration constants for the 4 layer
% retinal model. These expressions were generated using the MATLAB symbolic maths engine.
% Coefficients B1 and B2 govern dynamics in the GCL, coefficients C1 and C2 govern 
% dynamics in the NFL and D2 governs dynamics in the vitreous.
%
% See "help CellComp4Layer_Ve_Plane_Shaping" for a description of the 4-layer model geometry.

if strcmpi(Coeff,'B1')
    
    C = (eta_L.*m.*(eta_V.*sigma_V.*exp(Yii.*eta_L + 2.*eta_F.*h_F) + ...
        eta_F.*xi_T_f.*exp(Yii.*eta_L + 2.*eta_F.*h_F) - ...
        eta_V.*sigma_V.*exp(Yii.*eta_L) + ...
        eta_F.*xi_T_f.*exp(Yii.*eta_L)).*sigma_L_y.^2 - ...
        eta_F.*m.*xi_T_f.*(eta_V.*sigma_V.*exp(Yii.*eta_L + 2.*eta_F.*h_F) + ...
        eta_F.*xi_T_f.*exp(Yii.*eta_L + 2.*eta_F.*h_F) + ...
        eta_V.*sigma_V.*exp(Yii.*eta_L) - ...
        eta_F.*xi_T_f.*exp(Yii.*eta_L)).*sigma_L_y)./(exp(2.*eta_F.*h_F).*(eta_I.*sigma_I ...
        - eta_L.*sigma_L_y).*(eta_L.*sigma_L_y - eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) + exp(2.*Yii.*eta_L + 2.*eta_F.*h_F).*(eta_I.*sigma_I + ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y + eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) - ...
        2.*eta_V.*eta_L.^2.*sigma_L_y.^2.*sigma_V.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_I.*sigma_I.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 - 1./2) + ...
        2.*eta_F.*eta_L.^2.*sigma_L_y.^2.*xi_T_f.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_L.*sigma_L_y.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 + 1./2) - ...
        eta_I.*eta_V.*eta_L.*sigma_I.*sigma_L_y.*sigma_V.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_L.*sigma_I.*sigma_L_y.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_V.*sigma_I.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) - 1) + ...
        eta_F.*eta_V.*eta_L.*sigma_L_y.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1));
    
elseif strcmp(Coeff,'B2')
    
    C = -(m.*exp(-Yii.*eta_L).*(eta_I.*sigma_I - ...
        eta_L.*sigma_L_y).*(eta_F.^2.*xi_T_f.^2.*(exp(2.*Yii.*eta_L + 2.*eta_F.*h_F) + 1) ...
        - exp(2.*Yii.*eta_L).*(eta_L.*sigma_L_y - eta_F.*xi_T_f).*(eta_V.*sigma_V - ...
        eta_F.*xi_T_f) + exp(2.*eta_F.*h_F).*(eta_L.*sigma_L_y - ...
        eta_F.*xi_T_f).*(eta_V.*sigma_V + eta_F.*xi_T_f) + ...
        eta_V.*eta_L.*sigma_L_y.*sigma_V.*(exp(2.*Yii.*eta_L + 2.*eta_F.*h_F) - 1) + ...
        eta_F.*eta_L.*sigma_L_y.*xi_T_f.*(exp(2.*Yii.*eta_L + 2.*eta_F.*h_F) + 1) + ...
        eta_F.*eta_V.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L + 2.*eta_F.*h_F) - ...
        1)))./(2.*eta_L.*(exp(2.*eta_F.*h_F).*(eta_I.*sigma_I - ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y - eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) + exp(2.*Yii.*eta_L + 2.*eta_F.*h_F).*(eta_I.*sigma_I + ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y + eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) - ...
        2.*eta_V.*eta_L.^2.*sigma_L_y.^2.*sigma_V.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_I.*sigma_I.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 - 1./2) + ...
        2.*eta_F.*eta_L.^2.*sigma_L_y.^2.*xi_T_f.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_L.*sigma_L_y.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 + 1./2) - ...
        eta_I.*eta_V.*eta_L.*sigma_I.*sigma_L_y.*sigma_V.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_L.*sigma_I.*sigma_L_y.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_V.*sigma_I.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) - 1) + ...
        eta_F.*eta_V.*eta_L.*sigma_L_y.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1)));
    
elseif strcmpi(Coeff,'C1')
    
    C = -(2.*eta_L.*m.*sigma_L_y.^2.*exp(Yii.*eta_L).*(eta_V.*sigma_V - ...
        eta_F.*xi_T_f))./(exp(2.*eta_F.*h_F).*(eta_I.*sigma_I - ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y - eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) + exp(2.*Yii.*eta_L + 2.*eta_F.*h_F).*(eta_I.*sigma_I + ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y + eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) - ...
        2.*eta_V.*eta_L.^2.*sigma_L_y.^2.*sigma_V.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_I.*sigma_I.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 - 1./2) + ...
        2.*eta_F.*eta_L.^2.*sigma_L_y.^2.*xi_T_f.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_L.*sigma_L_y.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 + 1./2) - ...
        eta_I.*eta_V.*eta_L.*sigma_I.*sigma_L_y.*sigma_V.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_L.*sigma_I.*sigma_L_y.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_V.*sigma_I.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) - 1) + ...
        eta_F.*eta_V.*eta_L.*sigma_L_y.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1));
    
elseif strcmpi(Coeff,'C2')
    
    C = (2.*eta_L.*m.*sigma_L_y.^2.*exp(Yii.*eta_L + ...
        2.*eta_F.*h_F).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f))./(exp(2.*eta_F.*h_F).*(eta_I.*sigma_I - ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y - eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) + exp(2.*Yii.*eta_L + 2.*eta_F.*h_F).*(eta_I.*sigma_I + ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y + eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) - ...
        2.*eta_V.*eta_L.^2.*sigma_L_y.^2.*sigma_V.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_I.*sigma_I.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 - 1./2) + ...
        2.*eta_F.*eta_L.^2.*sigma_L_y.^2.*xi_T_f.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_L.*sigma_L_y.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 + 1./2) - ...
        eta_I.*eta_V.*eta_L.*sigma_I.*sigma_L_y.*sigma_V.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_L.*sigma_I.*sigma_L_y.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_V.*sigma_I.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) - 1) + ...
        eta_F.*eta_V.*eta_L.*sigma_L_y.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1));
    
elseif strcmpi(Coeff,'D2')
    
    C = (4.*eta_F.*eta_L.*m.*sigma_L_y.^2.*xi_T_f.*exp(Yii.*eta_L).*exp(eta_F.*h_F).* ...
        exp(eta_V.*h_F))./(exp(2.*eta_F.*h_F).*(eta_I.*sigma_I ...
        - eta_L.*sigma_L_y).*(eta_L.*sigma_L_y - eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) + exp(2.*Yii.*eta_L + 2.*eta_F.*h_F).*(eta_I.*sigma_I + ...
        eta_L.*sigma_L_y).*(eta_L.*sigma_L_y + eta_F.*xi_T_f).*(eta_V.*sigma_V + ...
        eta_F.*xi_T_f) - ...
        2.*eta_V.*eta_L.^2.*sigma_L_y.^2.*sigma_V.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_I.*sigma_I.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 - 1./2) + ...
        2.*eta_F.*eta_L.^2.*sigma_L_y.^2.*xi_T_f.*(exp(2.*Yii.*eta_L)./2 - 1./2) - ...
        2.*eta_F.^2.*eta_L.*sigma_L_y.*xi_T_f.^2.*(exp(2.*Yii.*eta_L)./2 + 1./2) - ...
        eta_I.*eta_V.*eta_L.*sigma_I.*sigma_L_y.*sigma_V.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_L.*sigma_I.*sigma_L_y.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1) + ...
        eta_F.*eta_I.*eta_V.*sigma_I.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) - 1) + ...
        eta_F.*eta_V.*eta_L.*sigma_L_y.*sigma_V.*xi_T_f.*(exp(2.*Yii.*eta_L) + 1));
    
end