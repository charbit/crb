function [CRB, Jacobav_k] = evalCRBwithLOC(xsensors_m, sigma2noise, C, polparam, fK_Hz)

[M,D]    = size(xsensors_m);
slowness_spm    = zeros(D,1);



switch D
    case 2
        cosa  = cos(polparam.a_deg*pi/180);
        sina  = sin(polparam.a_deg*pi/180);
        v_mps = polparam.v_mps;
        slowness_spm(1) = cosa/v_mps;
        slowness_spm(2) = sina/v_mps;
        delay_s         = xsensors_m * slowness_spm;
    case 3
        cosa  = cos(polparam.a_deg*pi/180);
        sina  = sin(polparam.a_deg*pi/180);
        cose  = cos(polparam.e_deg*pi/180);
        sine  = sin(polparam.e_deg*pi/180);
        c_mps = polparam.c_mps;
        v_mps = c_mps/cose;
        slowness_spm(1) = cosa*cose/c_mps;
        slowness_spm(2) = sina*cose/c_mps;
        slowness_spm(3) = sine/c_mps;
        delay_s         = xsensors_m * slowness_spm;
        
end
IM       = eye(M);
DK       = exp(-2j*pi*fK_Hz*delay_s');
K        = length(fK_Hz);
FIM_ik   = zeros(K,D,D);

for ik=1:K
    D_ik        = squeeze(DK(ik,:)) .';
    diagD_ik    = diag(D_ik);
    C_ik        = squeeze(C(ik,:,:));
    Gamma_ik    = diagD_ik*C_ik*diagD_ik'+sigma2noise*IM;
    rep3D_ik    = repmat(D_ik,1,D);
    d_ik        = (-2j*pi*fK_Hz(ik)*xsensors_m) .* rep3D_ik;
%     PiAortho = IM-D_ik*D_ik'/M;
%     totot = 2*real(d_ik'*PiAortho*d_ik)*real(D_ik'*(Gamma_ik\D_ik))/sigma2noise;
    for ell=1:D
        diagell = diag(d_ik(:,ell));
        dGamma_ik_ellplus = diagell *C_ik*diagD_ik';
        dGamma_ik_ellmoins = diagD_ik*C_ik* diagell';
        dGamma_ik_ell = dGamma_ik_ellplus+dGamma_ik_ellmoins;
        auxFIM        = Gamma_ik\dGamma_ik_ell;
        for ellp=1:D
            diagellp = diag(d_ik(:,ellp));
            dGamma_ik_ellpplus = diagellp *C_ik*diagD_ik';
            dGamma_ik_ellpmoins = diagD_ik*C_ik* diagellp';
            dGamma_ik_ellp = dGamma_ik_ellpplus+dGamma_ik_ellpmoins;
            auxFIMp = Gamma_ik\dGamma_ik_ellp';
            FIM_ik(ik,ell,ellp) = 0.5*real(trace(auxFIM * auxFIMp));
        end
    end
end
% factor 2 because we only add on K=N/2 frequencies dots
FIMslowness    = squeeze(sum(FIM_ik,1));
CRB.slowness   = pinv(FIMslowness);

% in 3D
% slowness_spm(1) = cosa*cose/c;
% slowness_spm(2) = sina*cose/c;
% slowness_spm(3) = sine/c;
switch D
    case 2
        Jacobav_k = ([...
            -sina/v_mps  -cosa/v_mps/v_mps; ...
            cosa/v_mps  -sina/v_mps/v_mps]);
        CRB.av =     (Jacobav_k\CRB.slowness)/Jacobav_k';
    case 3
        Jacobaec_k = ([...
            -sina*cose/c_mps -cosa*sine/c_mps -cosa*cose/c_mps/c_mps; ...
            cosa*cose/c_mps -sina*sine/c_mps -sina*cose/c_mps/c_mps;...
            0 cose/c_mps -sine/c_mps/c_mps]);
        CRB.aec =  (Jacobaec_k\CRB.slowness)/Jacobaec_k';
        
        % as v=c/cos(e)
        % slowness_spm(1) = cosa/v;
        % slowness_spm(2) = sina/v;
        % slowness_spm(3) = tane/v;
        
        Jacobav_k = ([...
            -sina/v_mps -cosa/v_mps/v_mps 0; ...
            cosa/v_mps -sina/v_mps/v_mps 0;...
            0 -sine/cose/v_mps/v_mps 1/cose/cose/v_mps ...
            ]);
        
        CRBaux = (Jacobav_k \ CRB.slowness) / Jacobav_k';
        CRB.av = CRBaux(1:2,1:2);
end

