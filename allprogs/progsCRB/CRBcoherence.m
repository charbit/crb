function CRB = CRBcoherence(xsensor, sigma2noise, C, aec, K, Fs_Hz)

M    = size(xsensor,1);

cosa = cos(aec.a_deg*pi/180);
sina = sin(aec.a_deg*pi/180);
cose = cos(aec.e_deg*pi/180);
sine = sin(aec.e_deg*pi/180);
c    = aec.c_mps;

slowness_spm    = zeros(3,1);
slowness_spm(1) = cosa*cose/c;
slowness_spm(2) = sina*cose/c;
slowness_spm(3) = sine/c;
delay_s         = xsensor * slowness_spm;

fK       = (1:K)'*Fs_Hz/2/K;
IM       = eye(M);
DK       = exp(-2j*pi*fK*delay_s');

FIM_ik   = zeros(K,3,3);

for ik=1:K
    D_ik        = squeeze(DK(ik,:)) .';
    diagD_ik    = diag(D_ik);
    C_ik        = squeeze(C(ik,:,:));
    Gamma_ik    = diagD_ik*C_ik*diagD_ik'+sigma2noise*IM;
    rep3D_ik = repmat(D_ik,1,3);
    d_ik = (-2j*pi*fK(ik)*xsensor) .* rep3D_ik;
%     PiAortho = IM-D_ik*D_ik'/M;
%     totot = 2*real(d_ik'*PiAortho*d_ik)*real(D_ik'*(Gamma_ik\D_ik))/sigma2noise;
    for ell=1:3
        diagell = diag(d_ik(:,ell));
        dGamma_ik_ellplus = diagell *C_ik*diagD_ik';
        dGamma_ik_ellmoins = diagD_ik*C_ik* diagell';
        dGamma_ik_ell = dGamma_ik_ellplus+dGamma_ik_ellmoins;
        for ellp=1:3
            diagellp = diag(d_ik(:,ellp));
            dGamma_ik_ellpplus = diagellp *C_ik*diagD_ik';
            dGamma_ik_ellpmoins = diagD_ik*C_ik* diagellp';
            dGamma_ik_ellp = dGamma_ik_ellpplus+dGamma_ik_ellpmoins;
            auxFIM  = Gamma_ik\dGamma_ik_ell;
            auxFIMp = Gamma_ik\dGamma_ik_ellp';
            FIM_ik(ik,ell,ellp) = real(trace(auxFIM * auxFIMp));
        end
    end
end
% FIM_ik       = real(FIM_ik);
FIM.slowness = squeeze(sum(FIM_ik,1));
CRB.slowness = pinv(FIM.slowness);

% slowness_spm(1) = cosa*cose/c;
% slowness_spm(2) = sina*cose/c;
% slowness_spm(3) = sine/c;

Jacobaec_k = ([...
    -sina*cose/c -cosa*sine/c -cosa*cose/c/c; ...
    cosa*cose/c -sina*sine/c sina*cose/c/c;...
    0 cose/c sine/c/c]);
FIM.aec = Jacobaec_k*FIM.slowness*Jacobaec_k';
CRB.aec = pinv(FIM.aec);
