function [Sigma2_theta_s2pm2, logSkcp_pred, logresidue] = ...
    estimSigmatheta(Skcp, xsensor_m, fqnormrange, Fs_Hz, dimSlowness)
%======================================================
% Estimation of Sigma_theta in the gaussian model of
% loss of coherence:
% log-coherence propto (r_m-rk)'*Sigma_theta*(rm-rk)
% where Sigma_theta is 3 x 3 positive matrix.
% which is linear wrt the 6 entries of Sigma_theta.
% 
%
% Synopsis
% [Sigma_theta_spkm, logSkcp_pred, logresidue] = ...
%     estimSigmatheta(Skcp, xsensor_km, ...
%     fqnormrange, Fs_Hz, dimSlowness)
%======
% Inputs:
%    Skcp : C x N, coherence where C = M(M-1)/2
%           and M is the number of sensors
%           N is the number of frequency points
%    xsensor_m = M x 3 sensor locations
%    fqnormrange : sequence of normalized frequencies
%           between 0 and 1
%    Fs_Hz: sampling frequency in Hz
%    dimSlowness : 2 or 3, depending we take or not the altitude
%           of the sensors
%======
% Outputs:
%    Sigma_theta_spkm: positive matrix 3 x 3
%    logSkcp_pred : predicted log-coherence
%    logresidue: residue of prediction
%
%======================================================

combidim            = dimSlowness*(dimSlowness+1)/2;
M                   = size(xsensor_m,1);
combi2by2           = M*(M-1)/2;
difference_coord_cp = zeros(combi2by2,dimSlowness);
cp=0;
for i1=1:M-1
    for i2=i1+1:M
        cp         = cp+1;
        difference_coord_cp(cp,:) = ...
            xsensor_m(i2,1:dimSlowness)-xsensor_m(i1,1:dimSlowness);
    end
end
N           = length(fqnormrange);
regressors  = zeros(combi2by2*N,combidim);
logSkcp     = log(reshape(Skcp,combi2by2*N,1));

cq=0;
for i1=1:dimSlowness
    for i2=i1:dimSlowness
        cq   = cq+1;
        ps   = -4*pi*pi*difference_coord_cp(:,i1) .* ...
            difference_coord_cp(:,i2);
        aux  = ps*(fqnormrange .* fqnormrange)*Fs_Hz*Fs_Hz;
        regressors(:,cq) = reshape(aux,N*combi2by2,1);
    end
end
delta_theta_spm    = regressors\logSkcp;
delta_theta_spm_sq = zeros(dimSlowness);
cq=0;
for i1=1:dimSlowness
    for i2=i1:dimSlowness
        cq=cq+1;
        if i1==i2
            delta_theta_spm_sq(i1,i2) = delta_theta_spm(cq);
        else
            delta_theta_spm_sq(i1,i2) = delta_theta_spm(cq)/2;
        end
        delta_theta_spm_sq(i2,i1)     = delta_theta_spm_sq(i1,i2);
    end
end
Sigma2_theta_s2pm2   = delta_theta_spm_sq;
logSkcp_pred         = regressors*delta_theta_spm;
logresidue           = logSkcp-logSkcp_pred;
logSkcp_pred         = reshape(logSkcp_pred,combi2by2,N);
logresidue           = reshape(logresidue,combi2by2,N);
%======================================================
