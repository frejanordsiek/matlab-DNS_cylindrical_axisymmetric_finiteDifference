function [Gamma_t eta_t] = DNScylAxiFDtimeDerivative(Gamma,eta,psi,nu,R,Z,periodicZ)
%NSCYLAXIFDTIMEDERIVATIVE finds the time derivatives in the axisymmetric cylindrical NSE.
% [GAMMA_t ETA_t] = NSCYLAXIFDTIMEDERIVATIVE(GAMMA,ETA,PSI,NU,R,Z,PERIODICZ)
% Finds the time derivatives of GAMMA and ETA in axisymmetric cylindrical
% Navier-Stokes equations using finite difference methods, but does not
% apply the boundary conditions. GAMMA, ETA, and PSI are the angular
% velocity, theta component of vorticity, and the 2D streamfunction
% respectively at the grid points whose (z,r) coordinates are given by the
% coordinate matrices R and Z (produced by meshgrid). NU is the viscosity
% (note, this function works with the NSE in real units, not in
% dimensionless form). PERIODICZ is a flag indicating whether the axial
% boundary condition is periodic (true) or not (false).
%
% References:
%   1)  JM Lopez and J Shen. An Efficient Spectral-Projection Method for
%       the Navier-Stokes Equations in Cylindrical Geometries - Part I -
%       Axisymmetric Cases. Journal of Computational Physics. 139, 308-326
%       (1998)
%
%
%
% Copyright 2012 Freja Nordsiek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Since the eta equation uses a lot of derivatives of eta/r, it is easier
% and more efficient to just calculate it first instead of doing it each
% time individually.

Rinv = R.^-1;

etaOverR = eta.*Rinv;

% Do all the derivatives with respect to r that are needed making sure to
% use the non-uniform grid finite difference functions.

Gamma_r = dfdxFiniteDifferenceNonUniform(Gamma,R,2);
Gamma_rr = d2fdx2FiniteDifferenceNonUniform(Gamma,R,2);
psi_r = dfdxFiniteDifferenceNonUniform(psi,R,2);
% psi_rr = d2fdx2FiniteDifferenceNonUniform(psi,R,2);
etaOverR_r = dfdxFiniteDifferenceNonUniform(etaOverR,R,2);
etaOverR_rr = d2fdx2FiniteDifferenceNonUniform(etaOverR,R,2);

% How the axial derivatives are handled depends on whether we are doing
% periodic boundary conditions or not. If it is periodic, then the periodic
% uniform grid finite difference codes are used. Otherwise, the non-uniform
% grid codes are used.

if periodicZ
    dz = diff(Z(1:2,1));
    Gamma_z = dfdxFiniteDifferencePeriodic(Gamma,dz,1);
    Gamma_zz = d2fdx2FiniteDifferencePeriodic(Gamma,dz,1);
    psi_z = dfdxFiniteDifferencePeriodic(psi,dz,1);
%     psi_zz = d2fdx2FiniteDifferencePeriodic(psi,dz,1);
    etaOverR_z = dfdxFiniteDifferencePeriodic(etaOverR,dz,1);
    etaOverR_zz = d2fdx2FiniteDifferencePeriodic(etaOverR,dz,1);
else
    Gamma_z = dfdxFiniteDifferenceNonUniform(Gamma,Z,1);
    Gamma_zz = d2fdx2FiniteDifferenceNonUniform(Gamma,Z,1);
    psi_z = dfdxFiniteDifferenceNonUniform(psi,Z,1);
%     psi_zz = d2fdx2FiniteDifferenceNonUniform(psi,Z,1);
    etaOverR_z = dfdxFiniteDifferenceNonUniform(etaOverR,Z,1);
    etaOverR_zz = d2fdx2FiniteDifferenceNonUniform(etaOverR,Z,1);
end

% We want all the derivative matrices to be the same size as Gamma, eta,
% and psi; so we have to zero pad in a column of zeros before and after for
% all the r-derivative matrices since the finite difference codes only give
% the derivatives for the interior r-values.

pad = zeros(size(R,1),1);

Gamma_r = [pad,Gamma_r,pad];
Gamma_rr = [pad,Gamma_rr,pad];
psi_r = [pad,psi_r,pad];
% psi_rr = [pad,psi_rr,pad];
etaOverR_r = [pad,etaOverR_r,pad];
etaOverR_rr = [pad,etaOverR_rr,pad];

% If we are not axially periodic, then we need to pad a row of zeros before
% and after the z-derivative matrices (the finite difference codes only
% give derivatives for interior z-values) and zero out the top and bottom
% rows of the r-derivatives so as to not effect the boundary points.

if ~periodicZ
    
    pad = zeros(1,size(Z,2));
    Gamma_z = [pad;Gamma_z;pad];
    Gamma_zz = [pad;Gamma_zz;pad];
    psi_z = [pad;psi_z;pad];
%     psi_zz = [pad;psi_zz;pad];
    etaOverR_z = [pad;etaOverR_z;pad];
    etaOverR_zz = [pad;etaOverR_zz;pad];
    
    Gamma_r([1 end],:) = 0;
    Gamma_rr([1 end],:) = 0;
    psi_r([1 end],:) = 0;
%     psi_rr([1 end],:) = 0;
    etaOverR_r([1 end],:) = 0;
    etaOverR_rr([1 end],:) = 0;
    
end

% Calculate the time derivatives of Gamma and eta. Equations A.4 and A.5,
% after writing out the terms, going back to non-normalized units (also
% means 1/Re is replaced with nu), and a bit of simplifying leads to these
% equations.

Gamma_t = nu*(Gamma_rr + Gamma_zz) + ((psi_z - nu).*Gamma_r - psi_r.*Gamma_z).*Rinv;

eta_t = nu*R.*(etaOverR_rr + etaOverR_zz) + (psi_z + nu).*etaOverR_r ...
            - psi_r.*etaOverR_z + 2*Gamma.*Gamma_z.*(Rinv.^3);
