function [u v w xsi eta zeta] = DNScylAxiFDcalcVelocity(R,Z,Gamma,eta,psi,periodicZ)
%NSCYLAXIFDCALCVELOCITY calculates the velocity and vorticity from the (Gamma,eta,psi) in the axisymmetric cylindrical NSE.
% [U V W XSI ETA ZETA] = NSCYLAXIFDCALCVELOCITY(R,Z,GAMMA,ETA,PSI,PERIODICZ)
% Calculates the fluid velocity (U,V,W) and vorticity (XSI,ETA,ZETA) fields
% (r,theta,z order) from the GAMMA, ETA, PSI fields (angular velocity,
% theta component of vorticity, and the 2D streamfunction respectively) at
% all the points on the (z,r) grid specified by the coordinate matrices R
% and Z (produced by meshgrid). The r-mesh can be non-uniformly spaced. The
% z-mesh can either be non-uniformly spaced or be periodic uniformly
% spaced. PERIODICZ indicates hether the system is axially periodic (true)
% or not (false)
%
% References:
%   1)  JM Lopez and J Shen. An Efficient Spectral-Projection Method for
%       the Navier-Stokes Equations in Cylindrical Geometries - Part I -
%       Axisymmetric Cases. Journal of Computational Physics. 139, 308-326
%       (1998)
%       DOI: 10.1006/jcph.1997.5872
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

% Calculate the various r-derivatives.

Gamma_r = dfdxFiniteDifferenceNonUniform(Gamma,R,2);
psi_r = dfdxFiniteDifferenceNonUniform(psi,R,2);
psi_rr = d2fdx2FiniteDifferenceNonUniform(psi,R,2);

% Calculate the various z-derivatives. The method of doing so depends on
% whether the system is periodic or not.

if periodicZ
    dz = diff(Z(1:2,1));
    Gamma_z = dfdxFiniteDifferencePeriodic(Gamma,dz,1);
    psi_z = dfdxFiniteDifferencePeriodic(psi,dz,1);
    psi_zz = d2fdx2FiniteDifferencePeriodic(psi,dz,1);
else
    Gamma_z = dfdxFiniteDifferenceNonUniform(Gamma,Z,1);
    psi_z = dfdxFiniteDifferenceNonUniform(psi,Z,1);
    psi_zz = d2fdx2FiniteDifferenceNonUniform(psi,Z,1);
end

% We want all the derivative matrices to be the same size as Gamma, eta,
% and psi; so we have to zero pad in a column of zeros before and after for
% all the r-derivative matrices since the finite difference codes only give
% the derivatives for the interior r-values.

pad = zeros(size(R,1),1);
Gamma_r = [pad,Gamma_r,pad];
psi_r = [pad,psi_r,pad];
psi_rr = [pad,psi_rr,pad];

% If we are not axially periodic, then we need to pad a row of zeros before
% and after the z-derivative matrices (the finite difference codes only
% give derivatives for interior z-values) and zero out the top and bottom
% rows of the r-derivatives so as to not effect the boundary points.

if ~periodicZ
    
    pad = zeros(1,size(Z,2));
    Gamma_z = [pad;Gamma_z;pad];
    psi_z = [pad;psi_z;pad];
    psi_zz = [pad;psi_zz;pad];
    
    Gamma_r([1 end],:) = 0;
    psi_r([1 end],:) = 0;
    psi_rr([1 end],:) = 0;
    
end

% Calculate the velocity and vorticity fields.

u = -psi_z./R;
v = Gamma./R;
w = psi_r./R;

xsi = -Gamma_z./R;
eta = -(psi_rr - psi_r./R + psi_zz)./R;
zeta = Gamma_r./R;
