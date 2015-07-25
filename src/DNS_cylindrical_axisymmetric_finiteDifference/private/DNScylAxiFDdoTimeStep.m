function [Gamma eta psi] = DNScylAxiFDdoTimeStep(Gamma,eta,psi,t,dt,nu,R,Z,poissonEquation,boundaryInformation)
%NSCYLAXIFDDOTIMESTEP does one time step in the axisymmetric cylindrical NSE.
% [GAMMA ETA PSI] = NSCYLAXIFDDOTIMESTEP(GAMMA,ETA,PSI,T,DT,NU,R,Z,POISSONEQUATION,BOUNDARYINFORMATION)
% Does one time step of time from time T to time T + DT in axisymmetric
% cylindrical Navier-Stokes equations using finite difference methods
% including applying the boundary conditions using the second-order
% predictor-corrector scheme lined out in Appendix A (Lopez and Shen 1998).
% GAMMA, ETA, and PSI are the angular velocity, theta component of
% vorticity, and the 2D streamfunction respectively at the grid points
% whose (z,r) coordinates are given by the coordinate matrices R and Z
% (produced by meshgrid). NU is the viscosity (note, this function works
% with the NSE in real units, not in dimensionless form).
% POISSONEQUATION is a structure containing stuff for solving the poisson
% equation in equation A.6 (Lopez and Shen 1998). It contains the
% PoissonCoefficientMatrix which solves
% PSI(:) = PoissonCoefficientMatrix * (R(:).*ETA(:)) along with other
% fields for solving the equation by methods other than doing a
% straightforward and slow and singular linear solution with mldivide.
% BOUNDARYINFORMATION is a boundaryInformation structure containing fields
% describing the boundary conditions.
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


% Doing step 1 in Appendix A to find Gamma_star and eta_star, which are
% Gamma and eta advanced in time by dt/2 in an Euler method fashion.

[Gamma_t eta_t] = DNScylAxiFDtimeDerivative(Gamma,eta,psi, ...
            nu,R,Z,boundaryInformation.periodicZ);
        
Gamma_star = Gamma + (0.5*dt)*Gamma_t;
eta_star = eta + (0.5*dt)*eta_t;

% Find psi_star (step 2 in Appendix A) given the values of eta_star and the
% information in poissonEquation describing equation A.6. There are
% multiple ways the equation can be solved and many are present here with
% all but one commented out. The first method is to just use MATLAB's built
% in linear equation solver with the original coefficient matrix but this
% is problematic since it is singular. The second method is to do a matrix
% multiplication with the pseudo inverse of the original coefficient matrix
% which avoids the issue of being singular but this is slow. The other
% methods work solely with the interior points so psi_star must first be
% initialized to zero. The third and fourth methods are just the first and
% second methods but for the interior points and thus avoids singular
% behavior and are a lot faster. The fifth method uses the conjugate
% gradients iterative method using the incomplete LU factorization as
% preconditioners.

% psi_star = reshape(poissonEquation.CoefficientMatrix\(R(:).*eta_star(:)),size(Gamma));
% psi_star = reshape(poissonEquation.invCoefficientMatrix*(R(:).*eta_star(:)),size(Gamma));

psi_star = zeros(size(Gamma));

% psi_star(poissonEquation.interiorPoints) = poissonEquation.InteriorCoefficientMatrix ...
%             \ (R(poissonEquation.interiorPoints).*eta_star(poissonEquation.interiorPoints));
        
psi_star(poissonEquation.interiorPoints) = full(poissonEquation.invInteriorCoefficientMatrix ...
            * (R(poissonEquation.interiorPoints).*eta_star(poissonEquation.interiorPoints)));
        
% [psi_star(poissonEquation.interiorPoints)  flag relres iter resvec] = ...
%             cgs(poissonEquation.InteriorCoefficientMatrix, ...
%             (R(poissonEquation.interiorPoints).*eta_star(poissonEquation.interiorPoints)), ...
%             1e-16,100, ...
%             poissonEquation.InteriorL,poissonEquation.InteriorU, ...
%             psi(poissonEquation.interiorPoints));

% Apply the boundary conditions (step 3).

[Gamma_star eta_star psi_star] = DNScylAxiFDapplyBoundaryConditions(Gamma_star, ...
            eta_star, psi_star, t, boundaryInformation);

% Doing step 4 in Appendix A which is to advance Gamma and eta by dt using
% the derivatives gotten using star'ed fields.

[Gamma_t eta_t] = DNScylAxiFDtimeDerivative(Gamma_star,eta_star,psi_star, ...
            nu,R,Z,boundaryInformation.periodicZ);
        
Gamma = Gamma + dt*Gamma_t;
eta = eta + dt*eta_t;

% Find psi_star (step 2 in Appendix A) given the values of eta_star and the
% information in poissonEquation describing equation A.6. There are
% multiple ways the equation can be solved and many are present here with
% all but one commented out. The first method is to just use MATLAB's built
% in linear equation solver with the original coefficient matrix but this
% is problematic since it is singular. The second method is to do a matrix
% multiplication with the pseudo inverse of the original coefficient matrix
% which avoids the issue of being singular but this is slow. The other
% methods work solely with the interior points so psi_star must first be
% initialized to zero. The third and fourth methods are just the first and
% second methods but for the interior points and thus avoids singular
% behavior and are a lot faster. The fifth method uses the conjugate
% gradients iterative method using the incomplete LU factorization as
% preconditioners.

% psi = reshape(poissonEquation.CoefficientMatrix\(R(:).*eta(:)),size(Gamma));
% psi = reshape(poissonEquation.invCoefficientMatrix*(R(:).*eta(:)),size(Gamma));

psi = zeros(size(Gamma));

% psi(poissonEquation.interiorPoints) = poissonEquation.InteriorCoefficientMatrix ...
%             \ (R(poissonEquation.interiorPoints).*eta(poissonEquation.interiorPoints));
        
psi(poissonEquation.interiorPoints) = full(poissonEquation.invInteriorCoefficientMatrix ...
            * (R(poissonEquation.interiorPoints).*eta(poissonEquation.interiorPoints)));
        
% [psi(poissonEquation.interiorPoints) flag relres iter resvec] =  ...
%             cgs(poissonEquation.InteriorCoefficientMatrix, ...
%             (R(poissonEquation.interiorPoints).*eta(poissonEquation.interiorPoints)), ...
%             1e-16,100, ...
%             poissonEquation.InteriorL,poissonEquation.InteriorU, ...
%             psi_star(poissonEquation.interiorPoints));

% Apply the boundary conditions (step 6).

[Gamma eta psi] = DNScylAxiFDapplyBoundaryConditions(Gamma, ...
            eta, psi, t, boundaryInformation);
