function [Gammas etas psis ts] = DNScylAxiFDsimulation(R,Z,parameters,initialConditions,boundaryInformation)
%NSCYLAXIFDSIMULATION does a complete simulation of the axisymmetric cylindrical NSE.
% [GAMMAS ETAS PSIS TS] = NSCYLAXIFDSIMULATION(R,Z,PARAMETERS,INITIALCONDITIONS,BOUNDARYINFORMATION)
% Does a complete simulation of the axisymmetric cylindrical NSE in real
% units as opposed to dimensionless units. The (z,r) grid is specified by
% the coordinate matrices R and Z (produced by meshgrid). The r-mesh can be
% non-uniformly spaced. The z-mesh can either be non-uniformly spaced or be
% periodic uniformly spaced. The fields Gamma, eta, and psi (angular
% velocity, theta component of vorticity, and the 2D streamfunction
% respectively) at each grid point over time are found using the
% second-order predictor-corrector scheme lined out in Appendix A (Lopez
% and Shen 1998). INITIALCONDITIONS is a structure holding the initial
% values of Gamma, eta, and psi in fields 'Gamma', 'eta', and 'psi'
% respectively. PARAMETERS is a structure holding the simulation parameters
% such as the number of time steps to simulate 'NumberTimeSteps', the time
% between each step 'dt', how often in steps should the fields be saved
% 'KeepEveryNSteps', and the kinematic viscosity 'nu'. BOUNDARYINFORMATION
% is a boundaryInformation structure that contains the boundary conditions
% and whose fields are described after this paragraph. The fields at the
% times steps that they are saved at (includes the initial and final
% fields) are returned in GAMMAS, ETAS, and PSIS with each page of the 3D
% matrices being a different time. The times of the saved fields are stored
% in the vector TS. Now for the fields of a boundaryInformation structure.
%
% periodicZ : Whether the system is axially periodic (true) or not (false)
%
% rMinBoundary, rMaxBoundary, zMinBoundary, zMaxBoundary : structures that
%   describe the boundaries at the minimum r-value, maximum r-value,
%   minimum z-value, and maximum z-value respectively. Each has a mandatory
%   field 'Type' which tells what type of boundary it it. It must be 'wall'
%   for all boundaries except the minimum r boundary where it can be 'axis'
%   if it is the axis of rotation (its r-value had better be zero or funny
%   results will occur). If it is a 'wall', then the field 'GammaValues'
%   must be a function handle that takes a single argument, the time T, and
%   returns the Gamma values that are the wall is supposed to have at each
%   grid point on the wall at that time.
%
%
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

warningState = warning('query','all');
warning('off','all');

% Apply the boundary conditions to the initial conditions in case they
% presently violate them.

[initialConditions.Gamma initialConditions.eta initialConditions.psi] = ...
            DNScylAxiFDapplyBoundaryConditions(initialConditions.Gamma, ...
            initialConditions.eta,initialConditions.psi,0,boundaryInformation);
        
% Make a matrix indicating which points are interior (true) and boundary
% (false). So, starting with a true matrix, the first and last r-values are
% set to false and if we are not axially periodic, the first and last
% z-value points are set to false.

interiorPoints = true(size(R));
interiorPoints(:,[1 end]) = false;
if ~boundaryInformation.periodicZ
    interiorPoints([1 end],:) = false;
end

% For the poisson equation in A.6 in the time evolution, we need a
% structure to hold various needed pieces of information. First, we should
% grab the interior points. Then, we need to get the CoefficientMatrix that
% gives CoefficientMatrix * psi(:) = R(:).*eta(:). Then, various things can
% be calculated from it to make it easier to solve that equation (the exact
% method will be decided in DNScylAxiFDdoTimeStep.m). The inverse will be
% gotten with the pseudoinverse (can't do normal inverse because it is
% singular). Then, the coefficient matrix will all rows and columns not in
% the interior (they are all empty rows and columns respectively) are
% stripped out producing a new coefficient matrix solely for the interior
% points. Then the matrix inverse and the incomplete LU factorization are
% computed for the interior matrix.

poissonEquation.interiorPoints = interiorPoints(:);
poissonEquation.CoefficientMatrix = DNScylAxiFDmakePoissonCoefficientMatrix(R,Z,boundaryInformation.periodicZ);
% poissonEquation.invCoefficientMatrix = pinv(full(poissonEquation.CoefficientMatrix));
poissonEquation.InteriorCoefficientMatrix = poissonEquation.CoefficientMatrix(interiorPoints,interiorPoints);
poissonEquation.invInteriorCoefficientMatrix = inv(poissonEquation.InteriorCoefficientMatrix);
% [poissonEquation.InteriorL poissonEquation.InteriorU] = luinc(poissonEquation.InteriorCoefficientMatrix,0);

% Make a vector to hold the times of the returned Gamma, eta, and psi
% fields. There will be an entry for the initial conditions, the end
% fields, and one every parameters.KeepEverNSteps.

ts = zeros(1,2 + floor(parameters.NumberTimeSteps/parameters.KeepEveryNSteps));

% Initialize Gammas, etas, and psis to hold all the Gamma, eta, and psi
% fields that will be saved (same number of rows and columns as R and one
% page for each one to be saved).

Gammas = zeros([size(R),numel(ts)]);
etas = Gammas;
psis = Gammas;

% Grab the initial conditions.

Gamma = initialConditions.Gamma;
eta = initialConditions.eta;
psi = initialConditions.psi;

% Put the initial conditions as the first entries of Gammas, etas, and
% psis.

Gammas(:,:,1) = Gamma;
etas(:,:,1) = eta;
psis(:,:,1) = psi;

% Do each time step.

for ii=1:parameters.NumberTimeSteps
    
    % Calculate the time at the beginning of the step.
    
    t = (ii-1)*parameters.dt;
    
    % Evolve the fluid by one time step.
    
    [Gamma eta psi] = DNScylAxiFDdoTimeStep(Gamma,eta,psi,t, ...
                parameters.dt,parameters.nu,R,Z, ...
                poissonEquation,boundaryInformation);
            
    % If this is one of the time steps at which we are to save the fields,
    % save the stuff to ts, Gammas, etas, and psis. dt must be added to t
    % because we are now at the end of the time step.
            
    if ii/parameters.KeepEveryNSteps == floor(ii/parameters.KeepEveryNSteps)
        disp([num2str(100*ii/parameters.NumberTimeSteps),'%']);
        ts(1 + ii/parameters.KeepEveryNSteps) = t+parameters.dt;
        Gammas(:,:,1 + ii/parameters.KeepEveryNSteps) = Gamma;
        etas(:,:,1 + ii/parameters.KeepEveryNSteps) = eta;
        psis(:,:,1 + ii/parameters.KeepEveryNSteps) = psi;
        
        % Check to see if any of the fields have nan's in them and if so,
        % terminate the simulation now and inform the user.
        
        if any(isnan(Gamma(:))) || any(isnan(eta(:))) || any(isnan(psi(:)))
            disp('NaN''s detected: simulation has become unstable.');
            ts = ts(1:(1 + ii/parameters.KeepEveryNSteps));
            Gammas = Gammas(:,:,1:(1 + ii/parameters.KeepEveryNSteps));
            etas = etas(:,:,1:(1 + ii/parameters.KeepEveryNSteps));
            psis = psis(:,:,1:(1 + ii/parameters.KeepEveryNSteps));
            break;
        end
        
    end
    
end

% Put the final conditions as the last entries of ts, Gammas, etas, and
% psis and we are done.

ts(end) = t+parameters.dt;
Gammas(:,:,end) = Gamma;
etas(:,:,end) = eta;
psis(:,:,end) = psi;

warning(warningState);
