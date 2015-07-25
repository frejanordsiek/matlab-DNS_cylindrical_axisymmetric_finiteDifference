function D2fDx2 = d2fdx2FiniteDifference(f,h,dim)
%D2FDX2FINITEDIFFERENCE does second finite difference on a uniform grid.
% D2FDX2 = D2FDX2FINITEDIFFERENCE(F, H, DIM)
% Does the second centered finite difference on a uniform grid of F
% along direction DIM (1 is row, 2 is column, and 3 is page) where H is the
% grid spacing. The formula is, for a one dimensional F (extends to two and
% three dimensions by just re-application of it to each
% row/col/page/whatever), given by
%
% D2f_i/Dx2 = (f_(i+1) - 2*f_i + f_(i-1)) / h^2
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

% f_plus, f_0, and f_minus (i+1, i, and i-1 respectively) are contructed in
% the appropriate directions.

switch dim
    case 1
        f_plus = f(3:end,:,:);
        f_0 = f(2:(end-1),:,:);
        f_minus = f(1:(end-2),:,:);     
    case 2
        f_plus = f(:,3:end,:);
        f_0 = f(:,2:(end-1),:);
        f_minus = f(:,1:(end-2),:);  
    case 3
        f_plus = f(:,:,3:end);
        f_0 = f(:,:,2:(end-1));
        f_minus = f(:,:,1:(end-2));
    otherwise
        error('DIM can only be 1, 2, or 3.');
end

% Now, calculate DfDx, which already has the appropriate direction stuff
% handled.

D2fDx2 = (f_plus - 2*f_0 + f_minus)*(1/h^2);
