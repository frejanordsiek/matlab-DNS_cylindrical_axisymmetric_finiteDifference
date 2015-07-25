function D2fDx2 = d2fdx2FiniteDifferenceNonUniform(f,x,dim)
%D2FDX2FINITEDIFFERENCENONUNIFORM does second finite difference on a non-uniform grid.
% D2FDX2 = D2FDX2FINITEDIFFERENCENONUNIFORM(F, X, DIM)
% Does the cecond centered finite difference on a non-uniform grid of F
% along direction DIM (1 is row, 2 is column, and 3 is page) where X is the
% vector containing the grid point positions in that direction. The formula
% is, for a one dimensional F (extends to two and three dimensions by just
% re-application of it to each row/col/page/whatever), given by
%
% D2f_i/Dx2 = [2/(x_(i+1) - x_(i-1))] * [(y_(i+1) - y_i)/(x_(i+1) - x_i) - (y_i - y_(i-1))/(x_i - x_(i-1))]
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

% We have to get the x_(i+1), x_i, x_(i-1), f_(i+1), f_i, and f(i-1) in the
% proper direction. To do this, x is first reshaped into a vector along the
% dimension dim and then repmatted to be the same size as f. Then x_plus,
% x_0, and x_minus (i+1, i, and i-1 respectively) are constructed in the
% appropriate directions. Then f_plus, f_0, and f_minus are contructed in
% the appropriate directions.

switch dim
    case 1
        
        if length(x) == numel(x)
            x = repmat(reshape(x,[numel(x),1,1]),[1,size(f,2),size(f,3)]);
        end
        
        x_plus = x(3:end,:,:);
        x_0 = x(2:(end-1),:,:);
        x_minus = x(1:(end-2),:,:);
        
        f_plus = f(3:end,:,:);
        f_0 = f(2:(end-1),:,:);
        f_minus = f(1:(end-2),:,:);
        
    case 2
        
        if length(x) == numel(x)
            x = repmat(reshape(x,[1,numel(x),1]),[size(f,1),1,size(f,3)]);
        end
        
        x_plus = x(:,3:end,:);
        x_0 = x(:,2:(end-1),:);
        x_minus = x(:,1:(end-2),:);
        
        f_plus = f(:,3:end,:);
        f_0 = f(:,2:(end-1),:);
        f_minus = f(:,1:(end-2),:);
        
    case 3
        
        if length(x) == numel(x)
            x = repmat(reshape(x,[1,1,numel(x)]),[size(f,1),size(f,2),1]);
        end
        
        x_plus = x(:,:,3:end);
        x_0 = x(:,:,2:(end-1));
        x_minus = x(:,:,1:(end-2));
        
        f_plus = f(:,:,3:end);
        f_0 = f(:,:,2:(end-1));
        f_minus = f(:,:,1:(end-2));
        
    otherwise
        error('DIM can only be 1, 2, or 3.');
end

% Now, calculate D2fDx2, which already has the appropriate direction stuff
% handled.

D2fDx2 = 2*((f_plus - f_0)./(x_plus - x_0) - (f_0 - f_minus)./(x_0 - x_minus))./(x_plus - x_minus);
