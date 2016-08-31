function R = rndexp(varargin)
%RNDEXP Random matrices from an Exponential distribution
% 
% CALL:  R = rndexp(m0,sz)
%        R = rndexp(phat,sz) 
%     
%        R = matrix of random numbers
%       m0 = mean of the Exponential distribution, m0>0
%    phat = Distribution parameter struct
%             as returned from FITEXP.  
%       sz = size(R)    (Default size(m0))
%            sz can be a comma separated list or a vector 
%            giving the size of R (see zeros for options).
%
% The random numbers are generated by the inverse method. 
%
% Example:
%   R=rndexp(5,1,100);
%   phat=plotweib(R),shg     % Exp is Weibull with b=1
%
% See also  pdfexp, cdfexp, invexp,  fitexp, momexp

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on: Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - made sizing of R more flexible
% added ms 10.08.2000

error(nargchk(1,inf,nargin))
Np = 1;
options = []; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
% if numel(options)>1
%   error('Multidimensional struct of distribution parameter not allowed!')
% end

m0 = params{1};

if isempty(rndsize),
  csize=size(m0);
else
  [csize] = comnsize(m0,zeros(rndsize{:}));
end 
if any(isnan(csize))
   error('m0 must be a scalar or comply to the size info given.');
end

%R = invexp(rand(csize),m0);
%R = -log1p(-rand(csize)).*m0;
R = -log(rand(csize)).*m0;

