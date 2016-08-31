function R = rndgenpar(varargin)
%RNDGENPAR Random matrices from a Generalized Pareto Distribution
% 
% CALL:  R = rndgenpar(k,s,m,sz,options);
%        R = rndgenpar(phat,sz,options);
%     
%        R = matrix of random numbers
%        k = shape parameter in the GPD  (see cdfgenpar)
%        s = scale parameter in the GPD    (default 1)
%        m = location parameter in the GPD (default 0)
%     phat = Distribution parameter struct
%             as returned from FITGENPAR.  
%       sz = size(R)    (Default common size of k, s and m)
%             sz can be a comma separated list or a vector 
%             giving the size of R (see zeros for options).
%   options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%
% The random numbers are generated by the inverse method. 
%
% Example:
%   R1=rndgenpar(2,1,0,1,100);  % GPD k=2
%   R2=rndgenpar(1,1,0,1,100);  % GPD k=1  ==>  Uniform
%   R3=rndgenpar(0,1,0,1,100);  % GPD k=0  ==>  Exponential
%   plot([R1 R2 R3],'.')
%
% See also  pdfgenpar, cdfgenpar, invgenpar, fitgenpar, momgenpar

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


% Tested on: Matlab 5.3
% History: 
% Revised by jr 22.12.1999
% revised ms 14.06.2000
% - updated header info
% - changed name to rndgenpar (from gpdrnd)
% - allowed 2 arguments
% revised pab 23.10.2000
%   - added default s,m0
%  - added comnsize, nargchk
%  - added greater flexibility on the sizing of R
% -revised pab sept 2007
% - more accurate simulation of upper tail.

error(nargchk(1,inf,nargin))
Np = 3;
options = struct('lowertail',false); 
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[k,s,m0] = deal(params{:});
if isempty(s),  s=1;end
if isempty(m0), m0=0;end

if isempty(rndsize)
  csize = comnsize(k,s,m0);
else
  csize = comnsize(k,s,m0, zeros(rndsize{:}));
end

if any(isnan(csize))
  error('k,s and m0 must be of common size or scalar.');
end
R = invgenpar(rand(csize),k,s,m0,options);

