function palm_core(varargin)
% This is the core PALM function.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2014
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Uncomment the line below for debugging:
%clear global plm opts; global plm opts;

% Take the arguments. Save a small log if needed.
ticI = tic;
[opts,plm] = palm_takeargs(varargin{:});
tocI = toc(ticI);
fprintf('Elapsed time parsing inputs: ~ %g seconds.\n',tocI);  
    
if opts.cca.do
    [opts,plm] = palm_prepcca(opts,plm)   % <- (new code)
    [opts,plm] = palm_permcca(opts,plm)   % <- (new code)
    palm_savecca(opts,plm)                % <- (new code)
else
    ticP = tic;
    [opts,plm] = palm_prepglm(opts,plm);
    [opts,plm] = palm_glm(opts,plm);      
    tocP = toc(ticP);
    fprintf('Elapsed time with permutations: ~ %g seconds.\n',tocP);
    
    % Save everything, except the few bits saved above in palm_glm
    ticS = tic;
    palm_saveglm(plm,opts);
    tocS = toc(ticS);                     
    fprintf('Elapsed time generating and saving results: ~ %g seconds.\n',tocS);
end

% Finished. :-)
fprintf('Overall elapsed time: ~ %g seconds.\n',tocI+tocP+tocS);
csvwrite(sprintf('%s_elapsed.csv',opts.o),[tocI tocP tocS]);              
fprintf('PALM finished at %s.\n',datestr(now));

