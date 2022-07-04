function [opts,plm] = palm_takeargs(varargin)
% Handle the inputs for PALM.
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
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Load the defaults
opts = palm_defaults;

% As varargin is actually from another function, fix it.
if nargin == 1
    if exist(varargin{1},'file')
        vararginx = palm_configrw(varargin{1});
    else
        error('Unknown option or file not found: %s',varargin{1});
    end
else
    vararginx = varargin;
    idxa = find(strcmpi(vararginx,'-o'));
    if isempty(idxa)
        otmp = opts.o;
    else
        otmp = vararginx{idxa+1};
    end
    if ~ strcmp(otmp(end),'_')
        otmp = horzcat(otmp,'_');
    end
    cfgname = horzcat(otmp,'palmconfig.txt');
    [opth,~,~] = fileparts(cfgname);
    if ~isempty(opth) && ~exist(opth,'dir')
        mkdir(opth);
    end
    palm_configrw(vararginx,cfgname);
end

% Number of input images/masks/surfaces
% These are NOT meant to be edited.
Ni     = sum(strcmp(vararginx,'-i'));        % number of data inputs
Nx     = sum(strcmp(vararginx,'-x'));        % number of data inputs (CCA left side)
Ny     = sum(strcmp(vararginx,'-y'));        % number of data inputs (CCA right side)
Nz     = sum(strcmp(vararginx,'-z'));        % number of nuisance inputs (CCA left side or both)
Nw     = sum(strcmp(vararginx,'-w'));        % number of nuisance inputs (CCA right side)
Nm     = sum(strcmp(vararginx,'-m'));        % number of masks
Ns     = sum(strcmp(vararginx,'-s'));        % number of surfaces
Nd     = sum(strcmp(vararginx,'-d'));        % number of design files
Nt     = sum(strcmp(vararginx,'-t'));        % number of t-contrast files
Nf     = sum(strcmp(vararginx,'-f'));        % number of F-test files
Ncon   = sum(strcmp(vararginx,'-con'));      % number of contrast files (t or F, mset format)
Nevd   = sum(strcmp(vararginx,'-evperdat')); % number of EV per datum inputs
Nzevd   = sum(strcmp(vararginx,'-zevperdat')); % number of ZEV per datum inputs (for CCA)
Nwevd   = sum(strcmp(vararginx,'-wevperdat')); % number of WEV per datum inputs (for CCA)
Nimiss = sum(strcmp(vararginx,'-imiss'));    % number of missing indicators for inputs
Ndmiss = sum(strcmp(vararginx,'-dmiss'));    % number of missing indicators for designs
opts.i      = cell(Ni,1);   % Input files (to constitute Y later)
opts.x      = cell(Nx,1);   % Input files (to constitute X for CCA later)
opts.y      = cell(Ny,1);   % Input files (to constitute Y for CCA later)
opts.z      = cell(Nz,1);   % Input files (to constitute Z for CCA later)
opts.w      = cell(Nw,1);   % Input files (to constitute W for CCA later)
opts.m      = cell(Nm,1);   % Mask file(s)
opts.s      = cell(Ns,1);   % Surface file(s)
opts.sa     = cell(Ns,1);   % Area file(s) or weight(s)
opts.d      = cell(Nd,1);   % Design file(s)
opts.imiss  = cell(Nd,1);   % Design file(s)
opts.dmiss  = cell(Nd,1);   % Design file(s)
opts.t      = cell(Nt,1);   % t contrast file(s)
opts.f      = opts.t;       % F contrast file(s)
opts.Ccon   = cell(Ncon,1); % Contrast file(s) (t or F, mset format)
opts.Dcon   = cell(Ncon,1); % Contrast file(s) (multivariate, mset format)
opts.eb       = [];       % File with definition of exchangeability blocks
opts.vg       = [];       % File with definition of variance groups
opts.EE       = false;    % To be filled below (don't edit this!)
opts.ISE      = false;    % To be filled below (don't edit this!)
opts.within   = false;    % To be filled below (don't edit this!)
opts.whole    = false;    % To be filled below (don't edit this!)
opts.conskipcount = 0;    % When saving the contrasts, skip how many from 1?
opts.singlevg = true;     % Make sure that sigle VG will be used if nothing is supplied (this is NOT a "default" setting, and it's not a setting at all, but hard coded. Don't edit it!)
opts.subjidx  = [];       % Filename of the indices of subjects to keep
opts.Nimiss   = Nimiss;   % Store Nimiss in orer pass it along to palm_prepglm.m
opts.Ndmiss   = Ndmiss;   % Store Ndmiss in orer pass it along to palm_prepglm.m
plm.subjidx   = [];       % Indices of subjects to keep
plm.nEVdat    = Nevd;     % Store number of EV per datum
plm.nZEVdat    = Nzevd;     % Store number of ZEV per datum (for CCA)
plm.nWEVdat    = Nwevd;     % Store number of WEV per datum (for CCA)

% These are to be incremented below
i = 1; m = 1; d = 1;
t = 1; s = 1;
con = 1; ev = 1;
imiss = 1; dmiss = 1;
x = 1; y = 1; z = 1; w = 1; zev = 1; wev = 1; % For CCA

% Remove trailing empty arguments. This is useful for some Octave versions.
while numel(vararginx) > 0 && isempty(vararginx{1})
    vararginx(1) = [];
end
narginx = numel(vararginx);

% Take the input arguments
a = 1;
while a <= narginx
    switch vararginx{a}
        case {'-help','-?','-basic','-advanced'}
            
            % Do nothing, as these options are parsed separately,
            % and should anyway be given without any other argument.
            a = a + 1;
            
        case '-i' % basic
            
            % Get the filenames for the data.
            opts.i{i} = vararginx{a+1};
            i = i + 1;
            a = a + 2;
            
            % If -i is provided, set CCA, PLS, CACIC to false
            opts.cca.do   = false;
            opts.pls.do   = false;
            opts.cacic.do = false; 
            
        case '-cca' % advanced
            
            % Doing CCA?
            opts.cca.do   = true;
            opts.pls.do   = false;
            opts.cacic.do = false;
            
            if nargin == a || ...
                    ischar(vararginx{a+1}) && ...
                    strcmpi(vararginx{a+1}(1),'-')
                opts.ccaorplsparm = 1;
                a = a + 1;
            elseif ischar(vararginx{a+1})
                opts.ccaorplsparm = eval(vararginx{a+1});
                a = a + 2;
            else
                opts.ccaorplsparm = vararginx{a+1};
                a = a + 2;
            end
            
            % Save additional parameters
            opts.idxout = true;
       
        case '-pls' % advanced
            
            % Doing PLS?
            opts.cca.do   = false;
            opts.pls.do   = true;
            opts.cacic.do = false;
            a = a + 1;
            
        case '-cacic' % advanced
            
            % Doing cacic?
            opts.cca.do   = false;
            opts.pls.do   = false;
            opts.cacic.do = true;
            a = a + 1;
            
        case '-x' % basic
            
            % Get the filenames for X (CCA)
            opts.x{x} = vararginx{a+1};
            x = x + 1;
            a = a + 2;
            
        case '-y' % basic
            
            % Get the filenames for the data.
            opts.y{y} = vararginx{a+1};
            y = y + 1;
            a = a + 2;   
            
        case '-z' % basic
            if z > 1
                error('At most one csv file for z is allowed. To supply vertex or voxelwise nuisance vars use -zevperdat')
            end
            
            % Get the filenames for nuisance (CCA left side or both)
            opts.z{z} = vararginx{a+1};
            z = z + 1;
            a = a + 2;
            
        case '-w' % basic
            if w > 1
                error('At most one csv file for w is allowed. To supply vertex or voxelwise nuisance vars use -wevperdat')
            end
            
            % Get the filenames for nuisance (CCA right side)
            opts.w{w} = vararginx{a+1};
            w = w + 1;
            a = a + 2;
            
        case {'-semipartial','-part'} % advanced
            
            if nargin == a || nargin > a && strcmp(vararginx{a+1}(1),'-')
                error('If using semipartial CCA, you must specify a side to adjust (i.e. left or right)');
               
            elseif nargin > a
                
                % Which side to adjust for nuisance variables Z?
                sidelist = {            ...
                    'left',             ...
                    'y',                ...
                    'right',            ...
                    'x'};
                
                sideidx = strcmpi(vararginx{a+1},sidelist);
                
                % Check if a valid side is specified, and raise error if not
                if ~any(sideidx)
                    error('For semipartial CCA, specified side must be either left (y), or right (x)');
                elseif strcmpi(vararginx{a+1},'left') || strcmpi(vararginx{a+1},'y')                 
                    opts.cca.semipartial.side  = 'left';
                    a = a + 2;
                    
                elseif strcmpi(vararginx{a+1},'right') || strcmpi(vararginx{a+1},'x')                
                    opts.cca.semipartial.side  = 'right';
                    a = a + 2;
                end
            end
       
        case '-theil' % advanced
            
            if nargin == a || nargin > a && strcmp(vararginx{a+1}(1),'-')
                opts.cca.theil.selectionmatrix=false;
                a = a + 1;
            else
                % Get the filename for the selection matrix (csv file)
                opts.cca.theil.selectionmatrix{1} = vararginx{a+1};
                a = a + 2;
            end
            
        case '-zevperdat' % advanced
            
            % Use one ZEV per datum?
            opts.zevperdat = true;
            opts.zevdatfile{zev} = vararginx{a+1};
            zev = zev + 1;
            a = a + 2;
            
        case '-wevperdat' % advanced
            
            % Use one ZEV per datum?
            opts.wevperdat = true;
            opts.wevdatfile{wev} = vararginx{a+1};
            wev = wev + 1;
            a = a + 2;
                 
        case '-m' % basic
            
            % Get the filenames for the masks, if any.
            opts.m{m} = vararginx{a+1};
            m = m + 1;
            a = a + 2;

        case {'-s','-surf'} % basic
            
            % Get the filenames for the surfaces, if any.
            opts.s{s}  = vararginx{a+1};
            if nargin == a+1 || (narginx>a+1 && strcmp(vararginx{a+2}(1),'-'))
                opts.sa{s} = [];
                a = a + 2;
            else
                opts.sa{s} = vararginx{a+2};
                a = a + 3;
            end
            s = s + 1;
            
        case '-d' % basic
            
            % Get the design matrix file.
            opts.d{d} = vararginx{a+1};
            d = d + 1;
            a = a + 2;
            
        case '-evperdat' % advanced
            
            % Use one EV per datum?
            opts.evperdat = true;
            opts.evdatfile{ev} = vararginx{a+1};
            if nargin == a + 1 || ...
                    ischar(vararginx{a+2}) && ...
                    strcmpi(vararginx{a+2}(1),'-')
                opts.evpos{ev}(1) = 1;  % EV position
                opts.evpos{ev}(2) = 1;  % Design number
                a = a + 2;
            elseif nargin == a + 2 || ...
                    ischar(vararginx{a+3}) && ...
                    strcmpi(vararginx{a+3}(1),'-')
                if ischar(vararginx{a+2}) % EV position
                    opts.evpos{ev}(1) = eval(vararginx{a+2});
                else
                    opts.evpos{ev}(1) = vararginx{a+2};
                end
                opts.evpos{ev}(2) = 1;  % Design number
                a = a + 3;
            else
                if ischar(vararginx{a+2}) % EV position
                    opts.evpos{ev}(1) = eval(vararginx{a+2});
                else
                    opts.evpos{ev}(1) = vararginx{a+2};
                end
                if ischar(vararginx{a+3}) % Design number
                    opts.evpos{ev}(2) = eval(vararginx{a+3});
                else
                    opts.evpos{ev}(2) = vararginx{a+3};
                end
                a = a + 4;
            end
            ev = ev + 1;
            
        case '-imiss' % basic
            
            % Get the filenames for the missing data indicators (inputs).
            opts.imiss{imiss} = vararginx{a+1};
            imiss = imiss + 1;
            a = a + 2;
            
        case '-dmiss' % basic
            
            % Get the filenames for the missing data indicators (designs).
            opts.dmiss{dmiss} = vararginx{a+1};
            dmiss = dmiss + 1;
            a = a + 2;
            
        case '-mcar'
            
            % For the missing data, treat as missing completely at random.
            opts.mcar = true;
            a = a + 1;
            
        case '-t' % basic
            
            % Get the t contrast files.
            opts.t{t} = vararginx{a+1};
            t = t + 1;
            a = a + 2;
            
        case '-f' % basic
            
            % Get the F contrast files.
            if t == 1
                error('The option "-f" cannot be specified before its respective "-t".');
            end
            opts.f{t-1} = vararginx{a+1};
            a = a + 2;
            
        case '-con' % advanced
            
            % Get the contrast files from an .mset file or
            % pair of files. If a pair, the 1st is for Cset
            % and the second for Dset.
            opts.Ccon{con} = vararginx{a+1};
            if nargin == a + 1 || ...
                    ischar(vararginx{a+2}) && ...
                    strcmpi(vararginx{a+2}(1),'-')
                opts.Dcon{con} = [];
                a = a + 2;
            else
                opts.Dcon{con} = vararginx{a+2};
                a = a + 3;
            end
            con = con + 1;
            
        case '-conskipcount' % advanced
            
            % Numbers to skip when saving the contrasts
            opts.conskipcount = vararginx{a+1};
            if ischar(opts.conskipcount)
                opts.conskipcount = str2double(opts.conskipcount);
            end
            a = a + 2;
            
        case '-tonly' % advanced
            
            % Run only the t-contrasts
            opts.tonly = true;
            a = a + 1;
            
        case '-fonly' % basic
            
            % Run only the F-contrasts
            opts.fonly = true;
            a = a + 1;
            
        case '-eb' % basic
            
            % Get the exchangeability blocks file.
            opts.eb = vararginx{a+1};
            a = a + 2;
            
        case '-vg' % basic
            
            % Get the variance groups file.
            opts.vg = vararginx{a+1};
            if     ischar(opts.vg) && ...
                    any(strcmpi(opts.vg,{'single'}))
                opts.vg = 'single';
                opts.singlevg = true;
            elseif ischar(opts.vg) && ...
                    any(strcmpi(opts.vg,{'auto','automatic'}))
                opts.vg = 'auto';
                opts.singlevg = false;
            else
                opts.singlevg = false;
            end
            a = a + 2;
            
        case '-swe' % advanced
            
            % Compute one (of various possible) sandwich estimators
            opts.SwE = true;
            a = a + 1;
            
        case '-o' % basic
            
            % Output prefix for the files to be saved.
            opts.o = vararginx{a+1};
            a = a + 2;
            
        case '-n' % basic
            
            % Number of permutations
            opts.nP0 = vararginx{a+1};
            if ischar(opts.nP0)
                opts.nP0 = str2double(opts.nP0);
            end
            a = a + 2;
            
        case '-C' % basic
            
            % Threshold for cluster extent, univariate, NPC and MV
            opts.cluster.uni.do = true;
            opts.cluster.uni.thr = vararginx{a+1};
            if ischar(opts.cluster.uni.thr)
                opts.cluster.uni.thr = str2double(opts.cluster.uni.thr);
            end
            opts.cluster.npc.do = true;
            opts.cluster.npc.thr = vararginx{a+1};
            if ischar(opts.cluster.npc.thr)
                opts.cluster.npc.thr = str2double(opts.cluster.npc.thr);
            end
            opts.cluster.mv.do  = true;
            opts.cluster.mv.thr = vararginx{a+1};
            if ischar(opts.cluster.mv.thr)
                opts.cluster.mv.thr = str2double(opts.cluster.mv.thr);
            end
            a = a + 2;
            
        case '-Cuni' % advanced
            
            % Threshold for cluster statistic, univariate
            opts.cluster.uni.do = true;
            opts.cluster.uni.thr = vararginx{a+1};
            if ischar(opts.cluster.uni.thr)
                opts.cluster.uni.thr = str2double(opts.cluster.uni.thr);
            end
            a = a + 2;
            
        case '-Cnpc' % advanced
            
            % Threshold for cluster statistic, NPC
            opts.NPC = true;
            opts.cluster.npc.do = true;
            opts.cluster.npc.thr = vararginx{a+1};
            if ischar(opts.cluster.npc.thr)
                opts.cluster.npc.thr = str2double(opts.cluster.npc.thr);
            end
            a = a + 2;
            
        case '-Cmv' % advanced
            
            % Threshold for cluster statistic, MV
            opts.MV = true;
            opts.cluster.mv.do = true;
            opts.cluster.mv.thr = vararginx{a+1};
            if ischar(opts.cluster.mv.thr)
                opts.cluster.mv.thr = str2double(opts.cluster.mv.thr);
            end
            a = a + 2;
            
        case '-Cstat' % advanced
            
            % Type of cluster statistic
            opts.cluster.stat = vararginx{a+1};
            if ~ any(strcmp(opts.cluster.stat,{'extent','mass','density','tippett','pivotal'}))
                error('Cluster statistic "%s" unknown.',opts.cluster.stat);
            end
            a = a + 2;
            
        case '-T' % basic
            
            % Do TFCE for univariate, NPC and MV?
            opts.tfce.uni.do = true;
            opts.tfce.npc.do = true;
            opts.tfce.mv.do = true;
            opts.tfce.stat = 'tfce';
            a = a + 1;
            
        case '-Tstat' % not in the help
            
            % Type of cluster statistic
            opts.tfce.stat = vararginx{a+1};
            if ~ any(strcmp(opts.tfce.stat,{'tfce','density','tippett'}))
                error('TFCE statistic "%s" unknown.',opts.tfce.stat);
            end
            a = a + 2;
            
        case '-Tuni' % advanced
            
            % Do TFCE for uni?
            opts.tfce.uni.do = true;
            a = a + 1;
            
        case '-Tnpc' % advanced
            
            % Do TFCE for NPC?
            opts.NPC = true;
            opts.tfce.npc.do = true;
            a = a + 1;
            
        case '-Tmv' % advanced
            
            % Do TFCE for MV?
            opts.MV = true;
            opts.tfce.mv.do = true;
            a = a + 1;
            
        case {'-tfce1D','-tfce1d'} % basic
            
            % Shortcut for -tfce_H 2 -tfce_E 2 -tfce_C 6,
            % i.e., parameters for TFCE in 2D mode
            opts.tfce.H      = 2;
            opts.tfce.E      = 2;
            opts.tfce.conn   = 6;
            a = a + 1;
            
        case {'-tfce2D','-tfce2d'} % basic
            
            % Shortcut for -tfce_H 2 -tfce_E 1 -tfce_C 26,
            % i.e., parameters for TFCE in 2D mode
            opts.tfce.H      = 2;
            opts.tfce.E      = 1;
            opts.tfce.conn   = 26;
            a = a + 1;
            
        case {'-tfce_H','-tfce_h'} % advanced
            
            % TFCE H parameter
            opts.tfce.H = vararginx{a+1};
            if ischar(opts.tfce.H)
                opts.tfce.H = str2double(opts.tfce.H);
            end
            a = a + 2;
            
        case {'-tfce_E','-tfce_e'} % advanced
            
            % TFCE E parameter
            opts.tfce.E = vararginx{a+1};
            if ischar(opts.tfce.E)
                opts.tfce.E = str2double(opts.tfce.E);
            end
            a = a + 2;
            
        case {'-tfce_C','-tfce_c'} % advanced
            
            % TFCE connectivity
            opts.tfce.conn = vararginx{a+1};
            if ischar(opts.tfce.conn)
                opts.tfce.conn = str2double(opts.tfce.conn);
            end
            a = a + 2;
            
        case '-tfce_dh' % advanced
            
            % TFCE delta-h parameter
            opts.tfce.deltah = vararginx{a+1};
            if ischar(opts.tfce.deltah)
                if strcmpi(opts.tfce.deltah,'auto')
                    opts.tfce.deltah = 0;
                else
                    opts.tfce.deltah = str2double(opts.tfce.deltah);
                end
            end
            a = a + 2;
            
        case '-tableasvolume' % basic
            
            % Treat tables (e.g., CSV inputs) as volume, such that TFCE can
            % be calculated. This is useful for TFCE over timeseries.
            opts.tableasvolume = true;
            a = a + 1;
            
        case '-within' % basic
            
            % Define whether should permute blocks as a whole or not
            opts.within = true;
            a = a + 1;
            
        case '-whole' % basic
            
            % Define whether should permute blocks as a whole or not
            opts.whole = true;
            a = a + 1;
            
        case '-ee' % basic
            
            % Exchangeable errors (EE)?
            % If yes, this means permutations.
            opts.EE = true;
            a = a + 1;
            
        case '-ise' % basic
            
            % Independent and symmetric errors (ISE)?
            % If yes, this means sign-flippings.
            opts.ISE = true;
            a = a + 1;
            
        case '-cmcp' % advanced
            
            % Define whether Conditional Monte Carlo should be used or not,
            % that is, ignoring repeated elements in the permutation set.
            opts.cmcp = true;
            a = a + 1;
            
        case '-cmcx' % advanced
            
            % Define whether repeated rows in X should be ignored or not
            % when defining the permutations, which constitutes another
            % form of CMC
            opts.cmcx = true;
            a = a + 1;
            
        case '-twotail' % basic
            
            % Do a two-tailed test for all t-contrasts?
            opts.twotail = true;
            a = a + 1;
            
        case '-concordant' % basic
            
            % For the NPC, favour alternatives with the same sign?
            opts.concordant = true;
            a = a + 1;
            
        case '-reversemasks' % basic
            
            % Reverse masks.
            opts.reversemasks = true;
            a = a + 1;
            
        case '-corrmod' % basic
            
            % Correct over modalities.
            opts.corrmod = true;
            a = a + 1;
            
        case '-corrcon' % basic
            
            % Correct over contrasts.
            opts.corrcon = true;
            a = a + 1;
            
        case '-saveparametric' % advanced
            
            % If the user wants to have also the parametric p-values.
            opts.savepara = true;
            a = a + 1;
            
        case '-saveglm' % advanced
            
            % If the user wants, save COPE and VARCOPEs in the 1st
            % permutation.
            opts.saveglm = true;
            a = a + 1;
            
        case {'-savecdf','-save1-p'} % basic
            
            % Save 1-p values (CDF) instead of the P-values
            opts.savecdf = true;
            a = a + 1;
            
        case '-logp'  % basic
            
            % Convert the P-values or (1-P)-values to -log10 before saving.
            opts.savelogp = true;
            a = a + 1;
            
        case '-savemask' % advanced
            
            % If the user wants to have also the masks used for each.
            % modality
            opts.savemask = true;
            a = a + 1;
            
        case '-rmethod' % advanced
            
            % Which method to use for the regression/permutation?
            if narginx > a
                methlist = {           ...
                    'Draper-Stoneman', ...
                    'Still-White',     ...
                    'Freedman-Lane',   ...
                    'terBraak',        ...
                    'Kennedy',         ... % Kennedy won't be in the help
                    'Manly',           ...
                    'Huh-Jhun',        ...
                    'Dekker'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx)
                    error('Regression/Permutation method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.rmethod = methlist{methidx};
            else
                error([...
                    'The option -rmethod requires a method to be specified.\n'...
                    'Consult the documentation.%s'],'');
            end
            
        case '-npc' % basic
            
            % This is a shortcut to enable NPC with the default settings.
            opts.NPC = true;
            opts.npcmod = true;
            a = a + 1;
            
        case '-npcmethod' % basic
            
            % Do the non-parametric combination?
            if nargin == a || (nargin > a && strcmp(vararginx{a+1}(1),'-'))
                error('The option "-npcmethod" requires a combining method to be indicated.');
                
            elseif nargin > a
                
                % Which combining function to use for the combination?
                methlist = {               ...
                    'Tippett',             ...
                    'Fisher',              ...
                    'Pearson-David',       ...
                    'Stouffer',            ...
                    'Wilkinson',           ...
                    'Winer',               ...
                    'Edgington',           ...
                    'Mudholkar-George',    ...
                    'Friston',             ...
                    'Darlington-Hayes',    ...
                    'Zaykin',              ...
                    'Dudbridge-Koeleman',  ...
                    'Dudbridge-Koeleman2', ...
                    'Taylor-Tibshirani',   ...
                    'Jiang'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if ~any(methidx)
                    error('Combining method "%s" unknown.',vararginx{a+1});
                elseif any(strcmpi(vararginx{a+1},{...
                        'Wilkinson',       ...
                        'Zaykin',          ...
                        'Jiang'}))
                    if ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-')
                        plm.npcparm = 0.05;
                        a = a + 2;
                    elseif ischar(vararginx{a+2})
                        a = a + 3;
                        plm.npcparm = eval(vararginx{a+2});
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif any(strcmpi(vararginx{a+1},{...
                        'Darlington-Hayes',   ...
                        'Dudbridge-Koeleman', ...
                        'Jiang'}))
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-')
                        plm.npcparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2})
                        plm.npcparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'Friston')
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-')
                        plm.npcparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2})
                        plm.npcparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'Dudbridge-Koeleman2')
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-')
                        plm.npcparm  = 1;
                        plm.npcparm2 = 0.05;
                        a = a + 2;
                    else
                        if ischar(vararginx{a+2})
                            plm.npcparm = eval(vararginx{a+2});
                        else
                            plm.npcparm = vararginx{a+2};
                        end
                        if nargin == a + 2 || ...
                                ischar(vararginx{a+3}) && ...
                                strcmpi(vararginx{a+3}(1),'-')
                            plm.npcparm2 = 0.05;
                        elseif ischar(vararginx{a+3})
                            plm.npcparm2 = eval(vararginx{a+3});
                        else
                            plm.npcparm2 = vararginx{a+3};
                        end
                        a = a + 4;
                    end
                else
                    a = a + 2;
                end
                opts.npcmethod = methlist{methidx};
            end
            
        case '-npcmod' % basic
            
            % NPC over modalities.
            opts.NPC    = true;
            opts.npcmod = true;
            a = a + 1;
            
        case '-npccon' % basic
            
            % NPC over contrasts -- that is, all contrasts, even contrasts
            % in different designs (if more than one -d is supplied).
            opts.NPC       = true;
            opts.npccon    = true;
            opts.syncperms = true;
            a = a + 1;
            
        case '-mv' % basic
            
            % Compute classic multivariate statistics
            if nargin == a
                opts.MV = true;
                a = a + 1;
                
            elseif nargin > a && strcmp(vararginx{a+1}(1),'-')
                opts.MV = true;
                a = a + 1;
                
            elseif nargin > a
                
                % Which multivariate statistic to use?
                methlist = {            ...
                    'auto',             ...
                    'HotellingTsq',     ...
                    'Wilks',            ...
                    'Lawley',           ...
                    'Lawley-Hotelling', ...
                    'Pillai',           ...
                    'Roy',              ...
                    'Roy-ii',           ...
                    'Roy-iii',          ...
                    'CCA',              ...
                    'PLS'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if ~any(methidx)
                    error('Multivariate statistic "%s" unknown.',vararginx{a+1});
                else
                    opts.MV  = true;
                    opts.CCA = false;
                    opts.PLS = false;
                    a = a + 2;
                end
                opts.mvstat = methlist{methidx};
            end
            
        case '-fdr' % basic
            
            % Compute FDR-adjusted p-values
            opts.FDR = true;
            a = a + 1;
            
        case {'-accel','-approx'} % advanced
            
            % Choose a method to do the approximation of p-values
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-')
                methlist = {   ...
                    'negbin',  ...
                    'tail',    ...
                    'noperm',  ...
                    'gamma',   ...
                    'lowrank'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~ any(methidx)
                    error('Approximation method "%s" unknown.',vararginx{a+1});
                end
                for mm = 1:numel(methlist)
                    opts.accel.(methlist{mm}) = methidx(mm);
                end
                
                % Extra parameters
                if opts.accel.negbin
                    
                    % Number of exceedances:
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-')
                        if ischar(vararginx{a+2})
                            opts.accel.negbin = str2double(vararginx{a+2});
                        else
                            opts.accel.negbin = vararginx{a+2};
                        end
                        a = a + 3;
                    else
                        opts.accel.negbin = opts.accel.negbin_nexced;
                        a = a + 2;
                    end
                    
                elseif opts.accel.tail
                    
                    % Define whether include or not the unpermuted stat:
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-')
                        if ischar(vararginx{a+2})
                            if     any(strcmpi(vararginx{a+2},{'out','G1out','T1out','true', '1'}))
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false; % defensive, as the uncorrected will be invalid here.
                            elseif any(strcmpi(vararginx{a+2},{'in', 'G1in', 'T1in', 'false','0'}))
                                opts.accel.G1out = false;
                            end
                        else
                            if vararginx{a+2}
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false; % defensive, as the uncorrected will be invalid here.
                            else
                                opts.accel.G1out = false;
                            end
                        end
                        a = a + 3;
                    else
                        a = a + 2;
                    end
                    
                elseif opts.accel.gamma
                    
                    % Define whether include or not the unpermuted stat:
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-')
                        if ischar(vararginx{a+2})
                            if     any(strcmpi(vararginx{a+2},{'out','G1out','T1out','true', '1'}))
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false;
                            elseif any(strcmpi(vararginx{a+2},{'in', 'G1in', 'T1in', 'false','0'}))
                                opts.accel.G1out = false; % defensive, as the uncorrected will be invalid here.
                            end
                        else
                            if vararginx{a+2}
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false; % defensive, as the uncorrected will be invalid here.
                            else
                                opts.accel.G1out = false;
                            end
                        end
                        a = a + 3;
                    else
                        a = a + 2;
                    end
                    
                elseif opts.accel.lowrank
                    
                    % Fraction of voxels to be sampled (if < 1) or actual
                    % number of voxels to be sampled.
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-')
                        if ischar(vararginx{a+2})
                            opts.accel.lowrank_val = str2double(vararginx{a+2});
                        else
                            opts.accel.lowrank_val = vararginx{a+2};
                        end
                        a = a + 3;
                    else
                        a = a + 2;
                    end
                    
                else
                    a = a + 2;
                end
            else
                error([...
                    'The options "-accel" and "-approx" require a method to.\n' ...
                    'be specified. Consult the documentation.%s'],'');
            end
            
        case {'-noniiclass','-nonifticlass'} % advanced
            
            % Disable using the NIFTI class
            opts.useniiclass = false;
            a = a + 1;
            
        case '-precision' % advanced
            
            % Precision to use?
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-')
                methlist = {'single','double'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx)
                    error('Precision "%s" unknown. Use "single" or "double".',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.precision = methlist{methidx};
            else
                error([...
                    'The option "-precision" requires a method to be specified.\n'...
                    'Use "-precision double" or "-precision single".']);
            end
            
        case '-saveperms' % advanced
            
            % Save the permutations
            opts.saveperms = true;
            a = a + 1;
            
        case '-savemax' % advanced
            
            % Save the permutations
            opts.savemax = true;
            a = a + 1;
            
        case '-savemetrics' % advanced
            
            % Save a file with the number of permutations, average
            % Hamming distance, etc.
            opts.savemetrics = true;
            a = a + 1;
            
        case '-inormal' % advanced
            
            % Inverse-normal transformation?
            opts.inormal = true;
            
            % Take the parameters given to -inormal
            parms = {};
            if narginx - a >= 1
                if ~strcmp(vararginx{a+1}(1),'-')
                    parms{1} = vararginx{a+1};
                end
            end
            if narginx - a >= 2
                if ~strcmp(vararginx{a+2}(1),'-')
                    parms{2} = vararginx{a+2};
                end
            end
            a = a + 1 + numel(parms);
            
            % Then modify the variables accordingly
            methlist = {   ...
                'Blom',    ...
                'Tukey',   ...
                'Bliss',   ...
                'Waerden', ...
                'SOLAR'};
            for p = 1:numel(parms)
                methidx = strcmpi(parms{p},methlist);
                if any(methidx)
                    opts.inormal_meth = parms{p};
                elseif any(parms{p},{'quali','qualitative','discrete'})
                    opts.inormal_quanti = false;
                elseif any(parms{p},{'quanti','quantitative','continuous'})
                    opts.inormal_quanti = true;
                else
                    error('Parameter "%s" unknown for the "-inormal" option.',parms{p});
                end
            end
            
        case '-probit' % advanced
            
            % Probit transformation?
            opts.probit = true;
            a = a + 1;
            
        case '-seed' % advanced
            
            % Seed for the random number generator
            opts.seed = vararginx{a+1};
            if ischar(opts.seed) && ...
                    ~any(strcmpi(opts.seed,{'shuffle','twist','reset'}))
                opts.seed = str2double(opts.seed);
            end
            a = a + 2;
            
        case '-demean' % basic
            
            % Demean data and design. Additionally, remove
            % a global intercept, if any, from the design.
            opts.demean = true;
            a = a + 1;
            
        case '-vgdemean' % advanced
            
            % Demean data and design within VG. Additionally, remove
            % a global intercept, if any, from the design.
            opts.vgdemean = true;
            a = a + 1;
            
        case '-ev4vg' % advanced
            
            % Add to the design matrix one EV for each variance group.
            opts.ev4vg = true;
            a = a + 1;
            
        case '-removevgbysize' % advanced
            
            % Remove from the analysis observations that are the only
            % in their variance group.
            opts.removevgbysize = vararginx{a+1};
            if ischar(opts.removevgbysize)
                opts.removevgbysize = str2double(opts.removevgbysize);
            end
            a = a + 2;
            
        case '-zstat' % advanced
            
            % Convert the statistic for each test to a z-score
            opts.zstat = true;
            a = a + 1;
            
        case '-pearson' % basic
            
            % Compute the Pearson's correlation coefficient (R^2 if rank(C)>1).
            opts.pearson = true;
            a = a + 1;
            
        case '-noranktest' % advanced
            
            % Don't test the rank(Y) before doing MANOVA/MANCOVA.
            opts.noranktest = true;
            a = a + 1;
            
        case '-transposedata' % advanced
            
            % Transpose the data if it's 2D?
            opts.transposedata = true;
            a = a + 1;
            
        case '-inputmv' % advanced
            
            % Treat the (sole) input as multivariate, that is,
            % each column is a variable in a multivariate model,
            % as opposed to independent univariate tests.
            % Useful with non-imaging data.
            opts.inputmv = true;
            a = a + 1;
            
        case '-verbosefilenames' % advanced
            
            % Don't simplify filenames when otherwise it would be possible
            opts.verbosefilenames = true;
            a = a + 1;
            
        case '-syncperms' % advanced
            
            % Synchronize permutations regardless of other options?
            opts.syncperms = true;
            a = a + 1;
            
        case {'-designperinput','-erie'} % advanced
            
            % Use one design matrix for each input dataset?
            % This enables
            % syncP regardless.
            opts.designperinput = true;
            opts.syncperms      = true;
            a = a + 1;
            
        case '-subjidx' % advanced
            
            % Indices of the subjects to keep in the design
            opts.subjidx = vararginx{a+1};
            a = a + 2;
            
        case '-quiet' % basic
            
            % Don't print lines indicating progress
            opts.showprogress = false;
            a = a + 1;
            
        case '-nounivariate' % advanced
            
            % Save or not univariate tests? Useful with MV/NPC/CCA
            opts.saveunivariate = false;
            a = a + 1;
            
        case '-nouncorrected'  % advanced
            
            % Save or not uncorrected p-vals
            opts.saveuncorrected = false;
            a = a + 1;
            
        case '-saveuncorrected'  % advanced
            
            % Save or not uncorrected p-vals
            opts.saveuncorrected = true;
            a = a + 1;
            
        case '-savedof' % advanced
            
            % Save a file with the degrees of freedom?
            opts.savedof = true;
            a = a + 1;
            
        case '-pmethodp' % advanced
            
            % Which method to use to partition the model when defining
            % the permutations?
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-')
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Ridgway',  ...
                    'none'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx)
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethodp = methlist{methidx};
            else
                error([...
                    'The option "-pmethodp" requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
        case '-pmethodr' % advanced
            
            % Which method to use to partition the model when defining
            % doing the actual regression?
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-')
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Ridgway',  ...
                    'none'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx)
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethodr = methlist{methidx};
            else
                error([...
                    'The option "-pmethodr" requires a method to be specified.\n' ...
                    'Consult the documentation.']);
            end
            
        case '-forceintersectionmask' % currently not in the help
            
            % Force the generation of an intersection mask, even if there
            % is no MV or NPC (useful for mediation analysis).
            opts.forcemaskinter = true;
            a = a + 1;
            
        otherwise
            error('Unknown option: "%s"',vararginx{a});
    end
end


% Check for the existence of other programs for input/output
palm_checkprogs;

% Run basic checks for input data
if ~opts.cca.do && Ni == 0
    error('Missing input data (missing "-i").');
elseif opts.cca.do && Nx == 0
    error('Missing input data (missing "-x" for cca).');
elseif opts.cca.do && Ny == 0
    error('Missing input data (missing "-y" for cca).');
end

% For CCA, the -semipartial option should not work if any -w is supplied
if isfield(opts.cca,'semipartial') && ~isempty(opts.w) 
    error('Option -semipartial should not be used if supplying w.');
end

% If doing CCA, determine if it's voxelwise or spatial 
if opts.cca.do && Nx > 1
    opts.cca.voxelwise = true;
else
    opts.cca.voxelwise = false;
end

% Make sure the NPC options make sense
% - if -npc only, it's -npcmod
% - if -npcmod only, it's -npcmod
% - if -npccon only, it's -npccon
% - if -npcmod and -npccon, it's both
if opts.NPC && ~ opts.npccon
    opts.npcmod = true;
end

% A quick check for the case of 1 EV per column of Y.
if opts.evperdat
    if any([...
            opts.ev4vg
            opts.pearson])
        error([...
            'The option "-evperdat" is incompatible with the options listed below:\n' ...
            '"-ev4vg"\n' ...
            '"-pearson"%s'],'');
    end
    if strcmpi(opts.rmethod,'terBraak')
        error('The option "-evperdat" cannot be used with the ter Braak method (not implemented)');
    end
end


% No spatial statistics if no univariate results will be saved anyway
if ~ opts.saveunivariate
    opts.cluster.uni.do = false;
    opts.tfce.uni.do    = false;
end

% This simplifies some tests later
opts.spatial.do  = false;
opts.spatial.uni = false;
opts.spatial.npc = false;
opts.spatial.mv  = false;
if any([ ...
        opts.cluster.uni.do  ...
        opts.cluster.npc.do  ...
        opts.cluster.mv.do   ...
        opts.tfce.uni.do     ...
        opts.tfce.npc.do     ...
        opts.tfce.mv.do]')
    opts.spatial.do = true;
    if any([ ...
            opts.cluster.uni.do  ...
            opts.tfce.uni.do]')
        opts.spatial.uni = true;
    end
    if any([ ...
            opts.cluster.npc.do  ...
            opts.tfce.npc.do]')
        opts.spatial.npc = true;
    end
    if any([ ...
            opts.cluster.mv.do   ...
            opts.tfce.mv.do]')
        opts.spatial.mv = true;
    end
end

% Adjust the delta-h.
if any([ ...
        opts.tfce.uni.do
        opts.tfce.npc.do
        opts.tfce.mv.do ]) && ...
        strcmpi(opts.tfce.deltah,'auto')
    opts.tfce.deltah = 0;
end

% Sanity checks for the acceleration modes.
if sum(logical([ ...
        opts.accel.negbin, ...
        opts.accel.tail,   ...
        opts.accel.noperm, ...
        opts.accel.gamma,  ...
        opts.accel.lowrank])) > 1
    error('Only one approximation method can be used for a given run.');
end
if opts.accel.negbin
    if opts.accel.negbin < 2
        error('The parameter r given to "-accel negbin <r>" must be >= 2.')
    end
    if opts.nP0 < 3
        error('The option "-accel negbin <r>" needs at least 3 permutations.')
    end
    if ~ opts.saveuncorrected
        error('The option "-nouncorrected" cannot be used with "-accel negbin".');
    end
    if (opts.corrmod || opts.corrcon) && ~ opts.FDR
        error('The option "-accel negbin" cannot be used with FWER-correction, only FDR.');
    end
    if opts.NPC
        error('The option "-accel negbin" cannot be used with NPC.');
    end
    if opts.spatial.do
        error('The option "-accel negbin" cannot be used with spatial statistics.');
    end
    if opts.saveperms
        error('The option "-saveperms" cannot be used together with "-accel negbin".');
    end
end
if opts.accel.noperm
    if ~ opts.saveuncorrected
        error('The option "-nouncorrected" cannot be used with "-accel noperm".');
    end
    if (opts.corrmod || opts.corrcon) && ~ opts.FDR
        error('The option "-accel noperm" cannot be used with FWER-correction, only FDR.');
    end
    if opts.NPC
        error('The option "-accel noperm" cannot be used with NPC.');
    end
    if opts.spatial.do
        error('The option "-accel noperm" cannot be used with spatial statistics.');
    end
    if opts.MV
        if ~ any(strcmpi(opts.mvstat,{'auto','pillai'}))
            warning([...
                'With multivariate tests, the option "-accel noperm" can\n' ...
                '         only be used with the Pillai'' trace statistic.\n' ...
                '         Changing automatically to Pillai.%s'],'');
        end
        opts.mvstat = 'pillai';
    end
    if Ni > 1 && ~ opts.MV
        error([...
            'The option "-accel noperm" needs to be used with a single modality\n'...
            '       modality or with "-mv".%s'],'');
    end
end
if opts.accel.lowrank
    if opts.pearson
        error('The option "-accel lowrank" cannot be used with "-pearson".');
    end
    if opts.MV
        error('The option "-accel lowrank" cannot be used with MV.');
    end
    if opts.NPC
        error('The option "-accel lowrank" cannot be used with NPC.');
    end
    if opts.spatial.do
        error('The option "-accel lowrank" cannot be used with spatial statistics.');
    end
    if opts.evperdat
        error('The option "-accel lowrank" cannot be used with "-evperdat".');
    end
    if opts.missingdata
        error('The option "-accel lowrank" cannot be used with missing data.');
    end
    if opts.saveglm
        error('The option "-accel lowrank" cannot be used with "-saveglm".');
    end
end

% Some options can't be used with missing data
if Nimiss || Ndmiss
    opts.missingdata = true;

    if opts.MV
        error('The option "-mv" cannot be used with missing data.');
    end
    if ~ opts.zstat && ~ opts.mcar
        warning([...
            'With missing data MAR/MNAR, the option "-zstat" is mandatory.\n' ...
            '         Adding it automatically.%s'],'');
        opts.zstat = true;
    end
    if ~ opts.cmcx && ~ opts.mcar
        warning([...
            'With missing data, the option "-cmcx" is mandatory.\n' ...
            '         Adding it automatically.%s'],'');
        opts.cmcx = true;
    end
    if opts.demean && ~ opts.mcar
        warning([...
            'With missing data, the option "-demean" must not be used.\n' ...
            '         Removing it automatically.%s'],'');
        opts.demean = false;
    end
    if ~ strcmpi(opts.pmethodp,'guttman') || ~ strcmpi(opts.pmethodr,'guttman')
        warning([...
            'With missing data, the partitioning must use the "Guttman".\n' ...
            '         method. Adding automatically the options\n' ...
            '         "-pmethodp Guttman" and "-pmethodr Guttman".%s'],'');
        opts.pmethodp = 'guttman';
        opts.pmethodr = 'guttman'; 
    end
    if opts.ev4vg || opts.removevgbysize > 0
        error('Missing data cannot be used with "-ev4vg" or "-removevgbysize".')
    end
end

% Some NPC methods don't have an analytical form for the parametric p-value
if opts.NPC && any(strcmpi(opts.npcmethod,{'Darlington-Hayes','Jiang'}))
    plm.nonpcppara = true;
    if opts.savepara
        warning([...
            'No parametric combination p-value will be saved for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
    end
    if opts.spatial.npc
        warning([ ...
            'No NPC spatial statistic will be produced for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
        opts.cluster.npc.do = false;
        opts.tfce.npc.do    = false;
        opts.spatial.npc    = false;
    end
else
    plm.nonpcppara = false;
end

% Likewise, some MV methods don't have an analytical form for the parametric p-value
if opts.MV && strcmpi(opts.mvstat,'Roy_iii')
    plm.nomvppara = true;
    if opts.savepara
        warning('No parametric p-value will be saved for the Roy_iii method.%s','');
    end
    if opts.spatial.mv
        warning([ ...
            'No multivariate cluster-level or TFCE statistic will be produced\n', ...
            '         for the Roy_iii statistic%s'],'');
        opts.cluster.mv.do = false;
        opts.tfce.mv.do    = false;
        opts.spatial.mv    = false;
    end
elseif (opts.CCA || opts.PLS) && opts.savepara
    plm.nomvppara = true;
    warning([...
        'No parametric p-value will be saved for CCA or PLS.\n', ...
        '         Use Wilks'' lambda instead.%s'],'');
else
    plm.nomvppara = false;
end

% Some more warnings and sanity checks
if opts.savecdf && opts.savelogp
    error('Should not use "-save1-p" together with "-logp"');
end
if ~ opts.inputmv && opts.designperinput && Ni ~= Nd
    error([
        'To use the option "-designperinput", the number of design files must\n' ...
        'match the number of inputs.\n%s'],'');
end
if (Nt || Nf) && Ncon
    error('Cannot mix options "-t" or "-f" with "-con".');
end
if Nt || Nf
    if Nt > Nd
        error('More t-contrast files (%d) than valid design files (%d) were supplied.',Nt,Nd);
    end
    if Nt ~= 1 && Nt ~= Nd
        error(['The number of supplied t-contrast files (option "-t") must be 1 or\n',...
            'the same as the number of design files (option "-d") (%d).'],Nd);
    end
    if Nf > Nt
        error('More F-contrast files (%d) than t-contrast files (%d) were supplied.',Nf,Nt);
    end
elseif Ncon
    if Ncon > Nd
        error('More contrast files (%d) than design files (%d) were supplied.',Nt,Nd);
    end
    if Ncon ~= 1 && Ncon ~= Nd
        error(['The number of supplied contrast files (option "-con") must be 1 or\n',...
            'the same as the number of design files (option "-d") (%d).'],Nd);
    end
end
if opts.pearson && (opts.NPC || opts.MV)
    error([ ...
        'It isn''t possible to compute the Pearson''s r or R^2 together with NPC or\n', ...
        '         multivariate methods.%s'],'');
end
if opts.pearson && ~ opts.demean
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the data and the design\n' ...
        '         must be mean centered. Adding option "-demean".%s'],'');
    opts.demean = true;
end
if opts.pearson && opts.savepara
    error([ ...
        'The option "-saveparametric" cannot be used with "-pearson".\n' ...
        '         For a parametric p-value, drop "-pearson" and use the\n' ...
        '         respective results for the t or F-statistics.%s'],'');
end
if opts.pearson && ~ any(strcmpi(opts.pmethodr,{'beckmann','ridgway'}))
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option "-pmethodr Beckmann".%s'],'');
    opts.pmethodr = 'beckmann';
end
if (opts.cca.do || opts.pls.do) && ~ opts.demean
    warning([ ...
        'To perform CCA or PLS, the data and the design\n' ...
        '         must be mean centered. Adding option "-demean".%s'],'');
    opts.demean = true;
end
if (opts.CCA || opts.PLS) && ~ any(strcmpi(opts.pmethodr,{'beckmann','ridgway'}))
    warning([ ...
        'To perform CCA or PLS, the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option "-pmethodr Beckmann".%s'],'');
    opts.pmethodr = 'beckmann';
end
if opts.demean && opts.vgdemean && ~opts.pearson && ~opts.CCA && ~opts.PLS
    warning([...
        'Cannot use the option "-demean" together with "-vgdemean"\n'...
        '         Ignoring the option "-vgdemean".%s'],'');
    opts.vgdemean = false;
end
if opts.ev4vg && opts.vgdemean
    error('Cannot use the option "-ev4vg" together with "-vgdemean"');
end
if opts.MV && ~ opts.singlevg
    error('Currently MV is only allowed with a single VG, that is, assuming homoscedasticity.');
end
if opts.designperinput && opts.MV
    error('It''s not possible to use the option "-designperinput" together with the option "-mv".');
end
if ~opts.saveunivariate && ~opts.MV && ~opts.NPC && ~opts.CCA && ~opts.PLS
    error('The option "-nounivariate" can only be used with "-mv" or "-npc".');
end
if opts.MV && (opts.cca.do || opts.pls.do)
    error('Cannot do classical MANOVA at the same time as CCA or PLS.');
end
if Ni > 1 && opts.inputmv
    error('Option "-inputmv" cannot be used with more than one "-i".')
end
if opts.concordant && ~ opts.NPC
    error('The option "-concordant" is for use with NPC only.');
end
if opts.concordant && opts.twotail
    error(['Cannot use "-concordant" together with "-twotail" (inadmissible).\n'...
        'Use either of these, but not both together.%s'],'');
end
if opts.tonly && opts.fonly
    error('Cannot use "-tonly" together with "-fonly".');
end
if opts.probit && opts.inormal
    error('Cannot use "-probit" together with "-inormal".');
end
if opts.FDR && ~ opts.saveuncorrected
    error(['Option "-fdr" cannot be used together with "-nouncorrected".\n'...
        'Use either of these, but not both together. If you are using tail or\n'...
        'gamma approximations, consider keeping "-nouncorrected", and leave\n'...
        '"-fdr" for a separate call in which tail or gamma are not used.%s'],'');
end

% Initialize the random number generator (if nP = 0, no need for that)
if opts.nP0
    if palm_isoctave
        if any(strcmpi(opts.seed,{'reset','shuffle','twist'}))
            opts.seed = 'reset';
        end
        rand('state',opts.seed); %#ok
    else
        if any(strcmpi(opts.seed,{'reset','shuffle','twist'}))
            opts.seed = 'shuffle';
        end
        rng('default');
        rng(opts.seed);
    end
end

% Read and organise the surfaces. If no surfaces have been loaded, but the
% user wants cluster extent, cluster mass, or TFCE, an error will be
% printed later down in the code.
% At this stage the average area from native geometry is also loaded if it
% was provided. If a weight was given, such as 1, use this weight, so that
% all faces are treated as if having the same size. If nothing was
% supplied, use the actual area of the surface.
if opts.spatial.do && Ns > 0
    plm.srf     = cell(Ns,1);
    plm.srfarea = cell(Ns,1);
    plm.srfadj  = cell(Ns,1);
    for s = 1:Ns
        
        % Load surface
        fprintf('Loading surface %d/%d: %s\n',s,Ns,opts.s{s});
        plm.srf{s} = palm_miscread(opts.s{s},[],[],[],true);
        
        % Load areas
        if isempty(opts.sa{s})
            % No area means area of the actual surface file
            plm.srfarea{s}.data = [];
        elseif exist(opts.sa{s},'file')
            % A file with the average areas from native geometry
            plm.srfarea{s} = palm_miscread(opts.sa{s},opts.useniiclass,opts.o,opts.precision,false);
        elseif ~ isnan(str2double(opts.sa{s}))
            % A weight (such as 1)
            plm.srfarea{s}.data = str2double(opts.sa{s});
        else
            % Otherwise gives a helpful error message:
            error('Unknown option given to "-s" or file doesn''t exist:\n%s',opts.sa{s});
        end
    end
end

% There should be no more masks than modalities, and the number of
% masks needs to be either 1 or the same number of modalities.
if Nm > Ni
    error([...
        'There are more masks supplied with "-m" (%d masks) than\n'...
        'modalities supplied with "-i" (%d modalities)'],Nm,Ni);
elseif Nm > 1 && Nm ~= Ni
    error([...
        'The number of masks supplied with "-m" (%d masks) is larger than 1,\n'...
        'but still not the same as the number of modalities supplied with\n'...
        'the option "-i" (%d modalities).'],Nm,Ni);
end

% Read and organise the masks. If there are no masks specified, one for
% each modality will be created after each modality is loaded.
plm.masks = cell(Ni,1);
for m = 1:Nm
    plm.masks{m} = palm_miscread(opts.m{m},opts.useniiclass,opts.o,opts.precision,false);
    if strcmp(plm.masks{m}.readwith,'nifticlass')
        plm.masks{m}.data = double(plm.masks{m}.data);
    end
    if opts.reversemasks
        plm.masks{m}.data(isnan(plm.masks{m}.data)) = 1;
        plm.masks{m}.data(isinf(plm.masks{m}.data)) = 1;
        plm.masks{m}.data = ~ logical(plm.masks{m}.data);
    else
        plm.masks{m}.data(isnan(plm.masks{m}.data)) = 0;
        plm.masks{m}.data(isinf(plm.masks{m}.data)) = 0;
        plm.masks{m}.data = logical(plm.masks{m}.data);
    end
end
if Nm == 1
    for i = 2:Ni
        plm.masks{i} = plm.masks{1};
    end
end

% Indices of the subjects that will be kept
if ~ isempty(opts.subjidx)
    plm.subjidx = palm_miscread(opts.subjidx);
    plm.subjidx = round(plm.subjidx.data);
    u = unique(plm.subjidx);
    if numel(u) == 1
        error('Subject indices are all identical in the -subjidx file.');
    end
end
