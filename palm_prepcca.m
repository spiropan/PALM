function [opts,plm] = palm_prepcca(opts,plm)
% This will prepare additional inputs required for PALM CCA.
%
% _____________________________________
% Anderson M. Winkler and Spiro P. Pantazatos
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

% Determine sizes of input data
Ny = size(opts.y,1);
Nx = size(opts.x,1);
Nz = size(opts.z,1);
Nw = size(opts.w,1);
Ns = size(opts.s,1);
Nm = size(opts.m,1);
Nd = 1;

% Read and organise the data.
% Initialize variables
plm.Yset     = cell(Ny,1);  % Regressands (Y)
plm.Yisvol   = false(Ny,1); % Is Y a volume image?
plm.Yissrf   = false(Ny,1); % Is Y a surface-based image (DPX)?
plm.Yisvtx   = false(Ns,1); % Is vertexwise?
plm.Yisfac   = false(Ns,1); % is facewise? (this is currently dichotomous with Yisvtx, but later there may be edges/connectivity too)
plm.Yarea    = cell(Ns,1);  % To store area per face or per vertex (used for cluster-level & TFCE).
plm.Ykindstr = cell(Ny,1);  % string to save the files later

% Read and organize Y
for i = 1:Ny
    % Read input file
    fprintf('Reading input %d/%d: %s\n',i,Ny,opts.y{i});
    if Nm == 0
        maskstruct = [];
    else
        maskstruct = plm.masks{i};
    end
    [plm.Yset{i},plm.Ymasks{i},plm.Yisvol(i),plm.Yissrf(i),plm.Ykindstr{i}] = ...
        palm_ready(opts.y{i},maskstruct,opts,true);
    
    % Select subjects
    if ~ isempty(plm.subjidx)
        plm.Yset{i} = plm.Yset{i}(plm.subjidx,:);
    end
    
    % For the first input data, keep the size to
    % compare with the others, then check the size
    if i == 1
        plm.N = size(plm.Yset{i},1);
    end
    if size(plm.Yset{i},1) ~= plm.N
        error([
            'At least two of the Y input data files do not have\n' ...
            'compatible sizes:\n' ...
            '- File %d (%s) has %d observations\n' ...
            '- File %d (%s) has %d observations'], ...
            1,opts.y{1},plm.N, ...
            i,opts.y{i},size(plm.Yset{i},1));
    end
end

% Read and organize X
for i = 1:Nx
    % Read input file
    fprintf('Reading input %d/%d: %s\n',i,Nx,opts.x{i});
    if Nm == 0
        maskstruct = [];
    else
        maskstruct = plm.masks{i};
    end
    [plm.Xset{i},plm.Xmasks{i},plm.Xisvol(i),plm.Xissrf(i),plm.Xkindstr{i}] = ...
        palm_ready(opts.x{i},maskstruct,opts,true);
    
    % Select subjects
    if ~ isempty(plm.subjidx)
        plm.Xset{i} = plm.Xset{i}(plm.subjidx,:);
    end
    
    % For the first input data, keep the size to
    % compare with the others, then check the size
    if i == 1
        plm.N = size(plm.Xset{i},1);
    end
    if size(plm.Xset{i},1) ~= plm.N
        error([
            'At least two of the X input data files do not have\n' ...
            'compatible sizes:\n' ...
            '- File %d (%s) has %d observations\n' ...
            '- File %d (%s) has %d observations'], ...
            1,opts.x{1},plm.N, ...
            i,opts.x{i},size(plm.Xset{i},1));
    end
end

plm.nmasksY = numel(plm.Ymasks);
plm.nmasksX = numel(plm.Xmasks);
plm.nY      = numel(plm.Yset); 
plm.nX      = numel(plm.Xset);

% Read the zevperdat input data (modeled after EV per datum in palm_prepglm)
if opts.zevperdat
    plm.ZEVset  = cell(plm.nZEVdat,1);
    
    % If there's one mask, use the same one
    % for all the input files. Otherwise, create them.
    if Nm == 1
        for ev = 1:plm.nZEVdat
            plm.masksZEV{ev} = plm.Ymasks{1};
        end
    else
        plm.masksZEV = cell(plm.nZEVdat,1);
    end
    
    % Read input file & select subjects
    for ev = 1:plm.nZEVdat
        fprintf('Reading zevperdat inputs %d/%d: %s\n',ev,plm.nZEVdat,opts.zevdatfile{ev});
        [plm.ZEVset{ev},plm.masksZEV{ev}] = palm_ready(opts.zevdatfile{ev},plm.masksZEV{ev},opts,opts.demean);
        if ~ isempty(plm.subjidx)
            plm.ZEVset{ev} = plm.ZEVset{ev}(plm.subjidx,:);
        end
    end
    plm.nmasksZEV = numel(plm.masksZEV);
    
    % Make the intersection of the ZEVs 
    newmask = true(size(plm.masksZEV{1}.data));
    for d = 1:plm.nmasksZEV,
        newmask = newmask & plm.masksZEV{d}.data;
        plm.ZEVset{d} = plm.ZEVset{d}(:,newmask(plm.masksZEV{d}.data(:)));
        plm.masksZEV{d}.data = newmask;
    end
    
    % Then make the intersections with the respective masks of the input files
    if (plm.nY == 1 && plm.nX == 1) % Spatial CCA
        ZEV = 1:plm.nZEVdat;
        newmask = zeros(numel(plm.Ymasks{1}.data),numel(ZEV)+1);
        c = 1;
        for ev = ZEV'
            newmask(:,c) = plm.masksZEV{ev}.data(:);
            c = c + 1;
        end
        newmask(:,c) = plm.Ymasks{1}.data(:);
        newmask = reshape(all(newmask,2),size(plm.Ymasks{1}.data));
        for ev = ZEV'
            plm.ZEVset{ev} = plm.ZEVset{ev}(:,newmask(plm.masksZEV{ev}.data(:)));
            plm.masksZEV{ev}.data = newmask;
        end
        plm.Yset{1} = plm.Yset{1}(:,newmask(plm.Ymasks{1}.data(:)));
        plm.Ymasks{1}.data    = newmask;
    else
        sizy = zeros(plm.nmasksY,1);
        for y = 1:plm.nmasksY
            sizy(y) = numel(plm.Ymasks{y}.data);
        end
        sizev = zeros(plm.nmasksZEV,1);
        for ev = 1:plm.nmasksZEV
            sizev(ev) = numel(plm.masksZEV{ev}.data);
        end
        if  numel(unique(sizy)) > 1 || ...
                numel(unique(sizev)) > 1 || ...
                sizy(1) ~= sizev(1)
            error([...
                'For multiple "-i" and/or "-evperdat", with "-npccon",\n',...
                'and without the option "-designperinput", the inputs, EVs \n'...
                'and masks need to be all of the same sizes.%s'],'');
        end
        newmask = true(size(plm.masksZEV{1}.data));
        for y = 1:plm.nmasksY
            newmask = newmask & plm.Ymasks{y}.data;
        end
        for ev = 1:plm.nmasksZEV
            newmask = newmask & plm.masksZEV{ev}.data;
        end
        for y = 1:plm.nmasksY
            plm.Yset{y} = plm.Yset{y}(:,newmask(plm.Ymasks{y}.data(:)));
            plm.Ymasks{y}.data = newmask;
        end
        for x = 1:plm.nmasksX
            plm.Xset{x} = plm.Xset{x}(:,newmask(plm.Xmasks{x}.data(:)));
            plm.Xmasks{x}.data = newmask;
        end
        for ev = 1:plm.nmasksZEV
            plm.ZEVset{ev} = plm.ZEVset{ev}(:,newmask(plm.masksZEV{ev}.data(:)));
            plm.masksZEV{ev}.data = newmask;
        end
    end
end
clear('newmask');

% Read the wevperdat input data (modeled after EV per datum in palm_prepglm)
if opts.wevperdat
    plm.WEVset  = cell(plm.nWEVdat,1);
    
    % If there's one mask, use the same one
    % for all the input files. Otherwise, create them.
    if Nm == 1
        for ev = 1:plm.nWEVdat
            plm.masksWEV{ev} = plm.Xmasks{1};
        end
    else
        plm.masksWEV = cell(plm.nWEVdat,1);
    end
    
    % Read input file & select subjects
    for ev = 1:plm.nWEVdat
        fprintf('Reading wevperdat inputs %d/%d: %s\n',ev,plm.nWEVdat,opts.wevdatfile{ev});
        [plm.WEVset{ev},plm.masksWEV{ev}] = palm_ready(opts.wevdatfile{ev},plm.masksWEV{ev},opts,opts.demean);
        if ~ isempty(plm.subjidx)
            plm.WEVset{ev} = plm.WEVset{ev}(plm.subjidx,:);
        end
    end
    plm.nmasksWEV = numel(plm.masksWEV);
    
    % Make the intersection of the ZEVs 
    newmask = true(size(plm.masksWEV{1}.data));
    for d = 1:plm.nmasksWEV,
        newmask = newmask & plm.masksWEV{d}.data;
        plm.WEVset{d} = plm.WEVset{d}(:,newmask(plm.masksWEV{d}.data(:)));
        plm.masksWEV{d}.data = newmask;
    end
    
    % Then make the intersections with the respective masks of the input files
    if (plm.nY == 1 && plm.nX == 1) % Spatial CCA
        WEV = 1:plm.nWEVdat;
        newmask = zeros(numel(plm.Ymasks{1}.data),numel(WEV)+1);
        c = 1;
        for ev = WEV'
            newmask(:,c) = plm.masksWEV{ev}.data(:);
            c = c + 1;
        end
        newmask(:,c) = plm.Xmasks{1}.data(:);
        newmask = reshape(all(newmask,2),size(plm.Xmasks{1}.data));
        for ev = WEV'
            plm.WEVset{ev} = plm.WEVset{ev}(:,newmask(plm.masksWEV{ev}.data(:)));
            plm.masksWEV{ev}.data = newmask;
        end
        plm.Yset{1} = plm.Yset{1}(:,newmask(plm.Ymasks{1}.data(:)));
        plm.Ymasks{1}.data    = newmask;
    else
        sizy = zeros(plm.nmasksY,1);
        for y = 1:plm.nmasksY
            sizy(y) = numel(plm.Ymasks{y}.data);
        end
        sizev = zeros(plm.nmasksWEV,1);
        for ev = 1:plm.nmasksWEV
            sizev(ev) = numel(plm.masksWEV{ev}.data);
        end
        if  numel(unique(sizy)) > 1 || ...
                numel(unique(sizev)) > 1 || ...
                sizy(1) ~= sizev(1)
            error([...
                'For multiple "-i" and/or "-evperdat", with "-npccon",\n',...
                'and without the option "-designperinput", the inputs, EVs \n'...
                'and masks need to be all of the same sizes.%s'],'');
        end
        newmask = true(size(plm.masksWEV{1}.data));
        for y = 1:plm.nmasksY
            newmask = newmask & plm.Ymasks{y}.data;
        end
        for ev = 1:plm.nmasksWEV
            newmask = newmask & plm.masksWEV{ev}.data;
        end
        for y = 1:plm.nmasksY
            plm.Yset{y} = plm.Yset{y}(:,newmask(plm.Ymasks{y}.data(:)));
            plm.Ymasks{y}.data = newmask;
        end
        for x = 1:plm.nmasksX
            plm.Xset{x} = plm.Xset{x}(:,newmask(plm.Xmasks{x}.data(:)));
            plm.Xmasks{x}.data = newmask;
        end
        for ev = 1:plm.nmasksWEV
            plm.WEVset{ev} = plm.WEVset{ev}(:,newmask(plm.masksWEV{ev}.data(:)));
            plm.masksWEV{ev}.data = newmask;
        end
    end
end
clear('newmask');

% Create an intersection mask if NPC or MV is to be done, and further apply
% to the data that was previously masked above, as needed.
if opts.cca.do || opts.forcemaskinter
    if plm.nmasksY > 1
        
        % If there is one mask per modality, make an instersection mask.
        maskinter = true(size(plm.Ymasks{1}.data));
        for m = 1:plm.nmasksY
            maskinter = maskinter & plm.Ymasks{m}.data;
        end
        if opts.zevperdat
            for ev = 1:plm.nZEVdat
                maskinter = maskinter & plm.masksZEV{ev}.data;
            end
        end
        
        % Note that this line below uses Ytmp, which is from the previous loop.
        % This can be used here because with NPC all data has the same size.
        plm.maskinter = palm_maskstruct(maskinter(:)',plm.Ymasks{1}.readwith,plm.Ymasks{1}.extra);
        
        % Apply it to further subselect data points
        for y = 1:plm.nY
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.Ymasks{y}.data));
        end
        if opts.zevperdat
            for ev = 1:plm.nZEVdat
                plm.ZEVset{ev} = plm.ZEVset{ev}(:,plm.maskinter.data(plm.masksZEV{ev}.data));
            end
        end
    else
        % If only one mask was given.
        plm.maskinter = plm.Ymasks{1};
        for y = 1:plm.nY
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.Ymasks{1}.data));
        end
    end
end

% Sizes for later
if opts.zevperdat
    plm.ZEVsiz = zeros(plm.nZEVdat,1);
    for ev = 1:plm.nZEVdat
        plm.ZEVsiz(ev) = size(plm.ZEVset{ev},2);
    end
end

% Sizes for later
if opts.wevperdat
    plm.WEVsiz = zeros(plm.nWEVdat,1);
    for ev = 1:plm.nWEVdat
        plm.WEVsiz(ev) = size(plm.WEVset{ev},2);
    end
end


% Make sure that all data have the same size
if opts.cca.do || opts.forcemaskinter
    % The plm.Ysiz is redefined below.
    for y = 1:plm.nY
        plm.Ysiz(y) = size(plm.Yset{y},2);
    end
    [usiz,uidx] = unique(plm.Ysiz);
    if numel(usiz) > 2 || (numel(usiz) == 2 && min(usiz) ~= 1)
        error('The sizes of some of the imaging modalities don''t match');
    elseif numel(usiz) == 2 && usiz(1) == 1
        for y = 1:plm.nY
            if plm.Ysiz(y) == 1
                fprintf('Expanding modality #%d to match the size of the others.\n',y);
                plm.Yset{y} = repmat(plm.Yset{y},[1 usiz(2)]);
                if plm.nmasks > 1
                    if numel(plm.masks{y}.data) == 1
                        plm.masks{y}.data = plm.masks{uidx(2)}.data;
                    else
                        error('Modality expansion is only allowed for single input variables.')
                    end
                end
            end
        end
    end
    clear('usiz');
end

% Make sure none of the data inputs is empty
for y = 1:plm.nY
    if any(size(plm.Yset{y}) == 0)
        error('Modality %d has no data.\n',y);
    end
end

for x = 1:plm.nX
    if any(size(plm.Xset{x}) == 0)
        error('Modality %d has no data.\n',x);
    end
end

% Variables with the cumulative sizes of all data inputs will be handy later
plm.Ysiz = zeros(plm.nY,1);
for y = 1:plm.nY
    plm.Ysiz(y) = size(plm.Yset{y},2);
end
plm.Ycumsiz = vertcat(0,cumsum(plm.Ysiz));

plm.Xsiz = zeros(plm.nX,1);
for y = 1:plm.nX
    plm.Xsiz(y) = size(plm.Xset{y},2);
end
plm.Xcumsiz = vertcat(0,cumsum(plm.Xsiz));


% Take this opportunity to save the masks if the user requested.
if opts.savemask
    for y = 1:plm.nmasks
        M = plm.masks{y};
        if plm.nY == 1 || plm.nmasks == 1
            M.filename = sprintf('%smask',opts.o);
        else
            M.filename = sprintf('%smask_m%d',opts.o,y);
        end
        M.data = double(M.data);
        palm_miscwrite(M,true);
    end
    if plm.nY > 1 && (opts.npcmod || opts.MV || opts.forcemaskinter)
        M          = plm.maskinter;
        M.filename = sprintf('%sintersection_mask',opts.o);
        M.data     = double(M.data);
        palm_miscwrite(M,true);
    end
end

% If CCA was selected, make sure that Y is full rank.
if opts.cca.do && ~ opts.noranktest
    fprintf('Testing rank of the data for the MV tests. To skip, use -noranktest.\n')
    Y = cat(3,plm.Yset{:});
    Y = permute(Y,[1 3 2]);
    failed = false(1,size(Y,3));
    for v = 1:size(Y,3)
        if rank(Y(:,:,v)) ~= plm.nY
            failed(v) = true;
        end
    end
    if any(failed)
        fname = sprintf('%smv_illconditioned',opts.o);
        palm_quicksave(double(failed),0,opts,plm,[],[],[],fname);
        error([
            'One or more datapoints have ill-conditioned data. It is\n' ...
            'not possible to run multivariate analyses as MANOVA/MANCOVA.\n' ...
            'Please, see these datapoints marked as 1 in the file:\n' ...
            '%s.*\n'],fname); %#ok
    end
end; clear Y;

% Applies an inverse-normal transformation to the modalities if the user requested
if opts.inormal
    for y = 1:plm.nY
        plm.Yset{y} = palm_inormal( ...
            plm.Yset{y},            ...
            opts.inormal_meth,      ...
            opts.inormal_quanti);
    end
end

% Applies a probit transformation to the modalities if the user requested
if opts.probit
    for y = 1:plm.nY
        if min(plm.Yset{y}(:)) < 0 || max(plm.Yset{y}(:)) > 1
            error([
                'Probit transformation can only be used with data in the interval [0 1].\n' ...
                'This fails for at least modality #%d.'],y);
        end
        plm.Yset{y} = erfinv(2*(plm.Yset{y}*0.999999999999999 + .5e-15)-1)*sqrt(2);
    end
end

% Make the adjustments for the EE and ISE options.
% - if the user gives nothing, its EE by default.
% - if the user gives ISE only, it's ISE only
% - if the user gives EE only, it's EE only
% - if the user gives both, it's both
if ~opts.EE && ~opts.ISE
    opts.EE  = true;
end

% Read and assemble the z inputs if provided
fprintf('Reading Z nuisance variables.\n');
if Nz > 0
    plm.Zset = cell(1,1);
    for m = 1:numel(plm.Zset)
        plm.Zset{m} = zeros(plm.N,1);
    end
end

if Nz == 0 && ~ opts.zevperdat
    plm.Zset{1} = []; %ones(plm.N,1);
%     opts.EE     = false;
%     opts.ISE    = true;
elseif Nz > 0
    plm.nZ = numel(plm.Zset);
    for z = 1:Nz % The maximum Nz is one
        Ztmp = palm_miscread(opts.z{z},[],[],opts.precision);
        plm.Zset{z} = Ztmp.data;
        if ~ isempty(plm.subjidx) && size(plm.Zset{z},1) ~= plm.N
            plm.Zset{z} = plm.Zset{z}(plm.subjidx,:);
        end
        if size(plm.Zset{z},1) ~= plm.N
            error([
                'The number of rows in the nuisance covariate matrix does\n' ...
                'not match the number of observations in the data.\n' ...
                '- Rows in the matrix: %d\n' ...
                '- Observations in the data: %d\n' ...
                'In file %s\n'], ...
                size(plm.Zset{z},1),plm.N,opts.z{z});
        end
        if any(isnan(plm.Zset{z}(:))) || any(isinf(plm.Zset{z}(:)))
            error([
                'The nuisance covariate matrix cannot contain NaN or Inf.\n' ...
                'In file %s\n'],opts.z{z});
        end
    end
end

% Include the ZEV per datum (assume Nz is one for now)
if opts.zevperdat
    for ev = 1:plm.nZEVdat
        disp('ndims')
        ndims(plm.Zset{1})
        nzcols = size(plm.Zset{1},2); % number of nuisance variables (columns in -z)
        if ndims(plm.Zset{1}) == 2 %#ok<ISMAT>
            plm.Zset{1} = ...
                repmat(plm.Zset{1},[1 1 plm.ZEVsiz(ev)]);
        end
        permute(plm.ZEVset{ev},[1 3 2])
        plm.Zset{1}(:,nzcols+ev,:) = ...
            permute(plm.ZEVset{ev},[1 3 2]);
    end
    plm = rmfield(plm,{'ZEVset','ZEVsiz'});
else 
    % Make Zset the right size if no zevperdat 
    if ndims(plm.Zset{1}) == 2 %#ok<ISMAT>
        plm.Zset{1} = ...
            repmat(plm.Zset{1},[1 1 size(plm.Yset{1},2)]);
    end
end
    
% Read and assemble the w inputs if provided
fprintf('Reading W nuisance variables.\n');
if Nw > 0
    plm.Wset = cell(1,1);
    for m = 1:numel(plm.Wset)
        plm.Wset{m} = zeros(plm.N,1);
    end
end

if Nw == 0 && ~ opts.wevperdat
    plm.Wset{1} = []; % ones(plm.N,1);
%     opts.EE     = false;
%     opts.ISE    = true;
elseif Nw > 0
    plm.nW = numel(plm.Wset);
    for w = 1:Nw % Maximum Nw is one
        Wtmp = palm_miscread(opts.w{w},[],[],opts.precision);
        plm.Wset{w} = Wtmp.data;
        if ~ isempty(plm.subjidx) && size(plm.Wset{z},1) ~= plm.N
            plm.Wset{w} = plm.Wset{z}(plm.subjidx,:);
        end
        if size(plm.Wset{z},1) ~= plm.N
            error([
                'The number of rows in the W nuisance covariate matrix does\n' ...
                'not match the number of observations in the data.\n' ...
                '- Rows in the matrix: %d\n' ...
                '- Observations in the data: %d\n' ...
                'In file %s\n'], ...
                size(plm.Wset{w},1),plm.N,opts.w{w});
        end
        if any(isnan(plm.Wset{w}(:))) || any(isinf(plm.Wset{w}(:)))
            error([
                'The nuisance covariate matrix cannot contain NaN or Inf.\n' ...
                'In file %s\n'],opts.w{w});
        end
    end
end

% Include the WEV per datum (assume Nw is one for now)
if opts.wevperdat
    for ev = 1:plm.nWEVdat
        disp('ndims')
        ndims(plm.Wset{1})
        nwcols = size(plm.Wset{1},2); % number of nuisance variables (columns in -z)
        if ndims(plm.Wset{1}) == 2 %#ok<ISMAT>
            plm.Wset{1} = ...
                repmat(plm.Wset{1},[1 1 plm.WEVsiz(ev)]);
        end
        permute(plm.WEVset{ev},[1 3 2])
        plm.Wset{1}(:,nwcols+ev,:) = ...
            permute(plm.WEVset{ev},[1 3 2]);
    end
    plm = rmfield(plm,{'WEVset','WEVsiz'});
else 
    % Make Wset the right size if no wevperdat 
    if ndims(plm.Wset{1}) == 2 %#ok<ISMAT>
        plm.Wset{1} = ...
            repmat(plm.Wset{1},[1 1 size(plm.Xset{1},2)]);
    end
end

% pause(10)
% Some related sanity checks
% if opts.evperdat
%     for m = 1:plm.nM
%         if opts.designperinput, loopY = m; else, loopY = 1:plm.nY; end
%         for y = loopY
%             if size(plm.Yset{y},2) ~= size(plm.Mset{m},3)
%                 error([
%                     'The size of the data and the size of the EV per datum\n' ...
%                     'don''t match.%s'],'')
%             end
%         end
%     end
% end

% Make sure EVs of interest aren't represented also in the nuisance
% Note that some lines below the same is done for the cases in which 
% data and design are mean-centered
if ~ opts.demean && ~ opts.vgdemean && ~ opts.noranktest
    testrank(plm)
end

if ~ opts.cmcx
    seqtmp = zeros(plm.N,1);
    j = 1;
    plm.seq = cell(1,1);
    for m = 1:1
        plm.seq{m} = cell(1,1);
        if opts.cca.do
            Xtmp = plm.Xset{1}(:,:,1);
            Ztmp = plm.Zset{1}(:,:,1);
            Xtmp   = Xtmp(1:size(Xtmp,1)-rank(Ztmp),:);
            %seqtmp = seqtmp(1:size(Xtmp,1),:);
            [~,~,plm.seq{1}{1}] = unique(Xtmp,'rows');
%             seqtmp(:,j) = plm.seq{m}{c};
%             j = j + 1;
        end
    end
    tmp = sum(diff(seqtmp,1,2).^2,2);
    if (opts.corrcon || opts.npccon || opts.syncperms) && any(tmp(:) ~= 0)
        warning([ ...
            'You chose to correct over contrasts, or run NPC\n'             ...
            '         between contrasts, but with the design(s) and,\n'     ...
            '         contrasts given it is not possible to run\n'          ...
            '         synchronised permutations without ignoring repeated\n'...
            '         elements in the design matrix (or matrices). To\n'    ...
            '         solve this, adding the option "-cmcx" automatically.%s\n'],'');
        opts.cmcx = true;
    end
end

if opts.cmcx % note can't use 'else' here because opts.cmcx is modified in the 'if' above
    plm.seq = cell(plm.nM,1);
    for m = 1:plm.nM
        plm.seq{m} = cell(plm.nC(m),1);
        for c = 1:plm.nC(m)
            plm.seq{m}{c} = (1:plm.N)';
        end
    end
end

% Make sure not too many components are asked if CCA or PLS are used
if opts.cca.do || opts.pls.do
    if opts.ccaorplsparm > plm.nY
        error(['Cannot ask more canonical correlations (CCA) or \n', ...
            'score vectors (PLS) (k=%d) than the number of modalities (#(-i)=%d).\n'],...
            opts.ccaorplsparm,plm.nY);
    end
end

% Read the exchangeability blocks. If none is specified, all observations
% are assumed to be in the same large block. Also treat the legacy format of
% a single column for the EBs.
if isempty(opts.eb)
    plm.EB = [];
    if opts.within || opts.whole
        error([ ...
            'Options -within and/or -whole require a file defining\n' ...
            '         the exchangeability blocks (option -eb).\n%s'],'');
    end
else
    plm.EB = palm_miscread(opts.eb);
    plm.EB = plm.EB.data;
    if ~ isempty(plm.subjidx)
        plm.EB = plm.EB(plm.subjidx,:);
    end
    if isvector(plm.EB)
        if opts.within && opts.whole % within + whole block shuffling
            plm.EB = [+ones(plm.N,1) +plm.EB(:) (1:plm.N)'];
        elseif opts.whole             % whole-block shuffling
            plm.EB = [+ones(plm.N,1) -plm.EB(:) (1:plm.N)'];
        else                          % within-block shuffling (this is the default, and not meant to be changed)
            plm.EB = [-ones(plm.N,1) +plm.EB(:) (1:plm.N)'];
        end
    elseif opts.within || opts.whole
        warning([ ...
            'Options -within and/or -whole ignored, as the file defining\n' ...
            '         the exchangeability blocks (option -eb) already \n' ...
            '         defines how the data should be shuffled.%s'],'');
    elseif ~ isempty(plm.subjidx)
        error('Cannot use "-subjidx" with multi-level blocks.');
    end
    plm.EB = palm_reindex(plm.EB,'fixleaves');
end

% Mean center input data (-demean)
if opts.demean
    for m = 1:numel(plm.Xset)
        plm.Xset{m} = bsxfun(@minus,plm.Xset{m},mean(plm.Xset{m},1));
    end
    for y = 1:plm.nY
        plm.Yset{y} = bsxfun(@minus,plm.Yset{y},mean(plm.Yset{y},1));
    end
    if Nz >0
        for m = 1:Nz
            plm.Zset{m} = bsxfun(@minus,plm.Zset{m},mean(plm.Zset{m},1));
        end
    end
    if Nw > 0
        for m = 1:Nw
            plm.Wset{m} = bsxfun(@minus,plm.Wset{m},mean(plm.Wset{m},1));
        end
    end
end

% Make sure EVs of interest aren't represented also in the nuisance
% Note that some lines above the same is done for the cases in which there
% is no mean-centering
if (opts.demean || opts.vgdemean) && ~ opts.noranktest
    testrank(plm)
end

% ==============================================================
function testrank(plm)
% Test the complicated case in which design is rank deficient, with
% redundant dimensions being equally represented in EVs of interest and in
% nuisance.
plm.nM=1; plm.nC=1;
% badcon = cell(plm.nM,1);
% for m = 1:plm.nM
%     badcon{m} = zeros(plm.nC(m),1);
%     for c = 1:plm.nC(m)
%         [Xgutt,Zgutt] = palm_partition(plm.Mset{m}(:,:,1),plm.Cset{m}{c},'Guttman');
%         if size(Zgutt,2) > 0
%             cc = palm_cca(Xgutt,Zgutt,plm.N);
%             if cc(1) == 1
%                 badcon{m}(c) = sum(cc == 1);
%             end
%         end
%     end
% end
% idxbad = vertcat(badcon{:});
% if any(idxbad)
%     badlist = zeros(2,sum(idxbad > 0));
%     j = 1;
%     for m = 1:plm.nM
%         for c = 1:plm.nC(m)
%             if badcon{m}(c)
%                 badlist(1,j) = m;
%                 badlist(2,j) = c;
%                 badlist(3,j) = badcon{m}(c);
%                 j = j + 1;
%             end
%         end
%     end
%     badmsg = sprintf('- Design %d, Contrast %d, at least %d regressors\n',badlist);
%     error([...
%         'The following contrasts try to test regressor(s) also fully represented by\n'...
%         'by nuisance variable(s), but such tests are not possible (rank deficiency):\n%s'],badmsg); %#ok<SPERR>
% end

% ==============================================================
function checkmiss(A,Afname,N)
% Check if the missing data indicators are sane.
for a = 1:size(A,2)
    U = unique(A(:,a));
    if size(A,1) ~= N
        error([ ...
            'The missing data indicators ("-imiss" and "-dmiss") must have\n', ...
            'the same number of observations as the data and design.'],'');
    elseif ...
            (numel(U) >  2) || ...
            (numel(U) == 2 && ~ any(U == 0)) || ...
            (numel(U) == 2 && ~ any(U(U~=0) == [-1 1 2])) || ...
            (numel(U) == 1 && U ~= 0)
        error([ ...
            'The missing data indicators ("-imiss" and "-dmiss") must have\n', ...
            'no more than two unique values per column, one being 0, the\n', ...
            'the other being either -1, 1, or 2.\n', ...
            'Consult the documentation for details.\n'...
            '- Filename: %s\n',...
            '- Column: %d (possibly also others)'],Afname,a);
    end
end
