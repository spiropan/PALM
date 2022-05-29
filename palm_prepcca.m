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

% Read and organise the data.
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
    [plm.Yset{i},plm.masks{i},plm.Yisvol(i),plm.Yissrf(i),plm.Ykindstr{i}] = ...
        palm_ready(opts.y{i},maskstruct,opts,true);
    
    
    % For the first input data, keep the size to
    % compare with the others, then check the size
    if i == 1
        plm.N = size(plm.Yset{i},1);
    end
    if size(plm.Yset{i},1) ~= plm.N
        error([
            'At least two of the input data files do not have\n' ...
            'compatible sizes:\n' ...
            '- File %d (%s) has %d observations\n' ...
            '- File %d (%s) has %d observations'], ...
            1,opts.y{1},plm.N, ...
            i,opts.y{i},size(plm.Yset{i},1));
    end
    
    % If this is a DPX/curvature file, and if one of the spatial
    % statistics has been invoked, check if surfaces are available
    % and with compatible size, then compute the area (dpv or dpf).
    % Also take this opportunity to compute the adjacency matrix.
    if opts.tableasvolume
        plm.Yisvol(i) = true;
    else
        if opts.spatial.do && Ns > 0
            
            
            % String defining the types, for the filenames and other tasks.
            if Ns == 1, s = 1; else, s = i; end
            if any(size(plm.srf{s}.data.vtx,1) == ...
                    size(plm.masks{i}.data))
                plm.Yisvol(i)   = false;
                plm.Yissrf(i)   = true;
                plm.Yisvtx(i)   = true;
                plm.Yisfac(i)   = false;
                plm.Ykindstr{i} = '_dpv';
            elseif any(size(plm.srf{s}.data.fac,1) == ...
                    size(plm.masks{i}.data))
                plm.Yisvol(i)   = false;
                plm.Yissrf(i)   = true;
                plm.Yisvtx(i)   = false;
                plm.Yisfac(i)   = true;
                plm.Ykindstr{i} = '_dpf';
            else
                error([...
                    'Surface file does not match the input data:\n' ...
                    '- Surface file %d has %d vertices and %d faces (%s)\n' ...
                    '- Input data file %d has %d points (%s)'],...
                    s,size(plm.srf{s}.data.vtx,1),size(plm.srf{s}.data.fac,1),opts.s{s},...
                    i,max(size(plm.masks{i}.data)),opts.i{i});
            end
            
            % Surface area, to be used by the spatial statistics
            if isempty(plm.srfarea{s}.data)
                
                % No area file given, use the actual surface area
                plm.Yarea{i} = palm_calcarea(plm.srf{s}.data,plm.Yisvtx(i));
                
            elseif numel(plm.srfarea{s}.data) == 1
                
                % A weight given (such as 1): use that for each vertex or face,
                % treating as if all had the same area.
                if plm.Yisvtx(i)
                    plm.Yarea{i} = plm.srfarea{s}.data .* ...
                        ones(size(plm.srf{s}.data.vtx,1),1);
                elseif plm.Yisfac(i)
                    plm.Yarea{i} = plm.srfarea{s}.data .* ...
                        ones(size(plm.srf{s}.data.fac,1),1);
                end
                
            else
                
                % Otherwise, just use the data from the file (already loaded).
                plm.Yarea{i} = plm.srfarea{s}.data;
            end
            
            % Compute the adjacency matrix
            plm.Yadjacency{i} = palm_adjacency(plm.srf{s}.data.fac,plm.Yisvtx(i));
            
        elseif opts.spatial.do && Ns == 0 && any(strcmpi(plm.masks{i}.readwith,{'load','fs_load_mgh','gifti'}))
            error([ ...
                'To use spatial statistics with vertexwise or facewise data it is\n'...
                'necessary to provide the surface files (with the option "-s").%s'],'');
        end
    end
end


plm.nmasks = numel(plm.masks);
plm.nY     = numel(plm.Yset); 
plm.nX     = numel(plm.Xset);

% Select subjects
if ~ isempty(plm.subjidx)
    plm.Yset{i} = plm.Yset{i}(plm.subjidx,:);
    plm.Xset{i} = plm.Xset{i}(plm.subjidx,:);
end

% Some extra packages for Octave
if palm_isoctave
    if opts.spatial.do && any(plm.Yisvol)
        pkg load image
    end
    if opts.accel.lowrank || opts.zstat || opts.corrcon
        pkg load statistics
    end
end

% Read the EV per datum:
if opts.evperdat
    
    % Some sanity check:
    opts.evpos = cat(1,opts.evpos{:});
    if size(unique(opts.evpos,'rows'),1) ~= size(opts.evpos,1)
        error([
            'Some EV per datum have been defined for the same\n'...
            'position in the same design matrices.%s'],'');
    end
    plm.EVset  = cell(Nevd,1);
    plm.nEVdat = Nevd;
    
    % If there's one design per input, use the same masks as
    % those of the input files. Otherwise, create them.
    if ~ opts.designperinput && Nm == 1
        for ev = 1:plm.nEVdat
            plm.masksEV{ev} = plm.masks{1};
        end
    elseif opts.designperinput || (plm.nY == 1 && Nd == 1)
        for ev = 1:plm.nEVdat
            plm.masksEV{ev} = plm.masks{opts.evpos(ev,2)};
        end
    else
        plm.masksEV = cell(plm.nEVdat,1);
    end
    
    % Read input file & select subjects
    for ev = 1:plm.nEVdat
        fprintf('Reading EV per datum %d/%d: %s\n',ev,plm.nEVdat,opts.evdatfile{ev});
        [plm.EVset{ev},plm.masksEV{ev}] = palm_ready(opts.evdatfile{ev},plm.masksEV{ev},opts,opts.demean);
        if ~ isempty(plm.subjidx)
            plm.EVset{ev} = plm.EVset{ev}(plm.subjidx,:);
        end
    end
    plm.nmasksEV = numel(plm.masksEV);
    
    % Make the intersection of the EVs that will go all in the same design
    Dlist = unique(opts.evpos(:,2),'rows');
    for d = 1:numel(Dlist,1)
        evidx = find(opts.evpos(:,2) == Dlist(d));
        if numel(evidx) > 1
            newmask = true(size(plm.masksEV{1}.data));
            for ev = evidx'
                newmask = newmask & plm.masksEV{ev}.data;
            end
            for ev = evidx'
                plm.EVset{ev} = plm.EVset{ev}(:,newmask(plm.masksEV{ev}.data(:)));
                plm.masksEV{ev}.data = newmask;
            end
        end
    end
    
    % Then make the intersections with the respective masks of the input files
    if (opts.designperinput || (plm.nY == 1 && Nd == 1)) && ...
            ~ opts.npcmod && ~ opts.npccon && ~ opts.MV
        for d = unique(opts.evpos(:,2))'
            EV = find(opts.evpos(:,2) == d);
            newmask = zeros(numel(plm.masks{d}.data),numel(EV)+1);
            c = 1;
            for ev = EV'
                newmask(:,c) = plm.masksEV{ev}.data(:);
                c = c + 1;
            end
            newmask(:,c) = plm.masks{d}.data(:);
            newmask = reshape(all(newmask,2),size(plm.masks{d}.data));
            for ev = EV'
                plm.EVset{ev} = plm.EVset{ev}(:,newmask(plm.masksEV{ev}.data(:)));
                plm.masksEV{ev}.data = newmask;
            end
            plm.Yset{d} = plm.Yset{d}(:,newmask(plm.masks{d}.data(:)));
            plm.masks{d}.data    = newmask;
        end
    else
        sizy = zeros(plm.nmasks,1);
        for y = 1:plm.nmasks
            sizy(y) = numel(plm.masks{y}.data);
        end
        sizev = zeros(plm.nmasksEV,1);
        for ev = 1:plm.nmasksEV
            sizev(ev) = numel(plm.masksEV{ev}.data);
        end
        if  numel(unique(sizy)) > 1 || ...
                numel(unique(sizev)) > 1 || ...
                sizy(1) ~= sizev(1)
            error([...
                'For multiple "-i" and/or "-evperdat", with "-npccon",\n',...
                'and without the option "-designperinput", the inputs, EVs \n'...
                'and masks need to be all of the same sizes.%s'],'');
        end
        newmask = true(size(plm.masksEV{1}.data));
        for y = 1:plm.nmasks
            newmask = newmask & plm.masks{y}.data;
        end
        for ev = 1:plm.nmasksEV
            newmask = newmask & plm.masksEV{ev}.data;
        end
        for y = 1:plm.nmasks
            plm.Yset{y} = plm.Yset{y}(:,newmask(plm.masks{y}.data(:)));
            plm.masks{y}.data = newmask;
        end
        for ev = 1:plm.nmasksEV
            plm.EVset{ev} = plm.EVset{ev}(:,newmask(plm.masksEV{ev}.data(:)));
            plm.masksEV{ev}.data = newmask;
        end
    end
end
clear('newmask');

% Create an intersection mask if NPC or MV is to be done, and further apply
% to the data that was previously masked above, as needed.
if opts.npcmod || opts.MV || opts.cca.do || opts.pls.do || opts.forcemaskinter
    if plm.nmasks > 1
        
        % If there is one mask per modality, make an instersection mask.
        maskinter = true(size(plm.masks{1}.data));
        for m = 1:plm.nmasks
            maskinter = maskinter & plm.masks{m}.data;
        end
        if opts.evperdat
            for ev = 1:plm.nEVdat
                maskinter = maskinter & plm.masksEV{ev}.data;
            end
        end
        
        % Note that this line below uses Ytmp, which is from the previous loop.
        % This can be used here because with NPC all data has the same size.
        plm.maskinter = palm_maskstruct(maskinter(:)',plm.masks{1}.readwith,plm.masks{1}.extra);
        
        % Apply it to further subselect data points
        for y = 1:plm.nY
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{y}.data));
        end
        if opts.evperdat
            for ev = 1:plm.nEVdat
                plm.EVset{ev} = plm.EVset{ev}(:,plm.maskinter.data(plm.masksEV{ev}.data));
            end
        end
    else
        % If only one mask was given.
        plm.maskinter = plm.masks{1};
        for y = 1:plm.nY
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{1}.data));
        end
    end
end

% Sizes for later
if opts.evperdat
    plm.EVsiz = zeros(plm.nEVdat,1);
    for ev = 1:plm.nEVdat
        plm.EVsiz(ev) = size(plm.EVset{ev},2);
    end
end

% Make sure that all data have the same size if NPC or MV are to be done
if opts.npcmod || opts.MV || opts.forcemaskinter
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

% Make sure none of the modalities is empty
for y = 1:plm.nY
    if any(size(plm.Yset{y}) == 0)
        error('Modality %d has no data.\n',y);
    end
end

% If the multiple columns of the (sole) input are to be treated
% in a multivariate fashion
if opts.inputmv
    tmp1 = plm.Yset{1};
    tmp2 = plm.Ykindstr{1};
    nY   = size(tmp1,2);
    plm.Yset     = cell(nY,1);
    plm.Ykindstr = cell(nY,1);
    for y = 1:nY
        plm.Yset{y}     = tmp1(:,y);
        plm.Ykindstr{y} = tmp2;
    end
    clear tmp1 tmp2 nY;
    plm.nY = numel(plm.Yset);
end

% A variable with the cumulative sizes of all modalities will be handy later
plm.Ysiz = zeros(plm.nY,1);
for y = 1:plm.nY
    plm.Ysiz(y) = size(plm.Yset{y},2);
end
plm.Ycumsiz = vertcat(0,cumsum(plm.Ysiz));

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

% If MV was selected, make sure that Y is full rank.
if opts.MV && ~ opts.noranktest
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

% Read and assemble the design matrices.
fprintf('Reading design matrix and contrasts.\n');
if opts.evperdat
    plm.Mset = cell(max(Nd,max(opts.evpos(:,2))),1);
    for m = 1:numel(plm.Mset)
        plm.Mset{m} = zeros(plm.N,1);
    end
else
    plm.Mset = cell(max(Nd,1),1);
end
plm.nM = numel(plm.Mset);
if Nd == 0 && ~ opts.evperdat
    plm.Mset{1} = ones(plm.N,1);
    opts.EE     = false;
    opts.ISE    = true;
elseif Nd > 0
    for m = 1:Nd
        Mtmp = palm_miscread(opts.d{m},[],[],opts.precision);
        plm.Mset{m} = Mtmp.data;
        if ~ isempty(plm.subjidx) && size(plm.Mset{m},1) ~= plm.N
            plm.Mset{m} = plm.Mset{m}(plm.subjidx,:);
        end
        if size(plm.Mset{m},1) ~= plm.N
            error([
                'The number of rows in the design matrix does\n' ...
                'not match the number of observations in the data.\n' ...
                '- Rows in the matrix: %d\n' ...
                '- Observations in the data: %d\n' ...
                'In file %s\n'], ...
                size(plm.Mset{m},1),plm.N,opts.d{m});
        end
        if any(isnan(plm.Mset{m}(:))) || any(isinf(plm.Mset{m}(:)))
            error([
                'The design matrix cannot contain NaN or Inf.\n' ...
                'In file %s\n'],opts.d{m});
        end
    end
end

% Include the EV per datum
if opts.evperdat
    for ev = 1:plm.nEVdat
        if ndims(plm.Mset{opts.evpos(ev,2)}) == 2 %#ok<ISMAT>
            plm.Mset{opts.evpos(ev,2)} = ...
                repmat(plm.Mset{opts.evpos(ev,2)},[1 1 plm.EVsiz(ev)]);
        end
        plm.Mset{opts.evpos(ev,2)}(:,opts.evpos(ev,1),:) = ...
            permute(plm.EVset{ev},[1 3 2]);
    end
    plm = rmfield(plm,{'EVset','EVsiz'});
end

% Some related sanity checks
if opts.evperdat
    for m = 1:plm.nM
        if opts.designperinput, loopY = m; else, loopY = 1:plm.nY; end
        for y = loopY
            if size(plm.Yset{y},2) ~= size(plm.Mset{m},3)
                error([
                    'The size of the data and the size of the EV per datum\n' ...
                    'don''t match.%s'],'')
            end
        end
    end
end

% Read and organise the contrasts for each design.
plm.Cset = cell(plm.nM,1);
plm.Dset = cell(plm.nM,1);
plm.nC   = zeros(plm.nM,1);
plm.nD   = zeros(plm.nM,1);
if Nt || Nf
    
    % Load FSL style t contrasts
    tcon = cell(Nt,1);
    for t = 1:Nt
        tmp = palm_miscread(opts.t{t},[],[],opts.precision);
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'}))
            tcon{t} = tmp.data;
        else
            error('Invalid t contrast file: %s',opts.t{t});
        end
    end
    
    % Load FSL style F contrasts
    fcon = cell(Nt,1);
    for t = 1:Nt
        if ~ isempty(opts.f{t})
            tmp = palm_miscread(opts.f{t},[],[],opts.precision);
            if any(strcmp(tmp.readwith,{'vestread','csvread','load'}))
                fcon{t} = tmp.data;
            else
                error('Invalid F contrast file: %s',opts.f{t});
            end
        end
    end
    
    % For each valid design, assemble the contrasts.
    for m = 1:plm.nM
        if Nt == 1
            t = 1;
        else
            t = m;
        end
        c = 1;
        for j = 1:size(tcon{t},1)
            plm.Cset{m}{c} = tcon{t}(j,:)';
            c = c + 1;
        end
        for j = 1:size(fcon{t},1)
            if ~ isempty(fcon{t})
                plm.Cset{m}{c} = tcon{t}(logical(fcon{t}(j,:)),:)';
                c = c + 1;
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        for c = 1:plm.nC(m)
            plm.Dset{m}{c} = eye(plm.nY);
        end
        plm.nD(m) = numel(plm.Dset{m});
    end
    
elseif Ncon
    
    % Load MSET style contrasts (all contrast pairs)
    Ccon = cell(Ncon,1);
    Dcon = cell(Ncon,1);
    for con = 1:Ncon
        tmp = palm_miscread(opts.Ccon{con},[],[],opts.precision);
        if strcmpi(tmp.readwith,'mset')
            Ccon{con} = tmp.data;
        else
            error(['Files given to the option "-con" must be in .mset format.\n' ...
                'For .csv/.con/.fts files, use "-t" or "-f".%s'],'');
        end
        if isempty(opts.Dcon{con})
            for c = 1:numel(Ccon{con})
                Dcon{con}{c} = eye(plm.nY);
            end
        else
            tmp = palm_miscread(opts.Dcon{con},[],[],opts.precision);
            if strcmpi(tmp.readwith,'mset')
                Dcon{con} = tmp.data;
            else
                error(['Files given to the option "-con" must be in .mset format.\n' ...
                    'For .csv/.con/.fts files, use "-t" or "-f".%s'],'');
            end
        end
    end
    
    % Assign the contrast sets to the design matrices
    for m = 1:plm.nM
        if Ncon == 1
            con = 1;
        else
            con = m;
        end
        plm.Cset{m} = Ccon{con};
        plm.Dset{m} = Dcon{con};
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
        for c = 1:plm.nC(m)
            plm.Cset{m}{c} = plm.Cset{m}{c}';
        end
        for d = 1:plm.nD(m)
            plm.Dset{m}{d} = plm.Dset{m}{d}';
        end
    end
else
    % If no constrasts were at all specified:
    for m = 1:plm.nM
        if size(plm.Mset{m},2) == 1
            
            % If there is only 1 regressor, test its effect both
            % positive and negative.
            % The statistic will be t or v, depending on the number of VGs.
            plm.Cset{m}{1} = 1;
            plm.Cset{m}{2} = -1;
            plm.Dset{m}{1} = eye(plm.nY);
            plm.Dset{m}{2} = eye(plm.nY);
        else
            % Otherwise, run an F-test over all regressors in the design matrix.
            % The statistic will be F or G, depending on the number of VGs.
            plm.Cset{m}{1} = eye(size(plm.Mset{m},2));
            plm.Dset{m}{1} = eye(plm.nY);
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
end

% Ranks of the contrasts
plm.rC = cell(plm.nM,1);
plm.rD = plm.rC;
for m = 1:plm.nM
    plm.rC{m} = zeros(plm.nC(m),1);
    plm.rD{m} = plm.rC{m};
    for c = plm.nC(m):-1:1
        rC = rank(plm.Cset{m}{c});
        rD = rank(plm.Dset{m}{c});
        if rC == 0 || rD == 0
            plm.Cset{m}(c) = [];
            plm.Dset{m}(c) = [];
            plm.rC{m}(c)   = [];
            plm.rD{m}(c)   = [];
            warning('For design %d, contrast %d has rank 0 and will be removed.',m,c);
        else
            plm.rC{m}(c) = rC;
            plm.rD{m}(c) = rD;
        end
    end
    plm.nC(m) = numel(plm.Cset{m});
    plm.nD(m) = numel(plm.Dset{m});
end
plm.rC0 = plm.rC; % the rC can be changed for z-scores, but not rC0.

% If only the t or F tests are to be performed
if opts.tonly
    for m = 1:plm.nM
        for c = plm.nC(m):-1:1
            if plm.rC{m}(c) > 1
                plm.Cset{m}(c) = [];
                plm.Dset{m}(c) = [];
                plm.rC{m}(c)   = [];
                plm.rC0{m}(c)  = [];
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
elseif opts.fonly
    for m = 1:plm.nM
        for c = plm.nC(m):-1:1
            if plm.rC{m}(c) == 1
                plm.Cset{m}(c) = [];
                plm.Dset{m}(c) = [];
                plm.rC{m}(c)   = [];
                plm.rC0{m}(c)  = [];
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
end

% Some more sanity checks
for m = 1:plm.nM
    for c = 1:plm.nC(m)
        if any(isnan(plm.Cset{m}{c}(:))) || any(isinf(plm.Cset{m}{c}(:)))
            error('The constrasts cannot contain NaN or Inf.');
        end
        if size(plm.Cset{m}{c},1) ~= size(plm.Mset{m},2)
            error('The size of one or more contrasts don''t match the size of the respective design matrix.')
        end
    end
    for c = 1:plm.nD(m)
        if any(isnan(plm.Dset{m}{c}(:))) || any(isinf(plm.Dset{m}{c}(:)))
            error('The constrasts cannot contain NaN or Inf.');
        end
    end
end
if opts.MV
    if any(plm.nC ~= plm.nD)
        error('The number of C and D contrasts must be the same');
    end
    for m = 1:plm.nM
        for d = 1:plm.nD(m)
            if size(plm.Dset{m}{d},1) ~= plm.nY
                error('The size of one or more MV contrasts don''t match the number of modalities.');
            end
        end
    end
end
for m = 1:plm.nM
    for c = 1:plm.nC(m)
        if opts.concordant && plm.rC{m}(c) > 1
            error(['Cannot use the "-concordant" option with F-tests (inadmissible).\n'...
                'Use "-tonly" to run just the t-tests, or remove "-concordant".%s'],'');
        end
    end
end

% Make sure EVs of interest aren't represented also in the nuisance
% Note that some lines below the same is done for the cases in which 
% data and design are mean-centered
if ~ opts.demean && ~ opts.vgdemean && ~ opts.noranktest
    testrank(plm)
end

% Check if the contrasts have all the same rank for correction over
% contrasts. If not, convert to zstat.
if opts.corrcon && ~ opts.zstat
    rC1 = plm.rC{1}(1);
    rD1 = plm.rD{1}(1);
    for m = 1:plm.nM
        brflag = false;
        for c = 1:plm.nC(m)
            if rC1 ~= plm.rC{m}(c) || rD1 ~= plm.rD{m}(c)
                warning([...
                    'Not all contrasts have the same rank, and the option "-corrcon" was used.\n' ...
                    '         Adding the option "-zstat" automatically.%s'],'');
                opts.zstat = true;
                brflag = true;
                break;
            end
        end
        if brflag
            break;
        end
    end
end

% Partition the model according to the contrasts and design matrix.
% The partitioning needs to be done now, because of the need for
% synchronised permutations/sign-flips
if opts.corrcon || opts.npccon
    opts.syncperms = true;
end
if ~ opts.cmcx
    seqtmp = zeros(plm.N,sum(plm.nC));
    j = 1;
    plm.seq = cell(plm.nM,1);
    for m = 1:plm.nM
        plm.seq{m} = cell(plm.nC(m),1);
        for c = 1:plm.nC(m)
            Xtmp = palm_partition(plm.Mset{m}(:,:,1),plm.Cset{m}{c},opts.pmethodp);
            [~,~,plm.seq{m}{c}] = unique(Xtmp,'rows');
            seqtmp(:,j) = plm.seq{m}{c};
            j = j + 1;
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
if opts.CCA || opts.PLS
    if opts.ccaorplsparm > plm.nY
        error(['Cannot ask more canonical correlations (CCA) or \n', ...
            'score vectors (PLS) (k=%d) than the number of modalities (#(-i)=%d).\n'],...
            opts.ccaorplsparm,plm.nY);
    end
    for m = 1:plm.nM
        for c = 1:plm.nC(m)
            if opts.ccaorplsparm > plm.rC{m}(c)
                error(['Cannot ask more canonical correlations (for CCA) or score \n', ...
                    'vectors (for PLS) than the rank of the contrast (k=%d > rank=%d).\n', ...
                    'Please check design %d, contrast %d.'],opts.ccaorplsparm,plm.rC{m}(c),m,c);
            end
        end
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

% Load/define the variance groups.
if opts.singlevg
    % If single VG, it's all ones
    plm.VG = ones(plm.N,1);
elseif strcmpi(opts.vg,'auto')
    if isempty(plm.EB)
        % If auto, but there are no exchangeability blocks, it's all ones too
        plm.VG = ones(plm.N,1);
    else
        % Generate an initial dependence tree, to be used to define variance groups.
        % The tree used for the permutations later require the design matrix, and
        % varies for each contrast -- all to be taken care of later.
        Ptree  = palm_tree(plm.EB,(1:plm.N)');
        plm.VG = palm_ptree2vg(Ptree);
    end
else
    % The automatic variance groups can be overriden if the user specified
    % a file with the custom definitions.
    plm.VG = palm_miscread(opts.vg);
    plm.VG = plm.VG.data;
end
if ~ isempty(plm.subjidx) && size(plm.VG,1) ~= plm.N
    plm.VG = plm.VG(plm.subjidx,:);
end
[tmp,~,plm.VG] = unique(plm.VG);
plm.nVG = numel(tmp);
if plm.nVG == 1, opts.singlevg = true; end

% MV can't yet be used if nVG>1, although NPC remains an option
if opts.MV && plm.nVG > 1
    error('There are more than one variance group. MV cannot be used (but NPC can).');
end

% There should be no more missing indicators than modalities, or designs.
% These need to be either 1 or the same as the corresponding numbers of
% modalities/designs.
if Nimiss > Ni
    error([...
        'There are more missing indicators supplied with "-imiss" (%d) than\n'...
        'modalities supplied with "-i" (%d)'],Nimiss,Ni);
elseif Nimiss > 1 && Nimiss ~= Ni
    error([...
        'The number of missing indicators supplied with "-imiss" (%d) is larger,\n'...
        'than 1, but still not the same as the number of modalities supplied with\n'...
        'the option "-i" (%d).'],Nimiss,Ni);
end
if Ndmiss > Nd
    error([...
        'There are more missing indicators supplied with "-dmiss" (%d) than\n'...
        'designs supplied with "-d" (%d)'],Nimiss,Ni);
elseif Ndmiss > 1 && Ndmiss ~= Nd
    error([...
        'The number of missing indicators supplied with "-dmiss" (%d) is larger,\n'...
        'than 1, but still not the same as the number of modalities supplied with\n'...
        'the option "-d" (%d).'],Ndmiss,Nd);
end

% Load the missing indicators for the data ("imiss"):
for i = 1:Nimiss
    if strcmpi(opts.imiss{i},'none')
        if isempty(plm.subjidx)
            tmp = zeros(plm.N,1);
        else
            tmp = zeros(size(plm.subjidx,1),1);
        end
    else
        tmp = palm_miscread(opts.imiss{i});
        tmpfname = tmp.filename;
        tmp = tmp.data;
        if ~ isempty(plm.subjidx) && size(tmp,1) ~= plm.N
            tmp = tmp(plm.subjidx,:);
        end
        checkmiss(tmp,tmpfname,plm.N);
    end
    plm.Ymiss{i} = tmp;
end
if Nimiss == 1
    for i = 2:Ni
        plm.Ymiss{i} = plm.Ymiss{1};
    end
end

% Load the missing indicators for the design ("dmiss"):
for d = 1:Ndmiss
    if strcmpi(opts.dmiss{d},'none')
        if isempty(plm.subjidx)
            tmp = zeros(plm.N,1);
        else
            tmp = zeros(size(plm.subjidx,1),1);
        end
    else
        tmp = palm_miscread(opts.dmiss{d});
        tmpfname = tmp.filename;
        tmp = tmp.data;
        if ~ isempty(plm.subjidx) && size(tmp,1) ~= plm.N
            tmp = tmp(plm.subjidx,:);
        end
        checkmiss(tmp,tmpfname,plm.N);
    end
    plm.Mmiss{d} = tmp;
end
if Ndmiss == 1
    for d = 2:Nd
        plm.Mmiss{d} = plm.Mmiss{1};
    end
end
for d = 1:Ndmiss
    if any(size(plm.Mmiss{d}) ~= size(plm.Mset{d}))
        if strcmpi(opts.dmiss{d},'none')
            plm.Mmiss{d} = repmat(plm.Mmiss{d},[1 size(plm.Mset{d},2)]);
        else
            error([ ...
                'The missing data indicator ("-dmiss") must have\n', ...
                'the same size as the respective design.%s'],'');
        end
    end
end

% If only data or design missing indicators are missing, fill the other
% with all-false indicators.
if     Nimiss && ~ Ndmiss
    for m = 1:plm.nM
        plm.Mmiss{m} = false(size(plm.Mset{m}));
    end
elseif Ndmiss && ~ Nimiss
    for y = 1:plm.nY
        plm.Ymiss{y} = false(size(plm.Yset{y}));
    end
end

% Remove the variance groups with tiny sample sizes?
if plm.nVG > 1 && ~ opts.removevgbysize && (opts.vgdemean || opts.ev4vg) && ...
        any(sum(bsxfun(@eq,plm.VG,unique(plm.VG)'),1) == 1)
    warning([...
        'The options "-vgdemean" and "-ev4vg" require that observations\n' ...
        '         in variance groups of size 1 are removed.\n' ...
        '         Enabling the option "-removevgbysize 1"%s.'],'');
    opts.removevgbysize = 1;
end
if ~ opts.removevgbysize
    tmpwarned = false;
    for u = 1:plm.nVG
        if sum((plm.VG == u),1) == 1
            if ~ tmpwarned
                warning([...
                    'There are variance groups with just one observation.\n' ...
                    '         Consider using the option "-removevgbysize 1" to improve the\n' ...
                    '         variance estimates (at the cost of reducing sample size).%s'],'');
                tmpwarned = true;
            end
        end
    end
end
if opts.removevgbysize > 0
    
    % Indices of the observations to keep
    uVG = unique(plm.VG)';
    idxvg = sum(bsxfun(@eq,plm.VG,uVG),1) <= opts.removevgbysize;
    idx   = any(bsxfun(@eq,plm.VG,uVG(~idxvg)),2);
    
    % Modify all data as needed
    for y = 1:plm.nY
        plm.Yset{y} = plm.Yset{y}(idx,:);
    end
    if ~ isempty(plm.EB)
        plm.EB = plm.EB(idx,:);
    end
    for m = 1:plm.nM
        plm.Mset{m} = plm.Mset{m}(idx,:);
    end
    for m = 1:plm.nM
        for c = 1:plm.nC(m)
            plm.seq{m}{c} = plm.seq{m}{c}(idx,:);
        end
    end
    plm.N = sum(idx);
    [tmp,~,plm.VG] = unique(plm.VG(idx));
    plm.nVG = numel(tmp);
end

% Add one regressor for each variance group if requested
if opts.ev4vg
    for m = 1:plm.nM
        Mvg = zeros(plm.N,plm.nVG);
        V = unique(plm.VG);
        for v = 1:plm.nVG
            Mvg(plm.VG == V(v),v) = 1;
        end
        rM   = round(sum(diag(plm.Mset{m}*pinv(plm.Mset{m}))));
        Mnew = horzcat(plm.Mset{m},Mvg);
        if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG)
            plm.Mset{m} = Mnew;
            nadded      = plm.nVG;
        else
            Mnew = Mnew(:,1:end-1);
            if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG - 1)
                plm.Mset{m} = Mnew;
                nadded      = plm.nVG - 1;
            else
                error([ ...
                    'It was not possible to add one regressor for each variance group\n' ...
                    'perhaps because they already exist in the design. Check your design\n' ...
                    'matrix and maybe consider including these regressors manually.%s'],'');
            end
        end
        for c = 1:plm.nC(m)
            plm.Cset{m}{c} = vertcat(plm.Cset{m}{c},...
                zeros(nadded,size(plm.Cset{m}{c},2)));
        end
    end
end

% Make sure that mean-centering won't damage a contrast that tests the intercept.
if ~opts.cca.do && (opts.demean || opts.vgdemean)
    for m = 1:plm.nM
        for c = 1:plm.nC(m)
            % The design may be 3D here if -evperdat was used, hence can't
            % do matrix multiplication of design by contrast
            for t = 1:size(plm.Cset{m}{c},2)
                Xtmp = sum(bsxfun(@times,plm.Mset{m},plm.Cset{m}{c}(:,t)'),2);
                siz  = size(Xtmp);
                isintercp = all(bsxfun(@eq,reshape(Xtmp(1,:),[1 siz(2:end)]),Xtmp),1);
                if isintercp
                    error([ ...
                        'Contrast %d for design %d (and perhaps others) tests the intercept.\n' ...
                        'This means that the options "-demean" and "-vgdemean" cannot be used.\n' ...
                        'If "-demean" was added to calculate Pearson''s "r" or the "R^2"\n' ...
                        'note that these statistics cannot be computed for constant variables.%s'],c,m,'');
                end
            end
        end
    end
end

% Mean center data and design (-demean)
if opts.demean
    for m = 1:plm.nM
        plm.Mset{m} = bsxfun(@minus,plm.Mset{m},mean(plm.Mset{m},1));
    end
    for y = 1:plm.nY
        plm.Yset{y} = bsxfun(@minus,plm.Yset{y},mean(plm.Yset{y},1));
    end
end

% Mean center data and design, within VG
if opts.vgdemean
    
    % For each VG
    V = unique(plm.VG);
    for v = 1:plm.nVG
        vidx = plm.VG == V(v);
        
        % Demean design within VG
        for m = 1:plm.nM
            plm.Mset{m}(vidx,:) = bsxfun(@minus,...
                plm.Mset{m}(vidx,:),mean(plm.Mset{m}(vidx,:),1));
        end
        
        % Demean data within VG
        for y = 1:plm.nY
            plm.Yset{y}(vidx,:) = bsxfun(@minus,...
                plm.Yset{y}(vidx,:),mean(plm.Yset{y}(vidx,:),1));
        end
    end
end

% Make sure EVs of interest aren't represented also in the nuisance
% Note that some lines above the same is done for the cases in which there
% is no mean-centering
if (opts.demean || opts.vgdemean) && ~ opts.noranktest
    testrank(plm)
end

% Number of tests to be selected for the low rank approximation
if opts.accel.lowrank
    if plm.nVG > 1
        error('The option "-accel lowrank" cannot be used with more than one variance group.');
    end
    if opts.nP0 == 0
        error('With lowrank approximation you must indicate a larger-than-zero number of permutations.');
    end
    if opts.nP0 < plm.N*(plm.N+1)/2
        error([ ...
            'Too few permutations selected to use with lowrank approximation.\n' ...
            'Use at least N*(N+1)/2 = %d to note a speed difference and have reasonably accurate results.\n'...
            'Otherwise, don''t bother using lowrank approximation.\n'],plm.N*(plm.N+1)/2);
    end
    if opts.spatial.do
        warning([ ...
            'There isn''t much benefit in using lowrank approximation with spatial statistics\n' ...
            '         like TFCE and cluster extent and/or mass. These cannot be accelerated with this\n' ...
            '         method, and the overall gain will be minimal. Consider other approximation methods,\n' ...
            '         or run the full permutation test, or just drop spatial statistics.%s'],'');
    end
    plm.nsel = zeros(plm.nY,1);
    if opts.accel.lowrank_val <= 1
        for y = 1:plm.nY
            plm.nsel(y) = ceil(opts.accel.lowrank_val*plm.Ysiz(y));
        end
    elseif opts.accel.lowrank_val > 1
        plm.nsel(1:end) = ceil(opts.accel.lowrank_val);
    elseif isnan(opts.accel.lowrank_val)
        plm.nsel(1:end) = plm.N*(plm.N+1)/2;
    end
end

% ==============================================================
function testrank(plm)
% Test the complicated case in which design is rank deficient, with
% redundant dimensions being equally represented in EVs of interest and in
% nuisance.
badcon = cell(plm.nM,1);
for m = 1:plm.nM
    badcon{m} = zeros(plm.nC(m),1);
    for c = 1:plm.nC(m)
        [Xgutt,Zgutt] = palm_partition(plm.Mset{m}(:,:,1),plm.Cset{m}{c},'Guttman');
        if size(Zgutt,2) > 0
            cc = palm_cca(Xgutt,Zgutt,plm.N);
            if cc(1) == 1
                badcon{m}(c) = sum(cc == 1);
            end
        end
    end
end
idxbad = vertcat(badcon{:});
if any(idxbad)
    badlist = zeros(2,sum(idxbad > 0));
    j = 1;
    for m = 1:plm.nM
        for c = 1:plm.nC(m)
            if badcon{m}(c)
                badlist(1,j) = m;
                badlist(2,j) = c;
                badlist(3,j) = badcon{m}(c);
                j = j + 1;
            end
        end
    end
    badmsg = sprintf('- Design %d, Contrast %d, at least %d regressors\n',badlist);
    error([...
        'The following contrasts try to test regressor(s) also fully represented by\n'...
        'by nuisance variable(s), but such tests are not possible (rank deficiency):\n%s'],badmsg); %#ok<SPERR>
end

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
