function [opts,plm] = palm_permcca(opts,plm)
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
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Uncomment the line below for debugging:
%clear global plm opts; global plm opts;

opts.evperdat=true;
plm.nM=1; plm.nC=1;

% Variables to store stuff for later.
if opts.missingdata,
    nY = plm.nY; else nY = 1;
    plm.Ymissp = cell(plm.nY,1);
end
tmp = cell(nY,1);
for y = 1:nY,
    tmp{y} = cell(plm.nM,1);
    for m = 1:plm.nM,
        tmp{y}{m} = cell(plm.nC(m),1);
    end
end; clear('y');
plm.X        = tmp; % effective regressors
plm.Z        = tmp; % nuisance regressors
%plm.eCm      = tmp; % effective contrast (for Mp)
%plm.eCx      = tmp; % effective contrast (for the effective regressors only)
%plm.eC       = tmp; % final effective contrast (depends on the method)
%plm.Mp       = tmp; % partitioned model, joined
%plm.nEV      = tmp; % number of regressors
%plm.Hm       = tmp; % hat (projection) matrix
%plm.Rm       = tmp; % residual forming matrix
%plm.dRm      = tmp; % diagonal elements of the residual forming matrix
%plm.rM       = tmp; % rank of the design matrix
%plm.Gname    = cell(plm.nM,1); % name of the statistic for each contrast
plm.nP       = cell(plm.nM,1); % number of permutations for each contrast

for m = 1:plm.nM,
    plm.Gname{m} = cell (plm.nC(m),1);
    plm.nP{m}    = zeros(plm.nC(m),1);
end; clear('m');


Q        = cell(plm.nM,1);  % to store MV G at each permutation
% Lower levels of these variables
for m = 1:plm.nM,
    Q{m}     = cell(plm.nC(m),1);
end

% Variables for CCA
if opts.cca.do || opts.pls.do,
    plm.mvstr = ''; % default string for the filenames.
end

% Inital strings to save the file names later.
plm.ystr = cell(plm.nY,1);
for y = 1:plm.nY,
    plm.ystr{y} = '';
end
plm.mstr = cell(plm.nM,1);
plm.cstr = cell(plm.nM,1);
for m = 1:plm.nM,
    plm.mstr{m} = '';
    plm.cstr{m} = cell(plm.nC(m));
    for c = 1:plm.nC(m),
        plm.cstr{m}{c} = '';
    end
end
clear y m c;

% For each design matrix and contrast:
prepglm = cell(plm.nM,1);
fastpiv = cell(plm.nM,1);
for m = 1:plm.nM,
%     prepglm{m} = cell(plm.nC(m),1);
%     fastpiv{m} = cell(plm.nC(m),1);
    for c = 1:plm.nC(m),
        % CCA
        %%% DOUBLE-CHECK THE DEGREES-OF-FREEDOM!!
        if opts.cca.do,
            % Output string base, statistic function, and side to test
            plm.Qname{m}{c} = sprintf('_cca')
            plm.qfun = @cca;
            plm.mvrev{m}{c} = false;
        elseif opts.pls.do,
            % Output string, statistic function, and side to test
            plm.Qname{m}{c} = sprintf('_pls%d',opts.ccaorplsparm);
            plm.qfun = @simpls;
            plm.mvrev{m}{c} = false;
        end
        
        % If there are voxelwise EVs:
%         if opts.evperdat,
%             fprintf('Doing maths for -evperdat before model fitting: [Design %d/%d, Contrast %d/%d] (may take several minutes)\n',m,plm.nM,c,plm.nC(m));
%         end
%         
        % Partition the model, now using the method chosen by the user
        %if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
%         if opts.missingdata,
%             if opts.showprogress,
%                 fprintf('Preparing designs for missing data [Design: %d/%d, Contrast %d/%d]\n',m,plm.nM,c,plm.nC(m));
%             end
%             % Partition the design
%             for y = loopY,
%                 [plm.X{y}{m}{c},plm.Z{y}{m}{c},...
%                     plm.eCm{y}{m}{c},plm.eCx{y}{m}{c},...
%                     plm.Ymissp{y},...
%                     plm.imov{y}{m}{c},plm.ifix{y}{m}{c},...
%                     plm.isdiscrete{y}{m}{c},plm.istwotail{y}{m}{c}] = ...
%                     palm_misspart(plm.Mset{m},plm.Cset{m}{c},...
%                     opts.pmethodr,plm.Ymiss{y},plm.Mmiss{m},opts.mcar,opts.rmethod);
%                 for o = 1:numel(plm.X{y}{m}{c}),
%                     plm.Mp{y}{m}{c}{o} = cat(2,plm.X{y}{m}{c}{o},plm.Z{y}{m}{c}{o});
%                 end
%             end
%         else % not missing data
        y = m; o = 1;
            
        
        %end
%         for y = loopY,
%             if opts.missingdata, loopO = 1:numel(plm.Mp{y}{m}{c}); else loopO = 1; end
%             for o = loopO,
%                 
%                 % To avoid rank deficiency issues after partitioning, remove
%                 % columns that are all equal to zero. This won't be done for
%                 % evperdat because it's too slow and can make the designs too
%                 % different if EVs are dropped from just some of the tests.
%                 if ~ opts.evperdat,
%                     idx = all(plm.X{y}{m}{c}{o}  == 0,1);
%                     plm.X{y}{m}{c}{o}(:,idx)   = [];
%                     plm.eCx{y}{m}{c}{o}(idx,:) = [];
%                     idx = all(plm.Z{y}{m}{c}{o}  == 0,1);
%                     plm.Z{y}{m}{c}{o}(:,idx)   = [];
%                     idx = all(plm.Mp{y}{m}{c}{o} == 0,1);
%                     plm.Mp{y}{m}{c}{o}(:,idx)  = [];
%                     plm.eCm{y}{m}{c}{o}(idx,:) = [];
%                 end
%             end
%         end
%         clear y o;
        
        % Some methods don't work well if Z is empty, and there is no point in
        % using any of them all anyway.
        %if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
%         for y = loopY,
%             if opts.missingdata, loopO = 1:numel(plm.Mp{y}{m}{c}); else loopO = 1; end
%             for o = loopO,
%                 if isempty(plm.Z{y}{m}{c}{o}),
%                     plm.rmethod{y}{m}{c}{o} = 'noz';
%                 else
%                     plm.rmethod{y}{m}{c}{o} = opts.rmethod;
%                 end
%             end
%         end
%         clear y o;
        
        % Compute initial CCA (Line 19 of algorithm in Winkler et al, 2020)
        % and initialize variables
        if opts.cca.voxelwise,
            %if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
            
            % Initialize and construct residual forming matrices (lines 2-16 in Algorithm 1 in Winkler et al.
            % 2020). % ToDo: Test with rank deficient matrices
            if ~isempty(plm.Zset{1})
                plm.Z{y}{m}{c}{o} = plm.Zset{1};
                plm.Qz{y}{m}{c}{o} = zeros(plm.N,plm.N-size(plm.Z{y}{m}{c}{o},2),plm.Ysiz(1));
                % Create residual forming matrix Qz due to Z
                % HJ here is simplified as in Winkler et al, 2020 (see the Appendix text of the paper)
                for t = 1:plm.Ysiz(1),
                    [Qz,Dz,~] = svd(null(plm.Z{y}{m}{c}{o}(:,:,t)'),'econ');
                    dlmwrite(['TEST_palm_Z_' int2str(t) '.txt'],plm.Z{y}{m}{c}{o}(:,:,t))
                    plm.Qz{y}{m}{c}{o}(:,:,t) = Qz*Dz;
                    disp('Qz dimensions are:')
                    size(plm.Qz{y}{m}{c}{o})
                end
            else
                plm.Qz{y}{m}{c}{o} = zeros(plm.N,plm.N,plm.Ysiz(1));
                plm.Z{y}{m}{c}{o} = [];
                for t = 1:plm.Ysiz(1),
                    plm.Qz{y}{m}{c}{o}(:,:,t) = eye(plm.N);
                    disp('Qz dimensions are:')
                    size(plm.Qz{y}{m}{c}{o})
                end
                %plm.Qz{y}{m}{c}{o} = bsxfun(@plus,eye(plm.N),plm.Rz{y}{m}{c}{o});
            end
 
            % If right-side nuisance W not given, and this is partial CCA, reuse Z as W
            if isempty(plm.Wset{1}) && ~isfield(opts.cca,'semipartial')
                plm.W{y}{m}{c}{o} = plm.Z{y}{m}{c}{o};
            elseif ~isempty(plm.Wset{1})
                plm.W{y}{m}{c}{o} = plm.Wset{1};
            end
            
            % If right side nuisance W has been defined so far
            if isfield(plm,'W')
                % Create residual forming matrix Qw due to W
                plm.Qw{y}{m}{c}{o} = zeros(plm.N,plm.N-size(plm.W{y}{m}{c}{o},2),plm.Xsiz(1));
                size(plm.Qw{y}{m}{c}{o})
                for t = 1:plm.Xsiz(1),
                    [Qw,Dw,~] = svd(null(plm.W{y}{m}{c}{o}(:,:,t)'),'econ');
                    dlmwrite(['TEST_palm_W_' int2str(t) '.txt'],plm.W{y}{m}{c}{o}(:,:,t))
                    plm.Qw{y}{m}{c}{o}(:,:,t) = Qw*Dw;
                    disp('Qw dimensions are:')
                    size(plm.Qw{y}{m}{c}{o})
                end
            else % Qw is identity
                plm.W{y}{m}{c}{o} = [];
                plm.Qw{y}{m}{c}{o} = zeros(plm.N,plm.N,plm.Xsiz(1));
                for t = 1:plm.Xsiz(1),
                    plm.Qw{y}{m}{c}{o}(:,:,t) = eye(plm.N);
                    disp('Qw dimensions are:')
                    size(plm.Qw{y}{m}{c}{o})
                end
            end
            
            % If this is semipartial, make sure Qz is apply to the right
            % side if specified, by swapping Qw (I) and Qz and Z and W
            if isfield(opts.cca,'semipartial')
                if strcmp(opts.cca.semipartial.side,'right')
                    tmpQw  = plm.Qw; tmpW = plm.W;
                    plm.Qw = plm.Qz;
                    plm.Qz = tmpQw;
                    plm.W  = plm.Z;
                    plm.Z  = tmpW;
                    
                    disp('Qw dimensions are now:')
                    size(plm.Qw{y}{m}{c}{o})
                    
                    disp('Qz dimensions are now:')
                    size(plm.Qz{y}{m}{c}{o})
                    
                    disp('Z dimensions are now:')
                    size(plm.Z{y}{m}{c}{o})
                    
                    disp('W dimensions are now:')
                    size(plm.W{y}{m}{c}{o})
                    
                end; clear tmpQw tmpW
            end
            
           
%             for y = loopY,
%                 if y ~= m,
%                     %plm.X{y}{m}{c}{o}   = plm.X{m}{m}{c}{o};
%                     plm.Z{y}{m}{c}{o}   = plm.Z{m}{m}{c}{o};
%                     plm.W{y}{m}{c}{o}   = plm.W{m}{m}{c}{o};
%                 end
%             end
%             clear y o;
            % CCA
            %if opts.cca.do
            y = 1; o = 1;
            
            % Make the 3D dataset & residualise wrt Z
            plm.Y{m}{c} = permute(cat(3,plm.Yset{:}),[1 3 2]);
            plm.X{y}{m}{c}{o} = permute(cat(3,plm.Xset{:}),[1 3 2]);
            
            % Initialize residualized X and Y variables
            plm.Yr{m}{c}=zeros(size(plm.Y{m}{c},1)-size(plm.Z{y}{m}{c}{o},2),...
                size(plm.Y{m}{c},2),size(plm.Y{m}{c},3));
            plm.Xr{m}{c}=zeros(size(plm.X{1}{m}{c}{1},1)-size(plm.W{y}{m}{c}{o},2),...
                size(plm.X{1}{m}{c}{1},2),size(plm.X{1}{m}{c}{1},3));
            
            disp('size of Yr:')
            size(plm.Yr{m}{c})
            
            disp('size of Xr:')
            size(plm.Xr{m}{c})
            
            % Create residualized X and Y variables
            for t = 1:size(plm.Yr{m}{c},3)
                plm.Yr{m}{c}(:,:,t) = plm.Qz{1}{m}{c}{1}(:,:,t)'*plm.Y{m}{c}(:,:,t);
                dlmwrite(['TEST_palm_Y_' int2str(t) '.txt'],permute(plm.Y{m}{c}(:,:,t),[1 3 2]))
                dlmwrite(['TEST_palm_Qz_' int2str(t) '.txt'],permute(plm.Qz{1}{m}{c}{1}(:,:,t),[1 3 2]))
            end; t = 1;
            for t = 1:size(plm.Xr{m}{c},3)
                plm.Xr{m}{c}(:,:,t) = plm.Qw{1}{m}{c}{1}(:,:,t)'*plm.X{1}{m}{c}{1}(:,:,t);
                dlmwrite(['TEST_palm_X_' int2str(t) '.txt'],permute(plm.X{1}{m}{c}{1}(:,:,t),[1 3 2]))
                dlmwrite(['TEST_palm_Qw_' int2str(t) '.txt'],permute(plm.Qw{1}{m}{c}{1}(:,:,t),[1 3 2]))
            end; t = 1;
            
            dlmwrite('TEST_palm_Yr.txt',plm.Yr{m}{c},'delimiter','\t')
            dlmwrite('TEST_palm_Xr.txt',plm.Xr{m}{c},'delimiter','\t')
            
            % Compute initial CCA (Line 19 of algorithm in Winkler et al, 2020)
            % initialize some variables
            yselq   = true(1,1,plm.Ysiz(1));
            % Upper case U and V will have N' and N" rows
            plm.U{m}{c} = zeros(size(plm.Yr{m}{c}));
            plm.V{m}{c} = zeros(size(plm.Xr{m}{c}));
            % Smaller case u and v will have N rows
            plm.u{m}{c} = zeros(size(plm.Y{m}{c}));
            plm.v{m}{c} = zeros(size(plm.X{1}{m}{c}{1}));
            
            % Compute U and V
            % [A,B,r] = cca(Qz*Y,Qw*X,R,S); For now assume R==S
            
            for t = find(yselq)',
                [A,B,r] = plm.qfun(plm.Qz{1}{m}{c}{1}(:,:,t)*plm.Yr{m}{c}(:,:,t),...
                    plm.Qw{1}{m}{c}{1}(:,:,t)*plm.Xr{m}{c}(:,:,t),size(plm.Z{1}{m}{c}{1},2),...
                    size(plm.Z{1}{m}{c}{1},2));
                
                % Store A, B, U, V, u and v. A, B, U and V will later be
                % written out by default. ToDo: Add Option like saveglm to save extra
                % outputs including u and v with N rows
                plm.A{m}{c}(:,:,t) = A;
                plm.B{m}{c}(:,:,t) = B;
                plm.U{m}{c}(:,:,t) = plm.Yr{m}{c}(:,:,t)*[A null(A')];
                plm.V{m}{c}(:,:,t) = plm.Xr{m}{c}(:,:,t)*[B null(B')];
                
                % change to
                plm.u{m}{c}(:,:,t) = plm.Qz{y}{m}{c}{1}(:,:,t)*plm.U{m}{c}(:,:,t);
                plm.v{m}{c}(:,:,t) = plm.Qw{y}{m}{c}{1}(:,:,t)*plm.V{m}{c}(:,:,t);
                
                % Initialise counter and lW. Test only the requested
                % number of canonical correlations (not all of them).
                % K = numel(r);
                K = opts.ccaorplsparm;
                plm.cnt{m}{c}(t,:) = zeros(1,K);
                plm.lW{m}{c}(t,:)  = zeros(1,K);
                
                % Create masks for A, B, U and V outputs
                if t==1
                    plm.maskA = plm.maskinter; plm.maskB = plm.maskinter;
                    plm.maskU = plm.maskinter; plm.maskV = plm.maskinter;
                    
                    % For 2D (csv) inputs
                    plm.maskA.data = repmat(plm.maskinter.data,size(plm.A{m}{c}(:,:,t),1),1);
                    plm.maskB.data = repmat(plm.maskinter.data,size(plm.B{m}{c}(:,:,t),1),1);
                    plm.maskU.data = repmat(plm.maskinter.data,size(plm.U{m}{c}(:,:,t),1),1);
                    plm.maskV.data = repmat(plm.maskinter.data,size(plm.V{m}{c}(:,:,t),1),1);
                end
            end; clear t A B r K
            %end
            clear y o;
        else % i.e., spatial CCA
            
            if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
            clear y o;
            
            if opts.cca.do
                y = 1; o = 1;
                % Residual forming matrix (Z only). ToDo: Test with rank
                % deficient matrix
                plm.Qz{y}{m}{c}{o} = zeros(plm.N,plm.N-size(plm.Z{y}{m}{c}{o},2),plm.Ysiz(1));
                
                if isempty(plm.Z{y}{m}{c}{o}),
                    %plm.Rz{y}{m}{c}{o} = eye(plm.N);
                    plm.Qz{y}{m}{c}{o} = eye(plm.N);
                else
                    % HJ here is simplified as in Winkler et al, 2020 (see the Appendix text of the paper)
                    [Qz,Dz,~] = svd(null(plm.Z{y}{m}{c}{o}'),'econ');
                    dlmwrite('TEST_palm_Z.txt',plm.Z{y}{m}{c}{o})
                    dlmwrite('TEST_palm_X.txt',plm.X{y}{m}{c}{o})
                    plm.Qz{y}{m}{c}{o} = Qz*Dz;
                    % For now make Qw the same as Qz
                    plm.Qw{y}{m}{c}{o}=plm.Qz{y}{m}{c}{o};
                    plm.Rz{y}{m}{c}{o} = eye(plm.N) - plm.Z{y}{m}{c}{o}*pinv(plm.Z{y}{m}{c}{o});
                end
                
                % Make the 3D dataset & residualise wrt Z
                plm.Yr{m}{c} = cat(3,plm.Yset{:});
                
                Y_tmp = permute(plm.Yr{m}{c},[1 3 2]);
                dlmwrite('TEST_palm_Y.txt',Y_tmp,'delimiter','\t')
                
                tmp=zeros(size(plm.Yr{m}{c},1)-size(plm.Z{y}{m}{c}{o},2),...
                    size(plm.Yr{m}{c},2),size(plm.Yr{m}{c},3));
                
                for y = 1:plm.nY,
                    %plm.Yq{m}{c}(:,:,y) = plm.Rz{1}{m}{c}{o}*plm.Yq{m}{c}(:,:,y);
                    tmp(:,:,y) = plm.Qz{1}{m}{c}{o}'*plm.Yr{m}{c}(:,:,y);
                end; y = 1;
                %plm.Yq{m}{c} = permute(plm.Yq{m}{c},[1 3 2]);
                plm.Yr{m}{c} = permute(tmp,[1 3 2]);
                dlmwrite('TEST_palm_Yr.txt',plm.Yr{m}{c},'delimiter','\t')
                disp('Y after residualization and swapping dimensions')
                
                % Compute initial CCA (Line 19 of algorithm in Winkler et al, 2020)
                % residualize X into new variable Xr with N-R rows
                plm.X{y}{m}{c}{1}
                plm.Xr{y}{m}{c}{1} = plm.Qw{y}{m}{c}{1}'*plm.X{y}{m}{c}{1};
                dlmwrite('TEST_palm_Xr.txt',plm.Xr{m}{c}{1},'delimiter','\t')
                
                % initialize some variables
                yselq   = true(1,1,plm.Ysiz(1));
                plm.U{m}{c} = zeros(size(plm.Yr{m}{c}));
                plm.V{m}{c} = zeros(size(plm.Xr{y}{m}{c}{1}));
                
                % Compute U and V
                % [A,B,r] = cca(Qz*Y,Qw*X,R,S); For now assume R==S
                
                for t = find(yselq)',
                    size(plm.Qz{y}{m}{c}{1})
                    size(plm.Yr{m}{c})
                    [A,B,r] = plm.qfun(plm.Qz{y}{m}{c}{1}*plm.Yr{m}{c}(:,:,t),...
                        plm.Qw{y}{m}{c}{1}*plm.Xr{y}{m}{c}{1},size(plm.Z{y}{m}{c}{o},2),...
                        size(plm.Z{y}{m}{c}{o},2));
                    
                    % Store A, B, U, V, u and v. A, B, U and V will later be
                    % written out by default. Add Option like saveglm to save extra
                    % outputs including u and v with N rows
                    plm.A{m}{c}(:,:,t) = A;
                    plm.B{m}{c}(:,:,t) = B;
                    plm.U{m}{c}(:,:,t) = plm.Yr{m}{c}(:,:,t)*[A null(A')];
                    plm.V{m}{c}(:,:,t) = plm.Xr{y}{m}{c}{1}*[B null(B')];
                    
                    % change to
                    %plm.u{m}{c}(:,:,t) = plm.Qz{y}{m}{c}{1}(:,:,t)*plm.U{m}{c}(:,:,t);
                    %plm.v{m}{c}(:,:,t) = plm.Qw{y}{m}{c}{1}(:,:,t)*plm.V{m}{c}(:,:,t);
                    
                    % Initialise counter and lW. Test only the requested
                    % number of canonical correlations (not all of them).
                    % K = numel(r);
                    K = opts.ccaorplsparm;
                    plm.cnt{m}{c}(t,:) = zeros(1,K);
                    plm.lW{m}{c}(t,:)  = zeros(1,K);
                    
                    % Create masks for A, B, U and V outputs
                    if t==1
                        plm.maskA = plm.maskinter; plm.maskB = plm.maskinter;
                        plm.maskU = plm.maskinter; plm.maskV = plm.maskinter;
                        
                        % For 2D (csv) inputs
                        plm.maskA.data = repmat(plm.maskinter.data,size(plm.A{m}{c}(:,:,t),1),1);
                        plm.maskB.data = repmat(plm.maskinter.data,size(plm.B{m}{c}(:,:,t),1),1);
                        plm.maskU.data = repmat(plm.maskinter.data,size(plm.U{m}{c}(:,:,t),1),1);
                        plm.maskV.data = repmat(plm.maskinter.data,size(plm.V{m}{c}(:,:,t),1),1);
                    end
                    
                end; clear t A B r
            end
        end
    end
end

% Create the permutation set, while taking care of the synchronized
% permutations (see the inner loop below)
if opts.syncperms,
    if ~ opts.accel.noperm,
        if isempty(plm.EB),
            if opts.savemetrics,
                [plm.Pset,plm.nP{1}(1),plm.metr{1}{1}] = ...
                    palm_shuffree(plm.seq{1}{1},opts.nP0, ...
                    opts.cmcp,opts.EE,opts.ISE,opts.idxout);
            else
                [plm.Pset,plm.nP{1}(1)] = ...
                    palm_shuffree(plm.seq{1}{1},opts.nP0, ...
                    opts.cmcp,opts.EE,opts.ISE,opts.idxout);
            end
        else
            if opts.savemetrics,
                [plm.Pset,plm.nP{1}(1),plm.metr{1}{1}] = ...
                    palm_shuftree(opts,plm,1,1);
            else
                [plm.Pset,plm.nP{1}(1)] = ...
                    palm_shuftree(opts,plm,1,1);
            end
        end
        fprintf('Building null distribution.\n');
    else
        fprintf('Doing the approximation without permutations.\n');
    end
    
    % This is for the negative binomial mode
    if opts.accel.negbin && ~ opts.saveunivariate,
        dothisY = false(plm.nY,1);
    else
        dothisY = true(plm.nY,1);
    end
    dotheMVorCCAorPLS = true;
    if opts.accel.noperm,
        P_outer = 1;
    else
        P_outer = 1:plm.nP{1}(1);
    end
else
    P_outer = 1;
end

% To calculate progress
if opts.syncperms,
    ProgressNum = 0;
    if opts.designperinput,
        ProgressDen = sum(plm.nC) * plm.nP{1}(1);
    else
        ProgressDen = sum(plm.nC) * plm.nP{1}(1) * plm.nY;
    end
else
    ProgressCon = 0;
end

% For each permutation (outer loop)
ticP = tic;
for po = P_outer,
    
    % For each design matrix
    for m = 1:plm.nM,
        
        % String with the counter
        if po == 1 && (plm.nM > 1 || opts.verbosefilenames),
            plm.mstr{m} = sprintf('_d%d',m);
        end
        
        % For each contrast
        for c = 1:plm.nC(m),
            
            % String with the counter
            if max(plm.nC) > 1 || opts.verbosefilenames,
                ctmp = c + opts.conskipcount;
                plm.cstr{m}{c} = sprintf('_c%d',ctmp);
            end
            
            % This is for the negative binomial mode
            if ~ opts.syncperms,
                if opts.accel.negbin && ~ opts.saveunivariate,
                    dothisY = false(plm.nY,1);
                else
                    dothisY = true(plm.nY,1);
                end
                dotheMVorCCAorPLS = true;
            end
            
            % Create the permutation set, while taking care of the synchronized
            % permutations (see the outer loop above)
            if opts.syncperms,
                P_inner = po;
                plm.nP{m}(c) = plm.nP{1}(1);
            else
                if ~ opts.accel.noperm,
                    if isempty(plm.EB),
                        if opts.savemetrics,
                            [plm.Pset,plm.nP{m}(c),plm.metr{m}{c}] = ...
                                palm_shuffree(plm.seq{m}{c},opts.nP0, ...
                                opts.cmcp,opts.EE,opts.ISE,opts.idxout);
                        else
                            [plm.Pset,plm.nP{m}(c)] = ...
                                palm_shuffree(plm.seq{m}{c},opts.nP0, ...
                                opts.cmcp,opts.EE,opts.ISE,opts.idxout);
                        end
                    else
                        if opts.savemetrics,
                            [plm.Pset,plm.nP{m}(c),plm.metr{m}{c}] = ...
                                palm_shuftree(opts,plm,m,c);
                        else
                            [plm.Pset,plm.nP{m}(c)] = ...
                                palm_shuftree(opts,plm,m,c);
                        end
                    end
                    fprintf('Building null distribution.\n');
                else
                    fprintf('Doing the approximation without permutations.\n');
                end
                P_inner = 1:plm.nP{m}(c);
            end
            
            if po == 1 && ~ opts.accel.noperm,
                % If the user wants to save the permutations, save the vectors now.
                % This has 3 benefits: (1) the if-test below will run just once, rather
                % than many times inside the loop, (2) if the user only wants the
                % vectors, not the images, he/she can cancel immediately after the
                % text file has been created and (3) having all just as a single big
                % file is more convenient than hundreds of small ones.
                if opts.saveperms,
                    % It's faster to write directly as below than using dlmwrite and
                    % palm_swapfmt.m
                    fid = fopen(sprintf('%s%s%s_permidx.csv',opts.o,plm.mstr{m},plm.cstr{m}{c}),'w');
                    for p = 1:plm.nP{m}(c),
                        if iscell(plm.Pset)
                            fprintf(fid,'%d,',palm_perm2idx(plm.Pset{p})');
                        else
                            fprintf(fid,'%d,',plm.Pset(:,p));
                        end
                        fseek(fid,-1,'eof');
                        fprintf(fid,'\n');
                    end
                    fclose(fid);
                end
                
                % If the user requests, save the permutation metrics
                if opts.savemetrics,
                    fid = fopen(sprintf('%s%s%s_metrics.csv',opts.o,plm.mstr{m},plm.cstr{m}{c}),'w');
                    fprintf(fid,[ ...
                        'Log of max number of permutations given the tree (W),%f\n' ...
                        'Log of max number of permutations if unrestricted (W0),%f\n' ...
                        'Huberman & Hogg complexity (tree only),%d\n' ...
                        'Huberman & Hogg complexity (tree & design),%d\n' ...
                        'Average Hamming distance (tree only),%f\n' ...
                        'Average Hamming distance (tree & design),%f\n' ...
                        'Average Euclidean distance (tree only),%f\n' ...
                        'Average Euclidean distance (tree & design),%f\n' ...
                        'Average Spearman correlation,%f\n'], plm.metr{m}{c});
                    fclose(fid);
                end
            end
            
            if ~ opts.syncperms,
                ProgressNum = 0;
            end
            
            % For each permutation (inner loop):
            for p = P_inner,
                
                % For each input dataset
                if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
                
                if opts.cca.do % ToDO: Can get rid of opt.PLS
                    
                    % This if is for the negative binomial mode.
                    if dotheMVorCCAorPLS,
                        if opts.showprogress,
                            if opts.syncperms,
                                fprintf('\t [Shuffling %d/%d, %s]\n', ...
                                    p,plm.nP{m}(c),m,upper(plm.Qname{m}{c}(2:4)));
                            else
                                fprintf('\t [Shuffling %d/%d, %s]\n', ...
                                    p,plm.nP{m}(c),upper(plm.Qname{m}{c}(2:4)));
                            end
                        end
                        
                        % Initialise some vars
                        if p == 1,
                            yselq   = true(1,1,plm.Ysiz(1));
                        end
                        
                        % Compute the CC coefficient
                        y = 1; m = 1; c = 1;
                        % Set S and R vars (ToDo: set these earlier) 
                        R=size(plm.Z{y}{m}{c}{1},2); 
                        S=size(plm.W{y}{m}{c}{1},2);
                        
                        if opts.cca.voxelwise, % Compute voxelwise CCA
                            for t = find(yselq)',
                                for k=1:length(plm.cnt{m}{c}(t,:)),
                                    % Shuffle either u or v
                                    if isfield(opts.cca,'semipartial') && strcmp(opts.cca.semipartial.side,'right')
                                        u=plm.Qz{y}{m}{c}{1}(:,:,t)*plm.U{m}{c}(:,k:end,t);
                                        v=plm.Qw{y}{m}{c}{1}(:,:,t)*plm.V{m}{c}(plm.Pset(:,p),k:end,t);
                                    else
                                        u=plm.Qz{y}{m}{c}{1}(:,:,t)*plm.U{m}{c}(plm.Pset(:,p),k:end,t);
                                        v=plm.Qw{y}{m}{c}{1}(:,:,t)*plm.V{m}{c}(:,k:end,t);
                                    end
                                    [~,~,r] = plm.qfun(u,v,R,S);
                                    lWtmp = -fliplr(cumsum(fliplr(log(1-r.^2))));
                                    lW(k) = lWtmp(1);
                                    Rtmp(k)=r(1);
                                end
                                if p == 1
                                    % save Wilk's statistic
                                    plm.lW{m}{c}(t,:) = lW;
                                    %Rtmp2(t,:)=Rtmp;
                                    %Q{m}{c}(t) = r(1);
                                end
                                % Store canonical correlations to Q
                                Q{m}{c}(t,:)=Rtmp;
                                plm.cnt{m}{c}(t,:) = plm.cnt{m}{c}(t,:) + (lW >= plm.lW{m}{c}(t,:));
                                
                            end; clear t r u v k lW
                            dlmwrite('TEST_palm_lW.txt',plm.lW{m}{c},'delimiter','\t')
                            dlmwrite('TEST_palm_r.txt',Q{m}{c},'delimiter','\t')
                            dlmwrite('TEST_palm_A.txt',plm.A{m}{c},'delimiter','\t')
                            
                            
                        else % Compute Spatial CCA
                            % Main loop cca (Line 27-32 of algorithm from
                            % Winkler et al 2020.
                            
                            for t = find(yselq)',
                                %Q{m}{c}(t) = plm.qfun(plm.Yq{m}{c}(:,:,t),M,opts.ccaorplsparm);
                                for k=1:length(plm.cnt{m}{c}(t,:)),
                                    % Shuffle just u
                                    u=plm.Qz{y}{m}{c}{1}*plm.U{m}{c}(plm.Pset(:,p),k:end,t);
                                    v=plm.Qw{y}{m}{c}{1}*plm.V{m}{c}(:,k:end,t);
                                    [~,~,r] = plm.qfun(u,v,R,S);
                                    lWtmp = -fliplr(cumsum(fliplr(log(1-r.^2))));
                                    lW(k) = lWtmp(1);
                                    Rtmp(k)=r(1);
                                end
                                if p == 1
                                    % save Wilk's test statistic for
                                    % unshuffled data
                                    plm.lW{m}{c}(t,:) = lW;
                                    %Rtmp2(t,:)=Rtmp;
                                    %Q{m}{c}(t) = r(1);
                                end
                                % store the canonical correlations into Q
                                Q{m}{c}(t,:)=Rtmp;
                                plm.cnt{m}{c}(t,:) = plm.cnt{m}{c}(t,:) + (lW >= plm.lW{m}{c}(t,:));
                                
                            end; clear t r u v k lW
                            dlmwrite('TEST_palm_lW.txt',plm.lW{m}{c},'delimiter','\t')
                            dlmwrite('TEST_palm_r.txt',Q{m}{c},'delimiter','\t')
                        end
                        
                        % Convert to zstat if that was asked
                        if opts.zstat,
                            Q{m}{c}(yselq) = atanh(Q{m}{c}(yselq));
                            if p == 1 && m == 1 && c == 1,
                                plm.mvstr = horzcat('_z',plm.mvstr(2:end));
                            end
                        end
                        
                        % Save the canonical correlations
                        if p == 1,
                            for nc=1:opts.ccaorplsparm
                                if opts.ccaorplsparm > 99
                                    ccname = [plm.Qname{m}{c} num2str(nc,'%03d')];
                                else
                                    ccname = [plm.Qname{m}{c} num2str(nc,'%02d')];
                                end
                                palm_quicksave(Q{m}{c}(:,nc),0,opts,plm,[],m,c, ...
                                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_r',plm.mstr{m},plm.cstr{m}{c}));
                            end; clear ccname nc
                        end
                        
                        % Save out permutations if requested
                        if opts.saveperms,
                            for nc=1:opts.ccaorplsparm
                                ccname = [plm.Qname{m}{c} int2str(nc)];
                                palm_quicksave(Q{m}{c}(:,nc),0,opts,plm,[],m,c, ...
                                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_r',plm.mstr{m},plm.cstr{m}{c},sprintf('_perm%06d',p)));
                            end; clear ccname
                            %palm_quicksave(Q{m}{c},0,opts,plm,[],m,c, ...
                            %    horzcat(sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{m}{c}),sprintf('_perm%06d',p)));
                        end
                        
                        % In the first permutation, keep Q and start the counter.
                        if p == 1,
                            plm.Q     {m}{c} = Q{m}{c};
                            plm.Qpperm{m}{c} = zeros(size(Q{m}{c}));
                        end
                        plm.Qpperm{m}{c}     = plm.Qpperm{m}{c} + (Q{m}{c} >= plm.Q{m}{c});
                        
                    end
                    if opts.accel.negbin && ~ any(yselq),
                        dotheMVorCCAorPLS = false;
                    end
                end
            end
        end
        
        % Output uncorrected p-values for now and check against permcca
        if opts.cca.do
            dlmwrite('TEST_palm_punc.txt',plm.cnt{m}{c}./p,'delimiter','\t')
        end
        
        if ~ opts.syncperms,
            ProgressCon = ProgressCon + 1;
        end
    end
end; clear po

tocP = toc(ticP);
fprintf('Elapsed time with permutations: ~ %g seconds.\n',tocP);
clear('M','Y','psi','res','G','df2','T','Q');

% Save everything, except the few bits saved above
% ticS = tic;
% palm_saveall(plm,opts);
% tocS = toc(ticS);
% fprintf('Elapsed time generating and saving results: ~ %g seconds.\n',tocS);
% fprintf('Overall elapsed time: ~ %g seconds.\n',tocI+tocP+tocS);
% csvwrite(sprintf('%s_elapsed.csv',opts.o),[tocI tocP tocS]);

% Finished.
fprintf('PALM finished at %s.\n',datestr(now));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  F U N C T I O N S  %%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==============================================================
% Below are the functions for each of the regression methods:
% ==============================================================
function [Mr,Y] = noz(P,Y,y,m,c,o,plm)
% This is equivalent to Draper-Stoneman, as when there is no Z
% Y remains unchanged.
Mr = P*plm.X{y}{m}{c}{o};
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = noz3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.X{y}{m}{c}{o}));
for t = 1:size(Y,2),
    Mr(:,:,t) = P*plm.X{y}{m}{c}{o}(:,:,t);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = nozm(P,Y,y,m,c,o,plm,ikeep)
Y  = Y(ikeep,:);
Mr = P*plm.X{y}{m}{c}{o};

% ==============================================================
function [Mr,Yr] = exact(P,Y,y,m,c,o,plm)
% The "exact" method, in which the coefficients for
% the nuisance are known.
Yr = Y - plm.Z{y}{m}{c}{o}*plm.g;
Mr = P*plm.X{y}{m}{c}{o};
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = exact3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.X{y}{m}{c}{o}));
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t)   = Y(:,t) - plm.Z{y}{m}{c}{o}(:,:,t)*plm.g;
    Mr(:,:,t) = P*plm.X{y}{m}{c}{o}(:,:,t);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = exactm(P,Y,y,m,c,o,plm,ikeep)
% The "exact" method, in which the coefficients for
% the nuisance are known.
Yr = Y(ikeep,:) - plm.Z{y}{m}{c}{o}(ikeep,:)*plm.g;
Mr = P*plm.X{y}{m}{c}{o};

% ==============================================================
function [Mr,Y] = draperstoneman(P,Y,y,m,c,o,plm)
% Draper and Stoneman (1966) method.
% Y remains unchanged
Mr = horzcat(P*plm.X{y}{m}{c}{o},plm.Z{y}{m}{c}{o});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = draperstoneman3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.Mp{y}{m}{c}{o}));
for t = 1:size(Y,2),
    Mr(:,:,t) = horzcat(P*plm.X{y}{m}{c}{o}(:,:,t),plm.Z{y}{m}{c}{o}(:,:,t));
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = draperstonemanm(P,Y,y,m,c,o,plm,ikeep)
Mr = horzcat(P*plm.X{y}{m}{c}{o},plm.Z{y}{m}{c}{o}(ikeep,:));
Yr = Y(ikeep,:);

% ==============================================================
function [Mr,Yr] = stillwhite(P,Y,y,m,c,o,plm)
% A method following the same logic as the one
% proposed by Still and White (1981)
Yr = plm.Rz{y}{m}{c}{o}*Y;
Mr = P*plm.X{y}{m}{c}{o};
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = stillwhite3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.X{y}{m}{c}{o}));
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t)   = plm.Rz{y}{m}{c}{o}*Y(:,t);
    Mr(:,:,t) = P*plm.X{y}{m}{c}{o}(:,:,t);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = stillwhitem(P,Y,y,m,c,o,plm,ikeep)
Yr = plm.Rz{y}{m}{c}{o}(ikeep,:)*Y;
Mr = P*plm.X{y}{m}{c}{o};

% ==============================================================
function [Mr,Yr] = freedmanlane(P,Y,y,m,c,o,plm)
% The Freedman and Lane (1983) method.
Mr = plm.Mp{y}{m}{c}{o};
Yr = (P'*plm.Rz{y}{m}{c}{o} + plm.Hz{y}{m}{c}{o})*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = freedmanlane3d(P,Y,y,m,c,o,plm)
Mr = plm.Mp{y}{m}{c}{o};
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t) = (P'*plm.Rz{y}{m}{c}{o}(:,:,t) + plm.Hz{y}{m}{c}{o}(:,:,t))*Y(:,t);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = freedmanlanem(P,Y,y,m,c,o,plm,ikeep)
Mr = plm.Mp{y}{m}{c}{o}(ikeep,:);
Yr = (P*plm.Rz{y}{m}{c}{o} + plm.Hz{y}{m}{c}{o}(ikeep,:))*Y;

% ==============================================================
function [Mr,Yr] = manly(P,Y,y,m,c,o,plm)
% The Manly (1986) method.
% There's no need for a 3D version of this method.
Mr = plm.Mp{y}{m}{c}{o};
Yr = P'*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = manlym(P,Y,y,m,c,o,plm,ikeep)
Mr = plm.Mp{y}{m}{c}{o}(ikeep,:);
Yr = P*Y;

% ==============================================================
function [Mr,Yr] = terbraak(P,Y,y,m,c,o,plm)
% The ter Braak (1992) method.
Mr = plm.Mp{y}{m}{c}{o};
Yr = (P'*plm.Rm{y}{m}{c}{o} + plm.Hm{y}{m}{c}{o})*Y; % original method
% Yr = P'*plm.Rm{y}{m}{c}{o}*Y; % alternative (causes unpermuted stat to be 0)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = terbraak3d(P,Y,y,m,c,o,plm)
Mr = plm.Mp{y}{m}{c}{o};
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t) = (P'*plm.Rm{y}{m}{c}{o}(:,:,t) + plm.Hm{y}{m}{c}{o}(:,:,t))*Y(:,t); % original method
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = terbraakm(P,Y,y,m,c,o,plm,ikeep)
Mr = plm.Mp{y}{m}{c}{o}(ikeep,:);
Yr = (P'*plm.Rm{y}{m}{c}{o} + plm.Hm{y}{m}{c}{o}(ikeep,:))*Y; % original method

% ==============================================================
function [Mr,Yr] = kennedy(P,Y,y,m,c,o,plm)
% The Kennedy (1996) method. This method should NEVER be used.
Mr = plm.Rz{y}{m}{c}{o}*plm.X{y}{m}{c}{o};
Yr = P'*plm.Rz{y}{m}{c}{o}*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = kennedy3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.X{y}{m}{c}{o}));
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Mr(:,:,t) = plm.Rz{y}{m}{c}{o}(:,:,t)*plm.X{y}{m}{c}{o}(:,:,t);
    Yr(:,t) = P'*plm.Rz{y}{m}{c}{o}(:,:,t)*Y(:,t);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = kennedym(P,Y,y,m,c,o,plm,ikeep)
Mr = plm.Rz{y}{m}{c}{o}(ikeep,:)*plm.X{y}{m}{c}{o};
Yr = P*plm.Rz{y}{m}{c}{o}*Y;

% ==============================================================
function [Mr,Yr] = huhjhun(P,Y,y,m,c,o,plm)
% The Huh and Jhun (2001) method, that fixes the issues
% with Kennedy's, but doesn't allow block permutation.

size(P')
size(plm.hj{y}{m}{c}{o}')
size(plm.Rz{y}{m}{c}{o})
size(Y)
pause(5)

Mr = plm.hj{y}{m}{c}{o}'*plm.Rz{y}{m}{c}{o}*plm.X{y}{m}{c}{o};
Yr = P'*plm.hj{y}{m}{c}{o}'*plm.Rz{y}{m}{c}{o}*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = huhjhun3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.hj{y}{m}{c}{o},1),size(plm.X{y}{m}{c}{o},2));
Yr = zeros(size(plm.hj{y}{m}{c}{o},1),size(Y,2));
for t = 1:size(Y,2),
    Mr(:,:,t) = plm.hj{y}{m}{c}{o}(:,:,t)'*plm.Rz{y}{m}{c}{o}(:,:,t)*plm.X{y}{m}{c}{o}(:,:,t);
    Yr(:,t)   = P'*plm.hj{y}{m}{c}{o}(:,:,t)'*plm.Rz{y}{m}{c}{o}(:,:,t)*Y(:,t);
end

% ==============================================================
function [Mr,Y] = dekker(P,Y,y,m,c,o,plm)
% The Dekker method, i.e., orthogonalization.
% Y remains unchanged
Mr = horzcat(P*plm.Rz{y}{m}{c}{o}*plm.X{y}{m}{c}{o},plm.Z{y}{m}{c}{o});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = dekker3d(P,Y,y,m,c,o,plm)
Mr = zeros(size(plm.Mp{y}{m}{c}{o}));
for t = 1:size(Y,2),
    Mr(:,:,t) = horzcat(P*plm.Rz{y}{m}{c}{o}(:,:,t)*plm.X{y}{m}{c}{o}(:,:,t),plm.Z{y}{m}{c}{o}(:,:,t));
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = dekkerm(P,Y,y,m,c,o,plm,ikeep)
Mr = horzcat(P*plm.Rz{y}{m}{c}{o}*plm.X{y}{m}{c}{o},plm.Z{y}{m}{c}{o}(ikeep,:));

% ==============================================================
% Below are the functions to compute univariate statistics:
% ==============================================================
% Reference:
% * Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE.
%   Permutation inference for the general linear model.
%   NeuroImage, 2014;92:381-397 (Open Access)
% ==============================================================
function G = fastr(M,psi,Y,y,m,c,o,plm)
% This only works if:
% - M and Y have zero mean.
% - rank(contrast) = 1
%
% Inputs:
% M   : design matrix (demeaned)
% psi : regression coefficients
% Y   : data (demeaned)
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : Pearson's correlation coefficient (r).
G = fastrsq(M,psi,Y,y,m,c,o,plm);
G = sign(plm.eC{y}{m}{c}{o}'*psi).*G.^.5;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function G = fastr3d(M,psi,Y,y,m,c,o,plm)
G = fastrsq3d(M,psi,Y,y,m,c,o,plm);
G = sign(plm.eC{y}{m}{c}{o}'*psi).*G.^.5;

% ==============================================================
function G = fastrsq(M,psi,Y,y,m,c,o,plm)
% This only works if:
% - M and Y have zero mean.
%
% Inputs:
% M   : design matrix (demeaned)
% psi : regression coefficients
% Y   : data (demeaned)
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : R^2, i.e., the coefficient of determination.
tmp = plm.mrdiv(plm.eC{y}{m}{c}{o},...
    plm.mrdiv(plm.eC{y}{m}{c}{o}',(M'*M))*plm.eC{y}{m}{c}{o})...
    *plm.eC{y}{m}{c}{o}';
G   = sum((tmp'*psi).*psi,1);
den = sum(Y.^2,1);
G   = G./den;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function G = fastrsq3d(M,psi,Y,y,m,c,o,plm)
for t = 1:size(psi,2),
    tmp = plm.mrdiv(plm.eC{y}{m}{c}{o},...
        plm.mrdiv(plm.eC{y}{m}{c}{o}',...
        (M(:,:,t)'*M(:,:,t)))*plm.eC{y}{m}{c}{o})*...
        plm.eC{y}{m}{c}{o}';
    G   = sum((tmp'*psi(:,t)).*psi(:,t),1);
end
den = sum(Y.^2,1);
G   = G./den;

% ==============================================================
function [G,df2] = fastt(M,psi,res,y,m,c,o,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups = 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : t statistic.
% df2 : Degrees of freedom. df1 is 1 for the t statistic.
df2 = size(M,1)-plm.rM{y}{m}{c}{o};
G   = plm.eC{y}{m}{c}{o}'*psi;
den = sqrt(plm.mrdiv(plm.eC{y}{m}{c}{o}',(M'*M))*plm.eC{y}{m}{c}{o}*sum(res.^2)./df2);
G   = G./den;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastt3d(M,psi,res,y,m,c,o,plm)
df2 = size(M,1)-plm.rM{y}{m}{c}{o};
G   = plm.eC{y}{m}{c}{o}'*psi;
S   = zeros(1,size(psi,2));
for t = 1:size(psi,2),
    S(t) = plm.mrdiv(plm.eC{y}{m}{c}{o}',(M(:,:,t)'*M(:,:,t)))*plm.eC{y}{m}{c}{o};
end
den = sqrt(S.*sum(res.^2,1)./df2);
G   = G./den;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fasttswe(M,psi,res,y,m,c,o,plm)
df2 = size(M,1)-plm.rM{y}{m}{c}{o};
G   = plm.eC{y}{m}{c}{o}'*psi;
S   = zeros(1,size(psi,2));
MtM = (M'*M);
Rm  = eye(plm.N) - M*pinv(M);
for t = 1:size(psi,2),
    V = res(:,t)*res(:,t)' ./ Rm;
    S(t) = plm.mrdiv(plm.eC{y}{m}{c}{o}',...
        plm.mrdiv(plm.mldiv(MtM,(M'*V*M)),MtM))*plm.eC{y}{m}{c}{o};
end
den = sqrt(S.*sum(res.^2,1)./df2);
G   = G./den;

% ==============================================================
function [G,df2] = fastf(M,psi,res,y,m,c,o,plm)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups = 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : F-statistic.
% df2 : Degrees of freedom 2. df1 is rank(C).
df2 = size(M,1)-plm.rM{y}{m}{c}{o};
cte = plm.mrdiv(plm.eC{y}{m}{c}{o},...
    plm.mrdiv(plm.eC{y}{m}{c}{o}',(M'*M))*plm.eC{y}{m}{c}{o})* ...
    plm.eC{y}{m}{c}{o}';
tmp = zeros(size(psi));
for j = 1:size(cte,2),
    tmp(j,:) = sum(bsxfun(@times,psi,cte(:,j)),1)';
end
G   = sum(tmp.*psi,1);
ete = sum(res.^2,1);
G   = G./ete*df2/plm.rC0{m}(c);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastf3d(M,psi,res,y,m,c,o,plm)
df2 = size(M,1)-plm.rM{y}{m}{c}{o};
nT = size(res,2);
cte = zeros(size(psi,1),size(psi,1),nT);
for t = 1:nT,
    cte(:,:,t) = plm.mrdiv(plm.eC{y}{m}{c}{o},plm.mrdiv(plm.eC{y}{m}{c}{o}', ...
        (M(:,:,t)'*M(:,:,t)))*plm.eC{y}{m}{c}{o})*plm.eC{y}{m}{c}{o}';
end
ppsi = permute(psi,[1 3 2]);
ppsi = sum(bsxfun(@times,ppsi,cte),1);
ppsi = permute(ppsi,[2 3 1]);
G    = sum(ppsi.*psi,1);
ete  = sum(res.^2,1);
G    = G./ete*df2/plm.rC0{m}(c);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastfswe(M,psi,res,y,m,c,o,plm)
df2 = size(M,1)-plm.rM{y}{m}{c}{o};
nT = size(res,2);
cte = zeros(size(psi,1),size(psi,1),nT);
MtM = (M'*M);
Rm = eye(plm.N) - M*pinv(M);
for t = 1:nT,
    V = res(:,t)*res(:,t)' ./ Rm;
    cte(:,:,t) = plm.mrdiv(plm.eC{y}{m}{c}{o},plm.mrdiv(plm.eC{y}{m}{c}{o}', ...
        plm.mrdiv(plm.mldiv(MtM,(M'*V*M)),MtM))* ...
        plm.eC{y}{m}{c}{o})*plm.eC{y}{m}{c}{o}';
end
ppsi = permute(psi,[1 3 2]);
ppsi = sum(bsxfun(@times,ppsi,cte),1);
ppsi = permute(ppsi,[2 3 1]);
G    = sum(ppsi.*psi,1);
ete  = sum(res.^2,1);
G    = G./ete*df2/plm.rC0{m}(c);

% ==============================================================
function [G,df2] = fastv(M,psi,res,y,m,c,o,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups > 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : Aspin-Welch v statistic.
% df2 : Degrees of freedom 2. df1 is 1.
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
den  = zeros(1,nT);
r    = size(M,2);
dRmb = zeros(plm.nVG,1);
cte  = zeros(r^2,nT);
for b = 1:plm.nVG,
    bidx    = plm.VG == b;
    dRmb(b) = sum(plm.dRm{y}{m}{c}{o}(bidx),1);
    W(b,:)  = dRmb(b)./sum(res(bidx,:).^2,1);
    Mb      = M(bidx,:)'*M(bidx,:);
    cte     = cte + Mb(:)*W(b,:);
    W(b,:)  = W(b,:)*sum(bidx);
end
for t = 1:nT,
    den(t) = plm.mrdiv(plm.eC{y}{m}{c}{o}',reshape(cte(:,t),[r r]))*plm.eC{y}{m}{c}{o};
end
G    = plm.eC{y}{m}{c}{o}'*psi./sqrt(den);
sW1  = sum(W,1);
bsum = sum(bsxfun(@rdivide,(1-bsxfun(@rdivide,W,sW1)).^2,dRmb),1);
df2  = 1/3./bsum;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastv3d(M,psi,res,y,m,c,o,plm)
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
den  = zeros(1,nT);
r    = size(M,2);
dRmb = zeros(plm.nVG,nT);
cte  = zeros(r,r,nT);
for b = 1:plm.nVG,
    bidx = plm.VG == b;
    dRmb(b,:) = sum(plm.dRm{y}{m}{c}{o}(bidx,:),1);
    W(b,:)    = dRmb(b,:)./sum(res(bidx,:).^2,1);
    for t = 1:nT,
        cte(:,:,t) = cte(:,:,t) + (M(bidx,:,t)'*M(bidx,:,t)).*W(b,t);
    end
    W(b,:) = W(b,:)*sum(bidx);
end
for t = 1:nT,
    den(t) = plm.mrdiv(plm.eC{y}{m}{c}{o}',cte(:,:,t))*plm.eC{y}{m}{c}{o};
end
G    = plm.eC{y}{m}{c}{o}'*psi./sqrt(den);
sW1  = sum(W,1);
bsum = sum((1-bsxfun(@rdivide,W,sW1)).^2./dRmb,1);
df2  = 1/3./bsum;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastvswe(M,psi,res,y,m,c,o,plm)
nT   = size(res,2);
den  = zeros(1,nT);
r    = size(M,2);
cte  = zeros(r,r,nT);
for b = 1:plm.nVG,
    bidx = plm.VG == b;
    MtM  = (M(bidx,:)'*M(bidx,:));
    Rm   = eye(sum(bidx)) - M(bidx,:)*pinv(M(bidx,:));
    for t = 1:nT,
        V = res(bidx,t)*res(bidx,t)' ./ Rm;
        cte(:,:,t) = cte(:,:,t) + plm.mrdiv(plm.mldiv(MtM,(M(bidx,:)'*V*M(bidx,:))),MtM);
    end
end
for t = 1:nT,
    den(t) = plm.mrdiv(plm.eC{y}{m}{c}{o}',cte(:,:,t))*plm.eC{y}{m}{c}{o};
end
G    = plm.eC{y}{m}{c}{o}'*psi./sqrt(den);
df2  = size(M,1)-plm.rM{y}{m}{c}{o};

% ==============================================================
function [G,df2] = fastg(M,psi,res,y,m,c,o,plm)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups > 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : Welch v^2 statistic.
% df2 : Degrees of freedom 2. df1 is rank(C).
r    = size(M,2);
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
dRmb = zeros(plm.nVG,1);
cte  = zeros(r^2,nT);
for b = 1:plm.nVG,
    bidx    = plm.VG == b;
    dRmb(b) = sum(plm.dRm{y}{m}{c}{o}(bidx));
    W(b,:)  = dRmb(b)./sum(res(bidx,:).^2);
    Mb      = M(bidx,:)'*M(bidx,:);
    cte     = cte + Mb(:)*W(b,:);
    W(b,:)  = W(b,:)*sum(bidx);
end
G = zeros(1,nT);
for t = 1:nT,
    A = psi(:,t)'*plm.eC{y}{m}{c}{o};
    G(t) = plm.mrdiv(A,plm.mrdiv(plm.eC{y}{m}{c}{o}',reshape(cte(:,t),[r r]))* ...
        plm.eC{y}{m}{c}{o})*A'/plm.rC0{m}(c);
end
sW1  = sum(W,1);
bsum = sum(bsxfun(@rdivide,(1-bsxfun(@rdivide,W,sW1)).^2,dRmb),1);
bsum = bsum/plm.rC0{m}(c)/(plm.rC0{m}(c)+2);
df2  = 1/3./bsum;
G    = G./(1 + 2*(plm.rC0{m}(c)-1).*bsum);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastg3d(M,psi,res,y,m,c,o,plm)
r    = size(M,2);
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
dRmb = zeros(plm.nVG,nT);
cte  = zeros(r,r,nT);
for b = 1:plm.nVG,
    bidx = plm.VG == b;
    dRmb(b,:) = sum(plm.dRm{y}{m}{c}{o}(bidx,:),1);
    W(b,:) = dRmb(b,:)./sum(res(bidx,:).^2,1);
    for t = 1:nT,
        cte(:,:,t) = cte(:,:,t) + (M(bidx,:,t)'*M(bidx,:,t))*W(b,t);
    end
    W(b,:) = W(b,:)*sum(bidx);
end
G = zeros(1,nT);
for t = 1:nT,
    A = psi(:,t)'*plm.eC{y}{m}{c}{o};
    G(t) = plm.mrdiv(A,plm.mrdiv(plm.eC{y}{m}{c}{o}',cte(:,:,t)) * ...
        plm.eC{y}{m}{c}{o})*A'/plm.rC0{m}(c);
end
sW1  = sum(W,1);
bsum = sum((1-bsxfun(@rdivide,W,sW1)).^2./dRmb,1);
bsum = bsum/plm.rC0{m}(c)/(plm.rC0{m}(c)+2);
df2  = 1/3./bsum;
G    = G./(1 + 2*(plm.rC0{m}(c)-1).*bsum);

% ==============================================================
function [Z,df2] = fastmiss(Y,M,y,m,c,plm,opts,fastpiv)
% Conputes the test statistic for missing data.
df2 = NaN;
persistent GPtmp df2tmp; % persistent so as to avoid re-allocing.
GPtmp  = zeros(numel(Y),size(Y{1},2)); % same var for G and P
df2tmp = GPtmp;
nO = numel(Y);
for o = nO:-1:1,
    if plm.isdiscrete{y}{m}{c}(o),
        GPtmp(o,:) = yates(Y{o},M{o});
        GPtmp(o,:) = palm_gpval(GPtmp(o,:),0);
    elseif testzeros(Y{o},M{o},y,m,c,o,plm),
        GPtmp(o,:) = []; df2tmp(o,:) = [];
    else
        psi = plm.mldiv(M{o},Y{o});
        res = Y{o} - M{o}*psi;
        [GPtmp(o,:),df2tmp(o,:)] = fastpiv(M{o},psi,res,y,m,c,o,plm);
        if plm.istwotail{y}{m}{c}(o) || opts.twotail,
            GPtmp(o,:) = abs(GPtmp(o,:));
        end
        GPtmp(o,:) = palm_gpval(GPtmp(o,:),plm.rC0{m}(c),df2tmp(o,:));
    end
end
G = -2*sum(log(GPtmp),1);
P = palm_gpval(G,-1,2*size(GPtmp,1));
Z = sqrt(2)*erfcinv(2*P);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function result = testzeros(Y,M,y,m,c,o,plm)
Mhaszero = any(all(M(:,any(plm.eCm{y}{m}{c}{o},2)) == 0,1),2);
Yhaszero = all(Y(:) == 0);
result = Yhaszero | Mhaszero;

% ==============================================================
% Below are the functions to compute multivariate statistics:
% ==============================================================
% Reference:
% * Winkler AM, Webster MA, Brooks JC, Tracey I, Smith SM, Nichols TE.
%   Non-Parametric Combination and related permutation tests for
%   neuroimaging. Hum Brain Mapp. 2016 Apr;37(4):1486-511. (Open Access)
% ==============================================================
function Q = fasttsq(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups = 1
% - psi and res are 3D
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% Q    : Hotelling's T^2 statistic.

% Swap dimensions so that dims 1 and 2 are subjects and variables
% leaving the voxels/tests as the 3rd.
res = permute(res,[1 3 2]);
psi = permute(psi,[1 3 2]);
nT  = size(res,3);
df0 = plm.N-plm.rM{1}{m}{c}{1};
S = spr(res)/df0;
Q = zeros(1,nT);
cte2 = plm.eC{1}{m}{c}{1}'/(M'*M)*plm.eC{1}{m}{c}{1};
for t = 1:nT,
    cte1 = plm.eC{1}{m}{c}{1}'*psi(:,:,t)*plm.Dset{m}{c};
    Q(1,t) = cte1/(plm.Dset{m}{c}'*S(:,:,t)*plm.Dset{m}{c})/cte2*cte1';
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Q = fasttsq3d(M,psi,res,m,c,plm)
res = permute(res,[1 3 2]);
psi = permute(psi,[1 3 2]);
nT  = size(res,3);
df0 = plm.N-plm.rM{1}{m}{c}{1};
S = spr(res)/df0;
Q = zeros(1,nT);
for t = 1:nT,
    cte1 = plm.eC{1}{m}{c}{1}'*psi(:,:,t)*plm.Dset{m}{c};
    cte2 = plm.eC{1}{m}{c}{1}'/(M(:,:,t)'*M(:,:,t))*plm.eC{1}{m}{c}{1};
    Q(1,t) = cte1/(plm.Dset{m}{c}'*S(:,:,t)*plm.Dset{m}{c})/cte2*cte1';
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = fasttsqp(Q,df2,p)
% P-value for Hotelling's T^2
P = palm_gpval(Q*(df2-p+1)/p/df2,p,df2-p+1);

% ==============================================================
function Q = fastq(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups = 1
% - psi and res are 3D
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_core.m' and 'palm_takeargs.m'
%
% Outputs:
% Q    : Multivariate (yet scalar) statistic.

% Swap dimensions so that dims 1 and 2 are subjects and variables
% leaving the voxels/tests as the 3rd.
res = permute(res,[1 3 2]);
psi = permute(psi,[1 3 2]);
nT   = size(res,3);
cte2 = plm.eC{1}{m}{c}{1}'/(M'*M)*plm.eC{1}{m}{c}{1};
E    = spr(res);
Q    = zeros(1,nT);
for t = 1:nT,
    cte1   = plm.Dset{m}{c}'*psi(:,:,t)'*plm.eC{1}{m}{c}{1};
    H      = cte1/cte2*cte1';
    Q(1,t) = plm.qfun(plm.Dset{m}{c}'*E(:,:,t)*plm.Dset{m}{c},H);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Q = fastq3d(M,psi,res,m,c,plm)
res  = permute(res,[1 3 2]);
psi  = permute(psi,[1 3 2]);
nT   = size(res,3);
E    = spr(res);
Q    = zeros(1,nT);
for t = 1:nT,
    cte1   = plm.Dset{m}{c}'*psi(:,:,t)'*plm.eC{1}{m}{c}{1};
    cte2   = plm.eC{1}{m}{c}{1}'/(M(:,:,t)'*M(:,:,t))*plm.eC{1}{m}{c}{1};
    H      = cte1/cte2*cte1';
    Q(1,t) = plm.qfun(plm.Dset{m}{c}'*E(:,:,t)*plm.Dset{m}{c},H);
end

% ==============================================================
function Q = wilks(E,H)
% Wilks' lambda.
Q = det(E)/det(E+H);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = wilksp(Q,df1,df2,p)
r = df2-(p-df1+1)/2;
u = (p*df1-2)/4;
cden = (p^2+df1^2-5);
if cden > 0,
    t = sqrt((p^2*df1^2-4)/cden);
else
    t = 1;
end
F = (r*t-2*u)*(1-Q.^(1/t))./(Q.^(1/t)*p*df1);
Fdf1 = p*df1;
Fdf2 = r*t-2*u;
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = lawley(E,H)
% Lawley-Hotelling's trace.
Q = trace(H/E);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = lawleyp(Q,df1,df2,p)
m = (abs(p-df1)-1)/2;
n = (df2-p-1)/2;
s = min(p,df1);
if n > 0,
    b = (p+2*n)*(df1+2*n)/(2*(2*n+1)*(n-1));
    c = (2+(p*df1+2)/(b-1))/(2*n);
    Fdf1 = p*df1;
    Fdf2 = 4+(p*df1+2)/(b-1);
    F = (Q/c)*Fdf2/Fdf1;
else
    Fdf1 = s*(2*m+s+1);
    Fdf2 = 2*(s*n+1);
    F = (Q/s)*Fdf2/Fdf1;
end
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = pillai(E,H)
% Pillai's trace.
Q = trace(H/(E+H));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = pillaip(Q,df1,df2,p)
m = (abs(p-df1)-1)/2;
n = (df2-p-1)/2;
s = min(p,df1);
F = (2*n+s+1)/(2*m+s+1)*(Q./(s-Q));
Fdf1 = s*(2*m+s+1);
Fdf2 = s*(2*n+s+1);
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = roy_ii(E,H)
% Roy's (ii) largest root (analogous to F).
Q = max(eig(H/E));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = roy_iip(Q,df1,df2,p)
Fdf1 = max(p,df1);
Fdf2 = df2-Fdf1+df1;
F = Q*Fdf2/Fdf1;
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = roy_iii(E,H)
% Roy's (iii) largest root (analogous to R^2).
% No p-vals for this (not even approximate or bound).
Q = max(eig(H/(E+H)));

% ==============================================================
% function cc = cca(Y,X,k)
% % Do CCA via QR & SVD.
% % The ranks of X and Y aren't checked for speed.
% % Inputs are assumed to have been mean-centered and be free
% % of nuisance (partial CCA) via Y=Rz*Y and X=Rz*X.
% % k is the k-th CC (typically we want the 1st).
% % Based on the algorithm proposed by:
% % * Bjorck A, Golub GH. Numerical methods for
% %   computing angles between linear subspaces.
% %   Math Comput. 1973;27(123):579-579.
% [Qy,~]  = qr(Y,0);
% [Qx,~]  = qr(X,0);
% [~,D,~] = svd(Qy'*Qx,0);
% cc      = max(min(D(k,k),1),0);

% ==============================================================
function [A,B,cc] = cca(Y,X,R,S)
% Compute CCA. From cca function (line 231) of permcca.m
N = size(Y,1);
[Qy,Ry,iY] = qr(Y,0);
[Qx,Rx,iX] = qr(X,0);
K  = min(rank(Y),rank(X));
[L,D,M] = svds(Qy'*Qx,K);
cc = min(max(diag(D(:,1:K))',0),1);
A  = Ry\L(:,1:K)*sqrt(N-R);
B  = Rx\M(:,1:K)*sqrt(N-S);
A(iY,:) = A;
B(iX,:) = B;

% ==============================================================
function rpls = simpls(X,Y,k)
% Compute the correlation among the k-th pair of
% Uses the SIMPLS algorithm for partial least squares to
% compute score vectors (T and U), then provide the
% correlation between the k-th pair.
% Based on the algorithm by:
% * de Jong S. SIMPLS: An alternative approach to
%   partial least squares regression.
%   Chemom Intell Lab Syst. 1993 Mar;18(3):251-63.
[N,nCx] = size(X);
nCy = size(Y,2);
T   = zeros(N,k); U = T;
V   = zeros(nCx,k);
z   = zeros(nCy,1);
S   = X'*Y;
StS = S'*S;
N1  = N - 1;
for j = 1:k,
    StS = StS - z*z';
    [evc,~] = eig(StS);
    q = evc(:,end);
    r = S*q;
    t = X*r;
    p = X'*t;
    if N > nCx, d = r'*p; else d = t'*t; end
    d = sqrt(d/N1);
    v = p - V(:,1:max(1,j-1))*(p'*V(:,1:max(1,j-1)))';
    v = v/sqrt(v'*v);
    z = S'*v;
    S = S - v*z';
    V(:,j) = v;
    T(:,j) = t/d;
    U(:,j) = Y*q;
end
while j > 1,
    U(:,j) = U(:,j) - T(:,1:j-1)*(U(:,j)'*T(:,1:j-1)/N1)';
    j = j - 1;
end
t = T(:,k); u = U(:,k);
rpls = t'*u/sqrt((t'*t)*(u'*u));

% ==============================================================
% Below are the functions to combine statistics:
% ==============================================================
% Reference:
% * Winkler AM, Webster MA, Brooks JC, Tracey I, Smith SM, Nichols TE.
%   Non-Parametric Combination and related permutation tests for
%   neuroimaging. Hum Brain Mapp. 2016 Apr;37(4):1486-511. (Open Access)
% ==============================================================
function T = tippett(G,df1,df2)
T = min(palm_gpval(G,df1,df2),[],1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = tippettp(T,nG)
%P = T.^plm.nY;
% Note it can't be simply P = 1-(1-T)^K when implementing
% because precision is lost if the original T is smaller than eps,
% something quite common. Hence the need for the Pascal
% triangle, etc, as done below.
pw  = nG:-1:1;
cf  = pascaltri(nG);
sgn = (-1)*(-1).^pw;
P   = sgn.*cf*bsxfun(@power,T,pw');

% ==============================================================
function T = fisher(G,df1,df2)
T = -2*sum(log(palm_gpval(G,df1,df2)),1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = fisherp(T,nG)
P = palm_gpval(T,-1,2*nG);

% ==============================================================
function T = stouffer(G,df1,df2)
T = sum(palm_gtoz(G,df1,df2),1)/sqrt(size(G,1));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = stoufferp(T,~)
P = palm_gpval(T,0);

% ==============================================================
function T = wilkinson(G,df1,df2,parma)
T = sum(palm_gpval(G,df1,df2) <= parma);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = wilkinsonp(T,nG,parma)
lfac    = palm_factorial(nG);
lalpha  = log(parma);
l1alpha = log(1-parma);
P = zeros(size(T));
for k = 1:nG,
    lp1 = lfac(nG+1) - lfac(k+1) - lfac(nG-k+1);
    lp2 = k*lalpha;
    lp3 = (nG-k)*l1alpha;
    P = P + (k>=T).*exp(lp1+lp2+lp3);
end

% ==============================================================
function T = winer(G,df1,df2)
df2 = bsxfun(@times,ones(size(G)),df2);
cte = sqrt(sum(df2./(df2-2),1));
gp  = palm_gpval(G,df1,df2);
gt  = sign(gp-.5).*sqrt(df2./betainv(2*min(gp,1-gp),df2/2,.5)-df2); % =tinv(gp,df2)
T   = -sum(gt)./cte;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = winerp(T,~)
P = palm_gcdf(-T,0);

% ==============================================================
function T = edgington(G,df1,df2)
T = sum(palm_gpval(G,df1,df2),1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = edgingtonp(T,nG)
lfac = palm_factorial(nG);
fT   = floor(T);
mxfT = max(fT(:));
P = zeros(size(T));
for j = 0:mxfT,
    p1  = (-1)^j;
    lp2 = - lfac(j+1) - lfac(nG-j+1);
    lp3 = nG*log(T-j);
    P = P + (j<=fT).*p1.*exp(lp2+lp3);
end

% ==============================================================
function T = mudholkargeorge(G,df1,df2)
nG = size(G,1);
mhcte = sqrt(3*(5*nG+4)/nG/(5*nG+2))/pi;
T = mhcte*sum(log(...
    palm_gcdf(G,df1,df2)./...
    palm_gpval(G,df1,df2)),1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = mudholkargeorgep(T,nG)
P = palm_gpval(T,1,5*nG+4);

% ==============================================================
function [T,Gpval] = friston(G,df1,df2)
Gpval = palm_gpval(G,df1,df2);
T = max(Gpval,[],1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = fristonp(T,nG,parmu)
P = T.^(nG - parmu + 1);

% ==============================================================
function T = darlingtonhayes(G,df1,df2,parmr)
df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= parmr;
G       = reshape(G(idx),horzcat(parmr,size(G,2)));
df2     = reshape(df2(idx),horzcat(parmr,size(df2,2)));
Z       = palm_gtoz(G,df1,df2);
T       = mean(Z,1);

% ==============================================================
function T = zaykin(G,df1,df2,parma)
P = -log10(palm_gpval(G,df1,df2));
P(P < -log10(parma)) = 0;
T = sum(P,1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = zaykinp(T,nG,parma)
lT      = -T;
lfac    = palm_factorial(nG);
lalpha  = log10(parma);
l1alpha = log10(1-parma);
P = zeros(size(lT));
for k = 1:plm.nY,
    lp1 = lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1);
    lp2 = (plm.nY-k)*l1alpha;
    Tsmall = lT <= k*lalpha;
    Tlarge = ~ Tsmall;
    p3 = 0;
    lnum = log10(k*lalpha - lT(Tsmall));
    for j = 1:k,
        p3 = p3 + 10.^(lT(Tsmall) + (j-1).*lnum - lfac(j));
    end
    lp3small = log10(p3);
    lp3large = k*lalpha;
    P(Tsmall) = P(Tsmall) + 10.^(lp1 + lp2 + lp3small);
    P(Tlarge) = P(Tlarge) + 10.^(lp1 + lp2 + lp3large);
end

% ==============================================================
function T = dudbridgekoeleman(G,df1,df2,parmr)
df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= parmr;
G       = reshape(G(idx),horzcat(parmr,size(G,2)));
df2     = reshape(df2(idx),horzcat(parmr,size(df2,2)));
P       = -log10(palm_gpval(G,df1,df2));
T       = sum(P,1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = dudbridgekoelemanp(T,nG,parmr)
lT = -T;
lfac = palm_factorial(nG);
P    = zeros(size(lT));
lp1  = lfac(nG+1)  - ...
    lfac(parmr+2)  - ...
    lfac(nG-parmr) + ...
    log10(parmr+2);
for v = 1:numel(lT);
    P(v) = quad(@(t)dkint(t,lp1,lT(v),nG,...
        parmr,lfac(1:parmr)),eps,1);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function T = dudbridgekoeleman2(G,df1,df2,parmr,parma)
df2 = bsxfun(@times,ones(size(G)),df2);
P = -log10(palm_gpval(G,df1,df2));
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
P(tmp > parmr) = 0;
P(P < -log10(parma)) = 0;
T = sum(P,1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = dudbridgekoeleman2p(T,nG,parmr,parma)
lT = -T;
lfac = palm_factorial(nG);
P    = zeros(1,size(T,2));
for k = 1:parmr,
    kk = (nG-k)*log(1-parma);
    if isnan(kk), kk = 0; end
    p1 = exp(lfac(nG+1) - lfac(k+1) - lfac(nG-k+1) + kk);
    p2 = awtk(lT,parma,k,lfac(1:k));
    P = P + p1.*p2;
end
if k < nG,
    lp1 = lfac(nG+1)   - ...
        lfac(parmr+2)  - ...
        lfac(nG-parmr) + ...
        log(parmr+2);
    for v = 1:numel(lT);
        P(v) = P(v) + ...
            quad(@(t)dkint(t,lp1,lT(v),nG,parmr, ...
            lfac(1:parmr)),eps,parma);
    end
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function q = dkint(t,lp1,lT,K,r,lfac)
lp2 = (K-r-1).*log(1-t);
ltr = r.*log(t);
L1  = real(lp1 + lp2 + ltr);
s1  = (lT > ltr).*exp(L1);
j   = (1:r)';
lp3 = lT + (j-1)*log(r*log(t)-lT) ...
    - repmat(lfac(j),[1 numel(t)]);
L2  = real(lp1 + repmat(lp2,[r 1]) + lp3);
s2  = (lT <= ltr).*sum(exp(L2));
q   = s1 + s2;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function A = awtk(lw,t,k,lfac)
ltk = k.*log(t);
tk = real(exp(ltk));
s = (1:k)';
L = bsxfun(@plus,lw,...
    bsxfun(@minus,(s-1)*log(k*log(t)-lw),lfac(s)));
S = sum(real(exp(L)),1);
A = (lw <= ltk).*S + (lw > ltk).*tk;

% ==============================================================
function T = taylortibshirani(G,df1,df2)
nG = size(G,1);
P = palm_gpval(G,df1,df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum(1-P.*(nG+1)./prank)/nG;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = taylortibshiranip(T,nG)
P = palm_gcdf(-T./sqrt(nG),0);

% ==============================================================
function T = jiang(G,df1,df2,parma)
nG = size(G,1);
P = palm_gpval(G,df1,df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum((P<=parma).*(1-P.*(nG+1)./prank))/nG;

% ==============================================================
function [fastnpc,pparanpc,npcrev,npcrel,npcextr] = ...
    npchandles(npcmethod,concordant)
% Create the function handles for the NPC.
if nargout == 2,
    concordant = false;
end
switch lower(npcmethod),
    
    case 'tippett',
        if concordant,
            fastnpc = @(G,df1,df2)min( ...
                tippett( G,df1,df2),   ...
                tippett(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)tippett(G,df1,df2);
        end
        pparanpc    = @(T,nG)tippettp(T,nG);
        npcrev      = true;
        
    case 'fisher',
        if concordant,
            fastnpc = @(G,df1,df2)max( ...
                fisher( G,df1,df2),    ...
                fisher(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)fisher(G,df1,df2);
        end
        pparanpc    = @(T,nG)fisherp(T,nG);
        npcrev      = false;
        
    case 'stouffer',
        if concordant,
            fastnpc = @(G,df1,df2)max( ...
                stouffer( G,df1,df2),  ...
                stouffer(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)stouffer(G,df1,df2);
        end
        pparanpc    = @(T,nG)stoufferp(T,nG);
        npcrev      = false;
        
    case 'wilkinson',
        if concordant,
            fastnpc = @(G,df1,df2)max(             ...
                wilkinson( G,df1,df2,npcparm), ...
                wilkinson(-G,df1,df2,npcparm));
        else
            fastnpc = @(G,df1,df2)wilkinson(G,df1,df2,npcparm);
        end
        pparanpc    = @(T,nG)wilkinsonp(T,nG,npcparm);
        npcrev      = false;
        
    case 'winer',
        if concordant,
            fastnpc = @(G,df1,df2)max( ...
                winer( G,df1,df2),     ...
                winer(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)winer(G,df1,df2);
        end
        pparanpc    = @(T,nG)winerp(T,nG);
        npcrev      = false;
        
    case 'edgington',
        if concordant,
            fastnpc = @(G,df1,df2)min( ...
                edgington( G,df1,df2), ...
                edgington(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)edgington(G,df1,df2);
        end
        pparanpc    = @(T,nG)edgingtonp(T,nG);
        npcrev      = true;
        
    case 'mudholkar-george',
        if concordant,
            fastnpc = @(G,df1,df2)max(       ...
                mudholkargeorge( G,df1,df2), ...
                mudholkargeorge(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)mudholkargeorge(G,df1,df2);
        end
        pparanpc    = @(T,nG)mudholkargeorgep(T,nG);
        npcrev      = false;
        
    case 'friston',
        if concordant,
            fastnpc = @(G,df1,df2)min( ...
                friston( G,df1,df2),   ...
                friston(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)friston(G,df1,df2);
        end
        pparanpc    = @(T,nG)fristonp(T,nG,npcparm);
        npcrev      = true;
        
    case 'darlington-hayes',
        if concordant,
            fastnpc = @(G,df1,df2)max(                   ...
                darlingtonhayes( G,df1,df2,npcparm), ...
                darlingtonhayes(-G,df1,df2,npcparm));
        else
            fastnpc = @(G,df1,df2)darlingtonhayes(G,df1,df2,npcparm);
        end
        pparanpc    = [];
        npcrev      = false;
        
    case 'zaykin',
        if concordant,
            fastnpc = @(G,df1,df2)max(          ...
                zaykin( G,df1,df2,npcparm), ...
                zaykin(-G,df1,df2,npcparm));
        else
            fastnpc = @(G,df1,df2)zaykin(G,df1,df2,npcparm);
        end
        pparanpc    = @(T,nG)zaykinp(T,nG,npcparm);
        npcrev      = false;
        
    case 'dudbridge-koeleman',
        if concordant,
            fastnpc = @(G,df1,df2)max(                     ...
                dudbridgekoeleman( G,df1,df2,npcparm), ...
                dudbridgekoeleman(-G,df1,df2,npcparm));
        else
            fastnpc = @(G,df1,df2)dudbridgekoeleman(G,df1,df2,npcparm);
        end
        pparanpc    = @(T,nG)dudbridgekoelemanp(T,nG,npcparm);
        npcrev      = false;
        
    case 'dudbridge-koeleman2',
        if concordant,
            fastnpc = @(G,df1,df2)max(                                   ...
                dudbridgekoeleman2( G,df1,df2,npcparm,npcparm2), ...
                dudbridgekoeleman2(-G,df1,df2,npcparm,npcparm2));
        else
            fastnpc = @(G,df1,df2)dudbridgekoeleman2(G,df1,df2,npcparm,npcparm2);
        end
        pparanpc    = @(T,nG)dudbridgekoeleman2p(T,nG,npcparm,npcparm2);
        npcrev      = false;
        
    case 'taylor-tibshirani',
        if concordant,
            fastnpc = @(G,df1,df2)max(        ...
                taylortibshirani( G,df1,df2), ...
                taylortibshirani(-G,df1,df2));
        else
            fastnpc = @(G,df1,df2)taylortibshirani(G,df1,df2);
        end
        pparanpc    = @(T,nG)taylortibshiranip(T,nG);
        npcrev      = false;
        
    case 'jiang',
        if concordant,
            fastnpc = @(G,df1,df2)max(         ...
                jiang( G,df1,df2,npcparm), ...
                jiang(-G,df1,df2,npcparm));
        else
            fastnpc = @(G,df1,df2)jiang(G,df1,df2,npcparm);
        end
        npcrev      = false;
end

% For the NPC methods in which the most significant stats are the
% smallest, rather than the largest, use reverse comparisons.
if npcrev,
    npcrel  = @le;
    npcextr = @min;
else
    npcrel  = @ge;
    npcextr = @max;
end

% ==============================================================
% Other useful functions:
% ==============================================================
function savedof(df1,df2,fname)
% Save the degrees of freedom.
% This is faster than dlmwrite.
fdof = fopen(fname,'w');
fprintf(fdof,'%g\n',df1);
fprintf(fdof,'%g,',df2);
fseek(fdof,-1,'cof');
fprintf(fdof,'\n');
fclose(fdof);

% ==============================================================
function S = spr(X)
% Compute the matrix with the sum of products.
% X is a 3D array, with the resilduals of the GLM.
% - 1st dimension are the subjects
% - 2nd dimension the modalities.
% - 3rd dimension would tipically be voxels
%
% S is the sum of products that make up the covariance
% matrix:
% - 1st and 3rd dimension have the same size as the number of
%   modalities and the 2nd dimension are typically the voxels.

% To make it faster, the check should be made just once, and
% the result kept throughout runs.
persistent useway1;
if isempty(useway1),
    
    % Test both ways and compute the timings.
    tic; S1 = way1(X); w1 = toc;
    tic; S2 = way2(X); w2 = toc;
    
    % The variables sp1 and sp2 should be absolutely
    % identical but they may have sightly different numerical
    % precisions so to be consistent, choose the same that will
    % be used for all permutations later
    if w1 < w2,
        useway1 = true;
        S = S1;
    else
        useway1 = false;
        S = S2;
    end
else
    if useway1,
        S = way1(X);
    else
        S = way2(X);
    end
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function sp = way1(X)
% Way 1: this tends to be faster for Octave and if
% the number of levels in X is smaller than about 5.
[~,nY,nT] = size(X);
sp = zeros(nY,nY,nT);
for y1 = 1:nY,
    for y2 = 1:y1,
        sp(y1,y2,:) = sum(X(:,y1,:).*X(:,y2,:),1);
        if y1 ~= y2,
            sp(y2,y1,:) = sp(y1,y2,:);
        end
    end
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function sp = way2(X)
% Way 2: This tends to be faster in Matlab or if
% there are many levels in X, e.g., more than about 7.
[~,nY,nT] = size(X);
sp = zeros(nY,nY,nT);
for t = 1:nT,
    sp(:,:,t) = (X(:,:,t)'*X(:,:,t));
end

% ==============================================================
function [B,S] = lowrankfac(eC,psi,res)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups = 1
%
% Inputs:
% eC  : effective contrast
% psi : regression coefficients
% res : residuals
%
% Outputs:
% B   : p-th row of B
% S   : p-th row of S
B   = eC'*psi;
S   = sum(res.^2);

% ==============================================================
function Q = mldiv(A,B)
% This is a slower version than mldivide, which that has no
% issues with rank deficiency. Useful for the regression in the
% missing data models.
Q = pinv(A)*B;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Q = mrdiv(A,B)
% This is a slower version than mrdivide, that has no
% issues with rank deficiency. Useful for computing the statistic
% in missing data models.
Q = A*pinv(B);

% ==============================================================
function Z = yates(Y,X)
% Compute a Chi^2 test in a 2x2 contingency table, using the
% Yates correction, then convert to a z-statistic.
% Reference:
% * Yates F. Contingency tables involving small numbers and the
%   Chi^2 test. Suppl to J R Stat Soc. 1934;1(2):217-35.

% Make sure it's all binary:
Y  = Y > 0;
iX = false(size(X,1),size(X,2)+1);
for x = 1:size(X,2),
    iX(:,x) = X(:,x) > 0;
    if size(X,2) > 1 && sum(iX(:,x),1) > size(X,1)/2,
        iX(:,x) = ~ iX(:,x);
    end
end
iX(:,end) = ~ any(iX,2);

% Contingency table:
Oconf = zeros(2,size(Y,2),size(iX,2));
for x = 1:size(iX,2),
    Oconf(1,:,x) = sum(bsxfun(@and, Y, iX(:,x)),1);
    Oconf(2,:,x) = sum(bsxfun(@and,~Y, iX(:,x)),1);
end

% Margins and expected values:
margH = sum(Oconf,1);
margV = sum(Oconf,3);
Econf = bsxfun(@times,margH,margV)./size(Y,1);

% Chi^2 staistic, p-value, and z-score:
X2 = (abs(Oconf-Econf)-.5).^2./Econf;
X2 = sum(sum(X2,1),3);
P  = palm_gammainc(X2/2,size(X,2)/2,'upper'); % division of P by 2 omitted.
Z  = sqrt(2)*erfcinv(P); % multiplication of P by 2 omitted.

% ==============================================================
function C = pascaltri(K)
% Returns the coefficients for a binomial expansion of
% power K, except the last term. This is used by the Tippett
% method to avoid issues with numerical precision.
persistent Cp;
if isempty(Cp),
    K = K + 1;
    if K <= 2,
        Cp = horzcat(ones(1,K),0);
    elseif K >= 3,
        Rprev = [1 1 0];
        for r = 3:K,
            Cp = horzcat(Rprev + fliplr(Rprev),0);
            Rprev = Cp;
        end
    end
end
C = Cp(1:end-2);

% Finished! :-)
