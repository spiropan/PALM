clear
% This file will test PALM CCA with various input data and options
% Test cases: A -   Y is imaging, X is also imaging using CSV (voxelwise)
%                   Run partial CCA with Z. Z has 3 columns (last one is zevperdat).

%             B -   Same as A except W is also defined. W has 3 columns (last
%                   one is wevperdat

%             C -   Same as A except Z is applied only to X (right side,
%                   i.e. semipartial CCA)

%             D -   Y is imaging, X is clinical variables using CSV (voxelwise)

% Add palm and permcca code to path
addpath ../../PALM
addpath ../../permcca

% Set options for testing 
new_data=0;     % Generate new test files? 
run_permcca=1;  % Run permcca to compare results to?
test='C';       % Set the test case
rng(1);         % Seed random number generator

% Options for generated CSV input files
N = 50; % Number of subjects 
nV = 8; % Number of voxels (for

% Delete old results files 
delete TEST_*.txt test_cca_dat_cca*.csv

% Generate new test data files if requested
if new_data
    % create imaging data cvs files
    Y1 = randn(N,nV); csvwrite('imaging_data1.csv',Y1);
    Y2 = randn(N,nV); csvwrite('imaging_data2.csv',Y2);
    Y3 = randn(N,nV); csvwrite('imaging_data3.csv',Y3);
    Y4 = randn(N,nV); csvwrite('imaging_data4.csv',Y4);
    Y5 = randn(N,nV); csvwrite('imaging_data5.csv',Y5);
    Y6 = randn(N,nV); csvwrite('imaging_data6.csv',Y6);
    Y7 = randn(N,nV); csvwrite('imaging_data7.csv',Y7);
    
    % create Z and W nuisance files
    Zmat = [randn(N,1) rand(N,1)>.5];
    csvwrite('Zmat.csv',Zmat);
    Z1 = randn(N,nV); csvwrite('imaging_data8.csv',Z1);
    
    Wmat = [randn(N,1) rand(N,1)>.5];
    csvwrite('Wmat.csv',Wmat);
    W1 = randn(N,nV); csvwrite('imaging_data9.csv',W1);
end

% Importdata for permcca if running permcca
if run_permcca
    Z=importdata('Zmat.csv'); W=importdata('Wmat.csv');
    Y1=importdata('imaging_data1.csv'); Y2=importdata('imaging_data2.csv');
    Y3=importdata('imaging_data3.csv'); Y4=importdata('imaging_data4.csv');
    Y5=importdata('imaging_data5.csv'); Y6=importdata('imaging_data6.csv');
    Y7=importdata('imaging_data7.csv'); Z1=importdata('imaging_data8.csv');
    W1=importdata('imaging_data9.csv');
end
               
switch test
    case 'A'     
        % Call PALM
        tic           
        palm -y imaging_data1.csv -y imaging_data2.csv -y imaging_data3.csv...
            -n 50 -o test_cca -cca 2 -z Zmat.csv -zevperdat imaging_data8.csv...
            -x imaging_data4.csv -x imaging_data5.csv -x imaging_data6.csv -x imaging_data7.csv
        palm_time = toc
        
        % Run permcca to compare to PALM if requested
        if run_permcca
            rmat=[]; wmat=[]; Amat=[]; Bmat=[]; Umat=[]; Vmat=[];
            tic
            for i=1:nV
                [punc(i,:),r{i},A{i},B{i},U{i},V{i},lW{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], [Y4(:,i) Y5(:,i) Y6(:,i) Y7(:,i)], 50, [Z Z1(:,i)]);
                rmat(i,:)=r{i}; lWmat(i,:)=lW{i}; Amat=[Amat A{i}]; Bmat=[Bmat B{i}]; Umat=[Umat U{i}]; Vmat=[Vmat V{i}];
            end
            permcca_time = toc
        end
        
        disp(['PALM took ' num2str(palm_time,'%0.2f') ' seconds'])
        if run_permcca, disp(['permcca took ' num2str(permcca_time,'%0.2f') ' seconds']); end
        
    case 'B'
        % Call PALM
        tic           
        palm -y imaging_data1.csv -y imaging_data2.csv -y imaging_data3.csv...
            -n 50 -o test_cca -cca 2 -z Zmat.csv -zevperdat imaging_data8.csv...
            -w Wmat.csv -wevperdat imaging_data9.csv...
            -x imaging_data4.csv -x imaging_data5.csv -x imaging_data6.csv -x imaging_data7.csv
        palm_time = toc
        
        % Run permcca to compare to PALM if requested
        if run_permcca
            rmat=[]; wmat=[]; Amat=[]; Bmat=[]; Umat=[]; Vmat=[]; lWmat=[];
            tic
            for i=1:nV
                %[pfwer,r(i),A,B,U,V] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)],X,20,Z,W,Sel,partial,Pset)
                [punc(i,:),r{i},A{i},B{i},U{i},V{i},lW{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], [Y4(:,i) Y5(:,i) Y6(:,i) Y7(:,i)],50,[Z Z1(:,i)],[W W1(:,i)]);
                % Store and write out outputs
                rmat(i,:)=r{i}; lWmat(i,:)=lW{i}; Amat=[Amat A{i}]; Bmat=[Bmat B{i}]; Umat=[Umat U{i}]; Vmat=[Vmat V{i}];
            end
            permcca_time = toc
        end
        
        disp(['PALM took ' num2str(palm_time,'%0.2f') ' seconds'])
        if run_permcca
            disp(['permcca took ' num2str(permcca_time,'%0.2f') ' seconds'])
        end
        
    case 'C'
        % Call PALM
        tic           
        palm -y imaging_data1.csv -y imaging_data2.csv -y imaging_data3.csv...
            -n 50 -o test_cca -cca 2 -z Zmat.csv -zevperdat imaging_data8.csv -semipartial x...
            -x imaging_data4.csv -x imaging_data5.csv -x imaging_data6.csv -x imaging_data7.csv
        palm_time = toc
        
        % Run permcca to compare to PALM if requested
        if run_permcca
            rmat=[]; wmat=[]; Amat=[]; Bmat=[]; Umat=[]; Vmat=[];
            tic
            for i=1:nV
                [punc(i,:),r{i},A{i},B{i},U{i},V{i},lW{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], [Y4(:,i) Y5(:,i) Y6(:,i) Y7(:,i)], 50, [], [Z Z1(:,i)],[],false);
                rmat(i,:)=r{i}; lWmat(i,:)=lW{i}; Amat=[Amat A{i}]; Bmat=[Bmat B{i}]; Umat=[Umat U{i}]; Vmat=[Vmat V{i}];
            end
            permcca_time = toc
        end
        
        disp(['PALM took ' num2str(palm_time,'%0.2f') ' seconds'])
        if run_permcca, disp(['permcca took ' num2str(permcca_time,'%0.2f') ' seconds']); end
        
    case 'D'
        
        % Call PALM
        % Y is imaging data, X is clinical variables, using csv files
%         palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -nounivariate...
%             -n 50 -t tcontrasts.csv -f fcontrasts.csv -mv CCA 1 -fonly -o test_cca
        
        palm -y imaging_data1.csv -x imaging_data2.csv -z Zmat.csv...
            -n 50 -o test_cca -semipartial x -cca 1 
        % palm -y fileY.ext -x fileX.ext -z fileZ1.csv -zperdat fileZ2.ext -w fileW1.csv -w fileW2.ext
        % Y_{NxP} ~ X_{NxQ}
        
        if run_permcca
            M=importdata('design.csv'); Y1=importdata('imaging_data1.csv');
            Y2=importdata('imaging_data2.csv'); Y3=importdata('imaging_data3.csv');
            
            X = M(:,[1:2]); Z = M(:,[3:4]); wmat=[]; Amat=[]; Bmat=[]; Umat=[]; Vmat=[]; punc=[]; rmat=[];
            for i=1:nV
                %[pfwer,r(i),A,B,U,V] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)],X,20,Z,W,Sel,partial,Pset)
                [punc(i,:),r{i},A{i},B{i},U{i},V{i},lW{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], X, 50, [Z Z1(:,i)]);
                
                % Store and write out outputs
                rmat(i,:)=r{i}; Wmat(i,:)=W{i}; Amat=[Amat A{i}]; Bmat=[Bmat B{i}];
                Umat=[Umat U{i}]; Vmat=[Vmat V{i}];
            end
            
        end
end

% Write out permcca files if permcca was requested
if run_permcca
    dlmwrite('TEST_permcca_punc.txt',punc,'delimiter','\t'); dlmwrite('TEST_permcca_r.txt',rmat,'delimiter','\t')
    dlmwrite('TEST_permcca_lW.txt',lWmat,'delimiter','\t'); dlmwrite('TEST_permcca_A.txt',Amat,'delimiter','\t')
    dlmwrite('TEST_permcca_B.txt',Bmat,'delimiter','\t'); dlmwrite('TEST_permcca_U.txt',Umat,'delimiter','\t')
    dlmwrite('TEST_permcca_V.txt',Vmat,'delimiter','\t')
end

%%
% Current:
% palm -> palm.m -> palm_takeargs.m  -> palm_core.m -> palm_saveall.m
% 
% Proposal:
% palm -> palm.m -> palm_takeargs.m* -> palm_prepglm.m**    -> palm_glm.m***   -> palm_saveglm.m****
%                                    -> palm_prepcca.m***** -> palm_cca.m***** -> palm_savecca.m*****
% 
% 
% *     Contains the first part of current palm_takeargs.m, and at the end would fork
% **    Contains the second part of current palm_takeargs.m (this new function can return to palm_takeargs.m)
% ***   Corresponds to current palm_core.m
% ****  Corresponds to current palm_saveall.m
% ***** Functions that need to be written (this new function can return to palm_takeargs.m)
% 
% To have the CCA functions invoked, we'd call something like this:
% palm -cca file1.ext file2.ext -transposedata -elementwise -o prefix -n 1000 -logp -eb EB.csv -cmcx -cmcp -ise -ee -nuisance nuisance1.ext nuisance2.ext  -testloadings

% For the '-testloadings' option (to be created), see 'permloads.m' on PermCCA on Github. ?

% ------ For vertexwise CCA -------------------
% palm -y fileY1.ext -y fileY2.ext -y fileY3.ext ... -x fileX1.ext -x fileX2.ext -x fileX3.ext ... -cca -z fileZ1.ext -z fileZ2.ext ... -w fileW1.ext -w fileW2.ext ... -semipartial [left,right,y,x] -theil selection.csv ... -n 10000 -o prefix -ee -ise -logp
% Each input given to -y, -x, -z, and -w, is a 4D image file, or a 2D csv file. The 4th dim (or the columns for .csv) form the variables that go in each of the internal Y, X, Z, and W.

% The -semipartial, if supplied, will use Z only with the side specified (left/y or right/x). 
% If -semipartial is not used, then Z is used for both sides (partial CCA). This -semipartial 
% option should not work if any -w is supplied (produce an error). Have synonym to -semipartial to be called -part.

% -theil is to provide a selection matrix. If -theil is used without the selection, the code should try an 
% "optimal" choice based on observations that are seldom permuted (if EB are provided) or simply make 
% a random choice (see PermCCA.m)

% ------- For spatial CCA (without transposition) --------
% palm -y fileY.ext -x fileX.ext -z fileZ1.csv -zperdat fileZ2.ext -w fileW1.csv -w fileW2.ext
% Y_{NxP} ~ X_{NxQ}

% If N are subjects and P and Q are image features (like voxels, can be images of different resolutions), 
% then spatial CCA will give mixtures of latent features that are maximally correlated, even if located in 
% non-overlapping parts of the brain.

% If -z is provided, it's treated as nuisance for all columns of Y. If -zperdat is provided, it's 
% treated as a nuisance for each column of Y (like -evperdat). Same goes for -w/-wperdat. If -zperdat 
% or -wperdat are used, we must also use -theil, not the default (Huh-Jhun). The option -semipartial is 
% as above.

% ----- For spatial CCA (with transposition) ------------
% If we transpose, then:
% Y_{PxN} ~ X_{QxN}

% P and Q must be identical. The N for each side can be different. Here CCA will give fuzzy sets of 
% subjects ("clusters") who share a similar spatial pattern across the two sides.

% Inputs similar to without transposition, but with the option -transposedata. All else the same.

%% palm_core.m pseudocode

% [opts,plm] = palm_takeargs(varargin{:}); <- current palm_takeargs lines 1-1179, 
%                                             and any others lines common to both GLM and CCA, 
%                                             including reading, organizing
%                                             masks, add case for -CCA to parse args as above
% if opts.CCA
%    [opts,plm] = palm_prepcca(opts,plm)   <- (new code)
%    [opts,plm] = palm_cca(opts,plm)       <- (new code)
%    palm_savecca(opts,plm)                <- (new code)
% else
%    [opts,plm] = palm_prepglm(opts,plm)   <- current palm_takeargs lines 1180-1553, 
%                                                     palm_core.m   lines 35-990 (remove lines about CCA and PLS)
%    [opts,plm] = palm_glm(opts,plm)       <- current palm_core.m   lines 990-2456 & 2478-end
%    palm_saveglm(opts,plm)                <- current palm_saveall.m
% end

%% June 9th notes
% Y ~ X
% nuisance variables: Z and W, respectively.
% 
% Each of these matrices can have its own number of columns.
% 
% For vertexwise CCA:
% -y fileY1.ext -y fileY3.ext -y fileY3.ext, becomes the 3 columns of Y. 
% 
% same applies to the options -x, -z, and -w. 
% 
% Each file becomes a column of the respective variable.
% 
% Since this is vertexwise, perhaps one way of doing this would be assembling internally 
% Y a 3D array of size N by nY by nV, where N is the number of subjects, nY is the number 
% of columns in Y (i.e., the number of times -y was used) and nV is the number of voxels/vertices 
% in each of the files supplied with -y. Do the same for X, Z, and W.
% 
% Then iterate over the levels (3rd dim) for each CCA, i.e., one CCA per vertex.
% 
% All files supplied with -y, -x, -z, and -w must have the same number of vertices, and must all 
% be in register (the users will take care of that during the preprocessing, we don't have to.)
% 
% Z and W can be vertexwise (say, cortical thickness), or can be constant for all vertices (say, age, sex). 
% In the latter, we can expand internally such that it becomes represented in the 3D array.

%% nov 19th items
% in palm_takeargs.m : 

% if opts.CCA
%     opts.idxout = true; % Set output of shuffree to array of permutation indeces
% ends

% If I set to false and use permutation matrices as in the Winkler et al. 2020 algorithm I don't
% get the same p-values as permCCA. Note permCCA.m uses randperm each iteration. 

% In model partition step, add the below to resize the permutation matrix/indeces 
% if opts.CCA
%     % To permute matrix with N-R rows
%     Xtmp   = Xtmp(1:size(Xtmp,1)-rank(Ztmp),:);
%     seqtmp = seqtmp(1:size(Xtmp,1),:);
% end

% in palm_core.m:
% Am only shuffling U, even through paper and permCCA shuffle both U and V
%u=plm.Qz{y}{m}{c}{1}*plm.U{m}{c}(plm.Pset(:,p),k:end,t);
%v=plm.Qw{y}{m}{c}{1}*plm.V{m}{c}(:,k:end,t);


% Mentioned including column of ones in Z
% However, this reduces rank of Z by one, since for CCA the matrix is demeaned and the
% column then becomes zero. This throws off the size calculation for R as
% written above. When including column of ones in X, I get the below
% conundrum:
% Error using palm_takeargs (line 2689)
% Contrast 1 for design 1 (and perhaps others) tests the intercept.
% This means that the options "-demean" and "-vgdemean" cannot be used.
% If "-demean" was added to calculate Pearson's "r" or the "R^2"
% note that these statistics cannot be computed for constant variables.
% Why do we want to include the colum of ones for CCA?

% note that adding -nounivariate to the function call only prevents 
% output of univariate being saved, not from being run. 
% To preclude the univariate from running I added 
% if ~ opts.CCA to line 1396 in palm_core.m
