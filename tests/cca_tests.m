addpath ../../PALM
addpath ../../permcca

N = 50;
nY = 8;

delete TEST_*.txt
new_data=1; run_permcca=0; % COME BACK AND RUN permcca with evperdat case

test='A'
switch test
    case 'A'
        % Create data 
        % create design with:
        % EV1: dummy (for imaging_data1, i.e., the first -evperdat)
        % EV2: dummy (for imaging_data2, i.e., the second -evperdat)
        % EV3: dummy (for imaging_data3, i.e., the third -evperdat)
        % EV4: meanFD  -  nuisance regressor
        % EV5: age     -  nuisance regressor
        % EV6: sex     -  nuisance regressor
        % save the above as design.csv
        if new_data
            M = [randn(N,1) randn(N,1) rand(N,1) randn(N,1) randn(N,1) rand(N,1) rand(N,1)>.5];
            csvwrite('design.csv',M);
            % create imaging data cvs files
            Y1 = randn(N,nY); csvwrite('imaging_data1.csv',Y1);
            Y2 = randn(N,nY); csvwrite('imaging_data2.csv',Y2);
            Y3 = randn(N,nY); csvwrite('imaging_data3.csv',Y3);
            Y4 = randn(N,nY); csvwrite('imaging_data4.csv',Y1);
            Y5 = randn(N,nY); csvwrite('imaging_data5.csv',Y2);
            Y6 = randn(N,nY); csvwrite('imaging_data6.csv',Y3);

            % create contrast files
            C=eye(6); csvwrite('tcontrasts.csv',C);
            C = [1 1 1 0 0 0]; % First 3 vars are of interest 
            csvwrite('fcontrasts.csv',C)
        end
        
        % Call PALM
        % Y is imaging data, X is also imaging data, using csv files
        palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -nounivariate...
        -n 20 -t tcontrasts.csv -f fcontrasts.csv -mv CCA 1 -fonly -demean -o test_cca... 
        -evperdat imaging_data4.csv 1 1 -evperdat imaging_data5.csv 2 1 -evperdat imaging_data6.csv 3 1
    
        if run_permcca
            M=importdata('design.csv'); Y1=importdata('imaging_data1.csv');
            Y2=importdata('imaging_data2.csv'); Y3=importdata('imaging_data3.csv');

            X = M(:,[3:4]); Z = M(:,[1:2]);
            for i=1:nY
              %[pfwer,r(i),A,B,U,V] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)],X,20,Z,W,Sel,partial,Pset)
              [punc(i,:),r{i},~,~,~,~,W{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], X, 50, Z);
              % Store first canonical correlations
              rmat(i,:)=r{i}; wmat(i,:)=W{i};
            end

            disp('Unc p-values:')
            dlmwrite('TEST_permcca_punc.txt',punc,'delimiter','\t')
            dlmwrite('TEST_permcca_r.txt',rmat,'delimiter','\t')
            dlmwrite('TEST_permcca_lW.txt',wmat,'delimiter','\t')
        end
    
    case 'B'
        if new_data
            M = [randn(N,1) randn(N,1) rand(N,1) rand(N,1)>.5];
            csvwrite('design.csv',M);

            % create imaging data cvs files
            Y1 = randn(N,nY); csvwrite('imaging_data1.csv',Y1);
            Y2 = randn(N,nY); csvwrite('imaging_data2.csv',Y2);
            Y3 = randn(N,nY); csvwrite('imaging_data3.csv',Y3);

            % create contrast files
            C=eye(4); csvwrite('tcontrasts.csv',C);
            C = [0 0 1 1];
            csvwrite('fcontrasts.csv',C)
        end
        
        % Call PALM
        % Y is imaging data, X is clinical variables, using csv files 
        palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -nounivariate...
        -n 50 -t tcontrasts.csv -f fcontrasts.csv -mv CCA 1 -fonly -o test_cca  
    
        if run_permcca
            M=importdata('design.csv'); Y1=importdata('imaging_data1.csv');
            Y2=importdata('imaging_data2.csv'); Y3=importdata('imaging_data3.csv');

            X = M(:,[3:4]); Z = M(:,[1:2]);
            for i=1:nY
              %[pfwer,r(i),A,B,U,V] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)],X,20,Z,W,Sel,partial,Pset)
              [punc(i,:),r{i},~,~,~,~,W{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], X, 50, Z);
              % Store first canonical correlations
              rmat(i,:)=r{i}; wmat(i,:)=W{i};
            end

            disp('Unc p-values:')
            dlmwrite('TEST_permcca_punc.txt',punc,'delimiter','\t')
            dlmwrite('TEST_permcca_r.txt',rmat,'delimiter','\t')
            dlmwrite('TEST_permcca_lW.txt',wmat,'delimiter','\t')
        end
end


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
