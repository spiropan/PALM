addpath ../../PALM
addpath ../../permcca

N = 50;
nV = 8;

delete TEST_*.txt test_cca_dat_cca*.csv
new_data=0; run_permcca=1; % COME BACK AND RUN permcca with evperdat case

% Test cases: A - Y is imaging, X is also imaging (evperdat)
%             B - Y is imaging, X is clinical variables

test='A';
switch test
    case 'A'
        % Create data 
        % create design with:
        % EV1: dummy (for imaging_data1, i.e., the first -evperdat)
        % EV2: dummy (for imaging_data2, i.e., the second -evperdat)
        % EV3: dummy (for imaging_data3, i.e., the third -evperdat)
        % EV4: dummy (for imaging_data4, i.e., the fourth -evperdat)
        % EV5: meanFD  -  nuisance regressor
        % EV6: age     -  nuisance regressor
        % EV7: sex     -  nuisance regressor
        % save the above as design.csv
        if new_data
            M = [randn(N,1) randn(N,1) rand(N,1) randn(N,1) randn(N,1) randn(N,1) rand(N,1)>.5];
            csvwrite('design.csv',M);
            % create imaging data cvs files
            Y1 = randn(N,nV); csvwrite('imaging_data1.csv',Y1);
            Y2 = randn(N,nV); csvwrite('imaging_data2.csv',Y2);
            Y3 = randn(N,nV); csvwrite('imaging_data3.csv',Y3);
            Y4 = randn(N,nV); csvwrite('imaging_data4.csv',Y4);
            Y5 = randn(N,nV); csvwrite('imaging_data5.csv',Y5);
            Y6 = randn(N,nV); csvwrite('imaging_data6.csv',Y6);
            Y7 = randn(N,nV); csvwrite('imaging_data7.csv',Y7);

            % create contrast files
            C=eye(7); csvwrite('tcontrasts.csv',C);
            C = [1 1 1 1 0 0 0]; % First 4 vars are of interest 
            csvwrite('fcontrasts.csv',C)
        end
        
        % Call PALM
        % Y is imaging data, X is also imaging data, using csv files
        tic
        palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -nounivariate...
            -n 50 -t tcontrasts.csv -f fcontrasts.csv -mv CCA 2 -fonly -demean -o test_cca... 
            -evperdat imaging_data4.csv 1 1 -evperdat imaging_data5.csv 2 1...
            -evperdat imaging_data6.csv 3 1 -evperdat imaging_data7.csv 4 1...
       
        palm_time = toc
        
        % save out of the files for all tested canonical correlations. For
        % CCA argument the number will be the number of corrs tested. 
        % in test_cca_data_cca1 pad the number with zeros depending on size of K 01, 02, 03
    
        if run_permcca
            clear rmat wmat
            M=importdata('design.csv'); 
            Y1=importdata('imaging_data1.csv'); Y2=importdata('imaging_data2.csv'); 
            Y3=importdata('imaging_data3.csv'); Y4=importdata('imaging_data4.csv'); 
            Y5=importdata('imaging_data5.csv'); Y6=importdata('imaging_data6.csv');
            Y7=importdata('imaging_data7.csv');
            
            Z = M(:,[5:7]); wmat=[]; Amat=[]; Bmat=[]; Umat=[]; Vmat=[];
            tic
            for i=1:nV
              %[pfwer,r(i),A,B,U,V] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)],X,20,Z,W,Sel,partial,Pset)
              [punc(i,:),r{i},A{i},B{i},U{i},V{i},W{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], [Y4(:,i) Y5(:,i) Y6(:,i) Y7(:,i)], 50, Z);
              % Store and write out outputs
              rmat(i,:)=r{i}; Wmat(i,:)=W{i}; Amat=[Amat A{i}]; Bmat=[Bmat B{i}];
              Umat=[Umat U{i}]; Vmat=[Vmat V{i}];
            end
            permcca_time = toc
            
            disp('Unc p-values:')
            dlmwrite('TEST_permcca_punc.txt',punc,'delimiter','\t')
            dlmwrite('TEST_permcca_r.txt',rmat,'delimiter','\t')
            dlmwrite('TEST_permcca_lW.txt',Wmat,'delimiter','\t')
            dlmwrite('TEST_permcca_A.txt',Amat,'delimiter','\t')
            dlmwrite('TEST_permcca_B.txt',Bmat,'delimiter','\t')
            dlmwrite('TEST_permcca_U.txt',Umat,'delimiter','\t')
            dlmwrite('TEST_permcca_V.txt',Vmat,'delimiter','\t')
        end
        
        disp(['PALM took ' num2str(palm_time,'%0.2f') ' seconds'])
        if run_permcca
            disp(['permcca took ' num2str(permcca_time,'%0.2f') ' seconds'])
        end
    
    case 'B'
        if new_data
            M = [randn(N,1) randn(N,1) rand(N,1) rand(N,1)>.5];
            csvwrite('design.csv',M);

            % create imaging data cvs files
            Y1 = randn(N,nV); csvwrite('imaging_data1.csv',Y1);
            Y2 = randn(N,nV); csvwrite('imaging_data2.csv',Y2);
            Y3 = randn(N,nV); csvwrite('imaging_data3.csv',Y3);

            % create contrast files
            C=eye(4); csvwrite('tcontrasts.csv',C);
            C = [1 1 0 0]; % First two vars are of interest
            csvwrite('fcontrasts.csv',C)
        end
        
        % Call PALM
        % Y is imaging data, X is clinical variables, using csv files 
        palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -nounivariate...
        -n 50 -t tcontrasts.csv -f fcontrasts.csv -mv CCA 1 -fonly -o test_cca  
    
        if run_permcca
            M=importdata('design.csv'); Y1=importdata('imaging_data1.csv');
            Y2=importdata('imaging_data2.csv'); Y3=importdata('imaging_data3.csv');

            X = M(:,[1:2]); Z = M(:,[3:4]);
            for i=1:nV
              %[pfwer,r(i),A,B,U,V] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)],X,20,Z,W,Sel,partial,Pset)
              [punc(i,:),r{i},A{i},B{i},U{i},V{i},W{i}] = permcca([Y1(:,i) Y2(:,i) Y3(:,i)], X, 50, Z);
              % Store first canonical correlations
              rmat(i,:)=r{i}; wmat(i,:)=W{i};
            end

            disp('Unc p-values:')
            dlmwrite('TEST_permcca_punc.txt',punc,'delimiter','\t')
            dlmwrite('TEST_permcca_r.txt',rmat,'delimiter','\t')
            dlmwrite('TEST_permcca_lW.txt',wmat,'delimiter','\t')
            dlmwrite('TEST_permcca_A.txt',Amat,'delimiter','\t')
            dlmwrite('TEST_permcca_B.txt',Bmat,'delimiter','\t')
            dlmwrite('TEST_permcca_U.txt',Umat,'delimiter','\t')
            dlmwrite('TEST_permcca_V.txt',Vmat,'delimiter','\t')
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
