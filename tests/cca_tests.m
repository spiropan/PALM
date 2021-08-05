addpath ../../PALM

N = 100;
nY = 16;

% create design with:
% EV1: dummy (for imaging_data1, i.e., the first -evperdat)
% EV2: dummy (for imaging_data2, i.e., the second -evperdat)
% EV3: dummy (for imaging_data3, i.e., the third -evperdat)
% EV4: meanFD
% EV5: age
% EV6: sex
% save the above as design.csv
M = [ones(N,1) ones(N,1) ones(N,1) randn(N,1) randn(N,1) rand(N,1)>.5];
%M = [randn(N,1) randn(N,1) rand(N,1)>.5];
csvwrite('design.csv',M);

% create imaging data cvs files
Y = randn(N,nY);
csvwrite('imaging_data1.csv',Y);
Y = randn(N,nY);
csvwrite('imaging_data2.csv',Y);
Y = randn(N,nY);
csvwrite('imaging_data3.csv',Y);

% create contrast files
C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];
csvwrite('tcontrasts.csv',C);

C = [1 1 1];
csvwrite('fcontrasts.csv',C)

% palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -t contrasts.csv -n 20 -corrcon -corrmod -o test_cca
% palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -n 20 -mv CCA -demean -o test_cca -nounivariate

palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv -nounivariate...
     -n 20 -t tcontrasts.csv -f fcontrasts.csv -mv CCA 2 -fonly -demean -o test_cca...
     -evperdat imaging_data1.csv 1 1 -evperdat imaging_data2.csv 2 1 -evperdat imaging_data3.csv 3 1

%delete('test_cca*')
