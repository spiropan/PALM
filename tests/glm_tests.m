addpath ../../PALM

N = 50;
nV = 8;

delete test_glm_dat*.csv
new_data=0;    

% Test cases: A - evperdat
%             B - not evperdat

test='B';
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
        tic
        palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv...
            -n 50 -t tcontrasts.csv -f fcontrasts.csv -fonly -demean -o test_glm... 
            -evperdat imaging_data4.csv 1 1 -evperdat imaging_data5.csv 2 1...
            -evperdat imaging_data6.csv 3 1 -evperdat imaging_data7.csv 4 1...
       
        palm_time = toc
        disp(['PALM took ' num2str(palm_time,'%0.2f') ' seconds'])
    
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
        palm -i imaging_data1.csv -i imaging_data2.csv -i imaging_data3.csv -d design.csv...
        -n 50 -t tcontrasts.csv -f fcontrasts.csv -fonly -o test_glm  
    
end


