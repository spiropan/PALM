function palm_savecca(plm,opts)
% Save most of the outputs from PALM.
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

% For the MV and NPC methods in which the most significant stats are the
% smallest, rather than the largest, use reverse comparisons.
opts.saveuncorrected = true
maskinter = plm.maskinter;
% Start with the uncorrected, but don't save them yet.
% They'll be used later for the FDR.
fprintf('Computing p-values.\n');
if opts.saveuncorrected,
   if opts.cca.do
        for m = 1:1,
            for c = 1:1,
                if ~ opts.accel.noperm, % the "noperm" case is already treated
                    plm.punc{m}{c} = plm.cnt{m}{c}/plm.nP{m}(c);
                end
            end
        end
    end             
end

if opts.cca.do
    
     % canonical weights and coefficients 
    %for t=1:size(plm.A{m}{c},3)
%         for nc=1:opts.ccaorplsparm
%             ccname = [plm.Qname{m}{c} int2str(nc)];
%             palm_quicksave(shiftdim(plm.A{m}{c}(:,nc,:),2),1,opts,plm,[],[],[], ...
%                 sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_A',plm.mstr{m},plm.cstr{m}{c}));  
%         end
    %end
    % cca r p-value
    for nc=1:opts.ccaorplsparm
        nc
        if opts.ccaorplsparm > 99
            ccname = [plm.Qname{m}{c} num2str(nc,'%03d')];
        else
            ccname = [plm.Qname{m}{c} num2str(nc,'%02d')];
        end
       
        % Uncorrected p-values
        plm.maskinter = maskinter;
        if opts.saveuncorrected,
            palm_quicksave(plm.punc{m}{c}(:,nc),1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_uncp',plm.mstr{m},plm.cstr{m}{c}));
        end
        
         % Corrected p-values
        
        palm_quicksave(cummax(plm.punc{m}{c}(:,nc)),1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_fwep',plm.mstr{m},plm.cstr{m}{c}));
        
        % Save out A, B, U and V (ToDo: replace 1st line with reshape to make faster)
        tmp=[]; for t=1:size(plm.A{m}{c},3), tmp=[tmp plm.A{m}{c}(:,nc,t)]; end
        size(tmp)
        plm.maskinter = plm.maskA;
        palm_quicksave(tmp,1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_A',plm.mstr{m},plm.cstr{m}{c}));  
    
        plm.maskinter = plm.maskB;
        tmp=[]; for t=1:size(plm.B{m}{c},3), tmp=[tmp plm.B{m}{c}(:,nc,t)]; end
        palm_quicksave(tmp,1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_B',plm.mstr{m},plm.cstr{m}{c}));  
        
        plm.maskinter = plm.maskU;
        tmp=[]; for t=1:size(plm.U{m}{c},3), tmp=[tmp plm.U{m}{c}(:,nc,t)]; end
        palm_quicksave(tmp,1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_U',plm.mstr{m},plm.cstr{m}{c}));  

        plm.maskinter = plm.maskV;
        tmp=[]; for t=1:size(plm.V{m}{c},3), tmp=[tmp plm.V{m}{c}(:,nc,t)]; end
        palm_quicksave(tmp,1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_V',plm.mstr{m},plm.cstr{m}{c})); 
        

%         palm_quicksave(plm.B{m}{c}(:,nc),1,opts,plm,[],[],[], ...
%             sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_B',plm.mstr{m},plm.cstr{m}{c}));
%         
%         palm_quicksave(plm.U{m}{c}(:,nc),1,opts,plm,[],[],[], ...
%             sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_U',plm.mstr{m},plm.cstr{m}{c}));
%         
%         palm_quicksave(plm.V{m}{c}(:,nc),1,opts,plm,[],[],[], ...
%             sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,ccname,'_V',plm.mstr{m},plm.cstr{m}{c}));
    end; clear ccname nc
end   


% ==============================================================
function padj = fastfdr(pval)
% Compute FDR-adjusted p-values.
V = numel(pval);
[pval,oidx] = sort(pval);
[~,oidxR]   = sort(oidx);
padj = zeros(size(pval));
prev = 1;
for i = V:-1:1,
    padj(i) = min(prev,pval(i)*V/i);
    prev = padj(i);
end
padj = padj(oidxR);

% ==============================================================
function pvals = approxgamma(G,Gdist,rev,prepl,G1out)
% Shortcut to the calls for Gamma approximation.
if G1out,
    Gdist = Gdist(2:end,:);
end
[mu,s2,gamm1] = palm_moments(Gdist);
pvals = palm_gamma(G,mu,s2,gamm1,rev,prepl);
