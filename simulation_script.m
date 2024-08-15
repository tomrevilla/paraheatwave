% Simulation script for the article:
%
% "Timing of heatwaves matters: extreme heat during parasitism disrupts
% top-down control" by Pardikes N. et al, 2024
%
% Code by Tomas A. Revilla
%
% This script loads mortality & attack data, computes parameters and simulates
% infected and non-infected host dynamics. The results are compared with data.

clear
rng(65138368);
% The line below is required by Octave, comment if using Matlab
%pkg load statistics

% Constant parameters
H0 = 50;
n_days = 8;

% Load data
[HOST,PARA,HEAT,FLYS,WASP] = textread('HWsimpleDataCleanNEW.csv', '%s %s %s %f %f', ...
                                      'delimiter', ',',  'headerlines', 1);

% Labels
lbl_host = unique(HOST); num_host = numel(lbl_host);
lbl_para = unique(PARA); num_para = numel(lbl_para);
lbl_heat = unique(HEAT); num_heat = numel(lbl_heat);

% Empirical host mortalities
mort = log(H0./max(1, FLYS))/n_days;



%%%%%% PARAMETER ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for hh = 1:num_host
  row = find(strcmp(HOST, lbl_host(hh)) & ...
             strcmp(PARA, 'NoWasp'      ) & ...
             strcmp(HEAT, 'NO'      ));
  mu0ave(hh) = mean(mort(row));
  mu0std(hh) = std(mort(row));
  ff0(hh) = mean(FLYS(row));
end

for hh = 1:num_host
  for ww = 1:num_heat
    row = find(strcmp(HOST, lbl_host(hh) ) & ...
               strcmp(PARA, 'NoWasp'       ) & ...
               strcmp(HEAT, lbl_heat(ww) ));
    mudave(hh, ww) = n_days*(mean(mort(row)) - mu0ave(hh));
    mudstd(hh, ww) = n_days*sqrt((std(mort(row)))^2 + (mu0std(hh))^2);
   end
end

for hh = 1:num_host  
  for pp = 1:num_para
    row = find(strcmp(HOST, lbl_host(hh)) & ...
               strcmp(PARA, lbl_para(pp)) & ...
               strcmp(HEAT, 'NO'      ));
    aa0ave(hh, pp) = n_days*(mean(mort(row)) - mu0ave(hh))/3;
    aa0std(hh, pp) = n_days*sqrt((std(mort(row)))^2 + (mu0std(hh))^2)/3;
    e0(hh, pp) = mean(WASP(row))/(ff0(hh) - mean(FLYS(row)));
  end
end

e0(isnan(e0)) = 0;
mudstd(:,4) = 0;
aa0std(:,4) = 0;

% Save parameters
save('parameters.mat', 'lbl_host', 'lbl_para', 'lbl_heat', ...
                       'mu0ave', 'mu0std', 'mudave', 'mudstd', ...
                       'aa0ave', 'aa0std', 'e0' )



%%%%%% SIMULATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%twave = [7, 1, 3, 0];
twave = [5, 1, 3, 0];
OUTPUT = cell(num_host, num_para, num_heat);
DATAFR = [];

for hh = 1:num_host
  for pp = 1:num_para
    for ww = 1:num_heat
      % Columns of the following matrix store the number of emerged
      % flies and wasps at the end of the experiment (n_days)
      row = find(strcmp(HOST, lbl_host(hh)) & ...
           strcmp(PARA, lbl_para(pp)) & ...
           strcmp(HEAT, lbl_heat(ww)));
      n_reps = numel(row);
      POPT = [];

      % Running the replicates
      for rr = 1:n_reps
        % Parameter sampling
        m0 = mu0ave(hh) + mu0std(hh)*randn;
        md = mudave(hh, ww) + mudstd(hh, ww)*randn;
        a0 = max(0, aa0ave(hh, pp) + aa0std(hh, pp)*randn );
        a0 = aa0ave(hh, pp) + aa0std(hh, pp)*randn;
        NVEC = [H0; 0]; % population vector of current replicate
        for tt = 1:n_days

          % Mortality rate: if day coincides with heatwave add
          % mortality differential (make result non-negative)
          mu = m0;
          if tt == twave(ww)
            mu = m0 + md;
          end
          mu = max(0, mu);
          
          % Survival test: number of surviving hosts of each type
          NVEC(1) = sum( rand(NVEC(1), 1) < exp(-mu) );
          NVEC(2) = sum( rand(NVEC(2), 1) < exp(-mu) );
        
          % Attack test: if day coincides with attack
          % copy uninfected to variable M,
          % reset uninfected,
          % add the difference to infected
          if tt == 3
            M = NVEC(1);
            NVEC(1) = sum( rand(NVEC(1), 1) < exp(-a0*3) );
            NVEC(2) = NVEC(2) + M - NVEC(1);
          end
        end
        % Concatenate the replicate column
        POPT = [POPT, [NVEC(1); poissrnd( e0(hh,pp)*NVEC(2) )] ];
        DATAFR = [DATAFR; [hh,mod(pp,4),twave(ww),POPT(1,end),POPT(2,end)] ];
      end
      % Assign replicates to appropriate experimental combination
      OUTPUT{hh, pp, ww} = POPT;
    end
  end
end
csvwrite('host_para_wave_flys_wasp.csv',DATAFR)

% Figures
deltax = 0.1;
widthx = 0.1;

f1 = figure(1);
clf
tiledlayout(3,4)
for hh = 1:num_host
  for pp = [4 1 2 3]
    X = [];
    grX =[];
    for ww = [4 2 3 1]
      Y = OUTPUT{hh, pp, ww}(1,:);
      X = [X, Y];
      grX = [grX, ww*ones(1,numel(Y))];
    end
    nexttile
    boxplot(X, grX, 'Labels',{'NoHW','HWB','HWD','HWA'}, 'Widths', widthx, 'boxstyle','filled')
    hold on
    xx = 1;
    for ww = [4 2 3 1]
      row = find(strcmp(HOST, lbl_host(hh)) & ...
                 strcmp(PARA, lbl_para(pp)) & ...
                 strcmp(HEAT, lbl_heat(ww)));
      plot(xx*ones(numel(row),1) + 2*deltax, FLYS(row), 'gx')
      xx = xx + 1;
    end
    ylabel('Number of flies')
    axis([0.5 4.5 0 50])
    title(strcat(lbl_host(hh), '-', lbl_para(pp)))
  end
end

f2 = figure(2);
clf
tiledlayout(3,3)
for hh = 1:num_host
  for pp = [1 2 3]
    X = [];
    grX =[];
    for ww = [4 2 3 1]
      Y = OUTPUT{hh, pp, ww}(1,:);
      X = [X, Y];
      grX = [grX, ww*ones(1,numel(Y))];
    end
    nexttile
    boxplot(X, grX, 'Labels',{'NoHW','HWB','HWD','HWA'}, 'Widths', widthx, 'boxstyle','filled')
    hold on
    xx = 1;
    for ww = [4 2 3 1]
      row = find(strcmp(HOST, lbl_host(hh)) & ...
                 strcmp(PARA, lbl_para(pp)) & ...
                 strcmp(HEAT, lbl_heat(ww)));
      plot(xx*ones(numel(row),1) + 2*deltax, WASP(row), 'gx')
      xx = xx + 1;
    end
    ylabel('Number of wasps')
    title(strcat(lbl_host(hh), '-', lbl_para(pp)))
    axis([0.5 4.5 0 50])
  end
end


f3 = figure(3);
clf
tiledlayout(3,4)
for hh = 1:num_host
  for pp = [4 1 2 3]
    nexttile
    hold on
    xx = 1;
    for ww = [4 2 3 1]
      X = mean(OUTPUT{hh, pp, ww}(1,:));
      Y = std(OUTPUT{hh, pp, ww}(1,:));
      plot(xx - deltax, X, 'bo')
      errorbar(xx - deltax, X, Y, 'b')
      
      row = find(strcmp(HOST, lbl_host(hh)) & ...
                 strcmp(PARA, lbl_para(pp)) & ...
                 strcmp(HEAT, lbl_heat(ww)));
      x = mean(FLYS(row));
      y = std(FLYS(row));
      plot(xx + deltax, x, 'gx')
      errorbar(xx + deltax, x, y, 'g')
      xx = xx + 1;
    end
    xticks([1 2 3 4])
    xticklabels({'NoHW','HWB','HWD','HWA'})
    ylabel('Number of flies')
    title(strcat(lbl_host(hh), '-', lbl_para(pp)))
    axis([0.5 4.5 0 50])
    box on
  end
end

f4 = figure(4);
clf
tiledlayout(3,3)
for hh = 1:num_host
  for pp = [1 2 3]
    nexttile
    hold on
    xx = 1;
    for ww = [4 2 3 1]
      X = mean(OUTPUT{hh, pp, ww}(2,:));
      Y = std(OUTPUT{hh, pp, ww}(2,:));
      plot(xx - deltax, X, 'bo')
      errorbar(xx - deltax, X, Y, 'b')
      
      row = find(strcmp(HOST, lbl_host(hh)) & ...
                 strcmp(PARA, lbl_para(pp)) & ...
                 strcmp(HEAT, lbl_heat(ww)));
      x = mean(WASP(row));
      y = std(WASP(row));
      plot(xx + deltax, x, 'gx')
      errorbar(xx + deltax, x, y, 'g')
      xx = xx + 1;
    end
    xticks([1 2 3 4])
    xticklabels({'NoHW','HWB','HWD','HWA'})
    ylabel('Number of wasps')
    title(strcat(lbl_host(hh), '-', lbl_para(pp)))
    axis([0.5 4.5 0 50])
    box on
  end
end

% Confidence interval plots for the host
% Store sample size, mean and semi-length in a file
fileIDhost = fopen('conf_int_host.csv','wt');
fprintf(fileIDhost,'%s, %s, %s, %s, %s, %s, %s\n', ...
                   'HOST', 'PARASITOID', 'HEATWAVE', ...
                   'CLASS', 'REPS', 'MEAN', 'SEMILEN');
f5 = figure(5);
clf
tiledlayout(3,4)
for hh = 1:num_host
  for pp = [4 1 2 3]
    nexttile
    hold on
    xx = 1;
    max_y = 1;
    for ww = [4 2 3 1]
      row = find(strcmp(HOST, lbl_host(hh)) & ...
           strcmp(PARA, lbl_para(pp)) & ...
           strcmp(HEAT, lbl_heat(ww)));
      n_reps = numel(row);
      X = mean(OUTPUT{hh, pp, ww}(1,:));
      Y = tinv(0.025, n_reps-1)*std(OUTPUT{hh, pp, ww}(1,:))/sqrt(n_reps);

      errorbar(xx - deltax, X, Y, 'b')
      % Write to file
      fprintf(fileIDhost, '%s, %s, %s, %s, %d, %6.4f, %6.4f\n', ...
              string(lbl_host(hh)), string(lbl_para(pp)), string(lbl_heat(ww)), ...
              'model', n_reps, X, abs(Y));

      row = find(strcmp(HOST, lbl_host(hh)) & ...
                 strcmp(PARA, lbl_para(pp)) & ...
                 strcmp(HEAT, lbl_heat(ww)));
      x = mean(FLYS(row));
      y = tinv(0.025, numel(row)-1)*std(FLYS(row))/sqrt(numel(row));
      errorbar(xx + deltax, x, y, 'g')
      % Write to file
      fprintf(fileIDhost, '%s, %s, %s, %s, %d, %6.4f, %6.4f\n', ...
              string(lbl_host(hh)), string(lbl_para(pp)), string(lbl_heat(ww)), ...
              'exper', n_reps, x, abs(y));

      text(xx + deltax, x, num2str(numel(row)))

      xx = xx + 1;
      max_y = max([max_y, abs(X)+abs(Y), abs(x)+abs(y)]);
    end
    xticks([1 2 3 4])
    xticklabels({'NoHW','HWB','HWD','HWA'})
    ylabel('Number of flies')
    plot(1, 1.05*ceil(max_y),'')
    title(strcat(lbl_host(hh), '-', lbl_para(pp)))
    axis([0.5 4.5 0 50])
    box on
  end
end

% Close the file
fclose(fileIDhost);

% Confidence interval plots for the host
% Store sample size, mean and semi-length in a file
fileIDpara = fopen('conf_int_para.csv','wt');
fprintf(fileIDpara,'%s, %s, %s, %s, %s, %s, %s\n', ...
                   'HOST', 'PARASITOID', 'HEATWAVE', ...
                   'CLASS', 'REPS', 'MEAN', 'SEMILEN');
f6 = figure(6);
clf
tiledlayout(3,3)
for hh = 1:num_host
  for pp = [1 2 3]
    nexttile
    hold on
    xx = 1;
    max_y = 1;
    for ww = [4 2 3 1]
      row = find(strcmp(HOST, lbl_host(hh)) & ...
           strcmp(PARA, lbl_para(pp)) & ...
           strcmp(HEAT, lbl_heat(ww)));
      n_reps = numel(row);
      X = mean(OUTPUT{hh, pp, ww}(2,:));
      Y = tinv(0.025, n_reps-1)*std(OUTPUT{hh, pp, ww}(2,:))/sqrt(n_reps);

      errorbar(xx - deltax, X, Y, 'b')
      % Write to file
      fprintf(fileIDpara, '%s, %s, %s, %s, %d, %6.4f, %6.4f\n', ...
              string(lbl_host(hh)), string(lbl_para(pp)), string(lbl_heat(ww)), ...
              'model', n_reps, X, abs(Y));

      row = find(strcmp(HOST, lbl_host(hh)) & ...
                 strcmp(PARA, lbl_para(pp)) & ...
                 strcmp(HEAT, lbl_heat(ww)));
      x = mean(WASP(row));
      y = tinv(0.025, numel(row)-1)*std(WASP(row))/sqrt(numel(row));

      errorbar(xx + deltax, x, y, 'g')
      % Write to file
      fprintf(fileIDpara, '%s, %s, %s, %s, %d, %6.4f, %6.4f\n', ...
              string(lbl_host(hh)), string(lbl_para(pp)), string(lbl_heat(ww)), ...
              'exper', n_reps, x, abs(y));

      text(xx + deltax, x, num2str(numel(row)))

      xx = xx + 1;
      max_y = max([max_y, abs(X)+abs(Y), abs(x)+abs(y)]);
    end
    xticks([1 2 3 4])
    xticklabels({'NoHW','HWB','HWD','HWA'})
    ylabel('Number of wasps')
    plot(1, 1.05*ceil(max_y),'')
    title(strcat(lbl_host(hh), '-', lbl_para(pp)))
    axis([0.5 4.5 0 50])
    box on
  end
end

% Close the file
fclose(fileIDpara);

%
