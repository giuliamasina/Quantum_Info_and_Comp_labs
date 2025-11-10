% this script will calculate the results for the CHSH evaluation

    folderPath = 'C:\Users\Hp\Documents\MATLAB\Examples\R2022b\matlab\Quantum\data_Old\';
    files = dir(fullfile(folderPath, '*.txt'));
    files_names = {files.name};
    numFiles = 16;

    final_results = table('Size', [numFiles, 2], ...
                     'VariableTypes', {'string', 'double'}, ...
                     'VariableNames', {'File', 'Coincidences'});
    
    k = 1;
    for x=[0,1]
        for y=[0,1]
            for a=[0,1]
                for b=[0,1]
                    file = sprintf('x%da%dy%db%d.txt', x, a, y, b);
                    disp(file);
                    fullFilePath = fullfile(folderPath, file);
                    coincidences = process_single_file(fullFilePath);
                    disp(coincidences);
                    final_results.File(k) = file;
                    final_results.Coincidences(k) = coincidences;
                    k = k + 1;
                end
            end
        end
    end

    disp('--- Final Coincidence Count Results Table ---');
    disp(final_results);


    S_values = [];
    dS_values = [];
    coincidence_counts = final_results.Coincidences;
    num_counts = length(coincidence_counts);

    for i = 1:4:num_counts
        block_coincidences = coincidence_counts(i : i+3);
        % Map the four counts (N_pp, N_pm, N_mp, N_mm)
        Npp = block_coincidences(1); 
        Npm = block_coincidences(2); 
        Nmp = block_coincidences(3); 
        Nmm = block_coincidences(4); 
        N_total = Npp + Npm + Nmp + Nmm;

        % CHSH parameters
        E_ab   = (Npp + Nmm - Nmp - Npm) / N_total;
        E_abp  = (Npp + Nmp - Npm - Nmm) / N_total; 
        E_apb  = (Npp + Npm - Nmp - Nmm) / N_total;
        E_apbp = (Npp + Nmm - Nmp - Npm) / N_total;
        E_values = [E_ab, E_abp, E_apb, E_apbp];

        % errors for each coincidence count (poisson errors)
        delta_Npp = sqrt(Npp);
        delta_Npm = sqrt(Npm);
        delta_Nmp = sqrt(Nmp);
        delta_Nmm = sqrt(Nmm);
        delta_N = sqrt(block_coincidences);
        % CHSH error propagation
        delta_E_ab = (2 / N_total) * sqrt( (Npp*Nmm*delta_Npp^2) + (Nmp*Npm*delta_Npm^2) ) / N_total;

        delta_E_values = sqrt( (1 - E_values.^2) / N_total );
        S = E_ab + E_abp + E_apb - E_apbp;
        delta_S = sqrt(sum(delta_E_values.^2));
        S_values = [S_values; S];
        dS_values = [dS_values; delta_S];
    end

    S_final = sum(S_values);
    dS_final = sqrt(sum(dS_values.^2)); 

    fprintf('S = %.3f ± %.3f\n', S_final, dS_final);
   