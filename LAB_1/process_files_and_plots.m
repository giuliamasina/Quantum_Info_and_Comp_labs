% this script will calculate the results for the CHSH evaluation

    folderPath = 'C:\Users\Hp\Documents\MATLAB\Examples\R2022b\matlab\Quantum\data_Corrupted\';
    files = dir(fullfile(folderPath, '*.txt'));

    for i = 1:length(files)
        disp(files(i).name);
        coincidences = process_single_file(fullfile(files(i).folder, files(i).name));
        disp(coincidences);
    end