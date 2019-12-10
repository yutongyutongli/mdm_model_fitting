function save_mat(Data, subjectNum, domain, fitbywhat, fitpar_out_path)

    if strcmp(domain, 'MON') ==1
        Datamon = Data;
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar.mat']), 'Datamon')
    elseif strcmp(domain, 'MED') ==1 && strcmp(fitbywhat, 'value') == 0
        Datamed = Data;
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar.mat']), 'Datamed')
    end