function [json_path, lab_found] = locateCalib(subject, date)
    % locates the path of calibration json file in the L-Drive backup folder

    labs = [4058, 4054, 4066, 4050];
    backup_path = '/media/u0139960/L-drive/GBW-0055_RVogels_Backup/Backup/Data/';

    % first locate the calibration directory
    calib_path = [];

    for lab = labs
        dpath = fullfile(backup_path, sprintf('setup-efysio-10-%d/calibvis/%s/%s*', lab, date, upcase(subject)));
        calib_path = dir(dpath);

        if ~isempty(calib_path)
            lab_found = lab;
            fprintf('Found calibration file in the lab %d \n', lab_found)
            break;
        end

    end

    % check 4054 lab data (has differenct directory sturcture)
    if isempty(calib_path)
        dpath = fullfile(backup_path, sprintf('setup-efysio-10-4054/calibvis/data/%s/%s*', date, upcase(subject)));
        calib_path = dir(dpath);

        if ~isempty(calib_path)
            lab_found = 4054;
            fprintf('Found calibration file in the lab %d \n', lab_found)
        end

    end

    % throw an error if the directory not found in the L-Drive
    if isempty(calib_path)
        error('No calibration file available in the L-Drive \n')
    end

    % check if the directory has multiple calibration instances
    if numel(calib_path) > 1
        warning('found %d calibration files, selecting the last one', numel(calib_path))
        %         calib_path = calib_path(end);
    end

    % look for the json file in the last calibration instance
    json_path = fullfile(calib_path(end).folder, calib_path(end).name, 'eyecalibration.json');

    % if unavailabe look for json file in earlier instances
    if ~isfile(json_path)
        json_path = fullfile(calib_path(end - 1).folder, calib_path(end - 1).name, 'eyecalibration.json');
    end

    % throw an error if json file not found in the L-Drive
    if ~isfile(json_path)
        error('No json file available \n');
    end

end
