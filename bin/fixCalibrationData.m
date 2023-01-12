function T = fixCalibrationData(json_path, subject, exp_date, lab)
    % fixCalibrationData -
    % 1. reads the data from json file
    % 2. centers the data
    % 3. fixes the axes filps by inveting the corresponding scaling factor
    % 4. converts the voltage data to degree data
    % 5. Puts them in a table

    % read the json file and convert it to a struct
    json_dile = fileread(json_path);
    caldata = jsondecode(json_dile);

    % fix the order of location names
    locat = cell(5, 1);

    for i = 1:numel(fieldnames(caldata.spec.locations))
        locs = fieldnames(caldata.spec.locations);
        iloc = locs{i};
        locat(caldata.spec.locations.(iloc) + 1) = locs(i);
    end

    % centering the data
    caldata_points = caldata.points - caldata.points(caldata.spec.locations.c + 1, :);

    x_v = caldata_points(:, 1); % in voltage
    y_v = caldata_points(:, 2); % in voltage

    % ceating table with all the data
    T = table(x_v, y_v, locat);

    % check if the data is flipped
    isYfliped = T.y_v(strcmp(T.locat, 'lt')) < T.y_v(strcmp(T.locat, 'lb')); % check if top value is less than bottom value
    isXfliped = T.x_v(strcmp(T.locat, 'rb')) < T.x_v(strcmp(T.locat, 'lb')); % check if right value is less than left value

    % original scaling factor from voltage to degree
    x_scale = caldata.toVisualDegrees(1);
    y_scale = caldata.toVisualDegrees(2);

    % multiplying scaling factor by -1 if the data are flipped
    if isYfliped
        fprintf('fliping y-values \n')
        y_scale = -y_scale;
    end

    if isXfliped
        fprintf('fliping x-values \n')
        x_scale = -x_scale;
    end

    % converting voltage to degree
    x = x_v * x_scale;
    y = y_v * y_scale;

    % adding the data to the table
    T.('x_d') = x;
    T.('y_d') = y;
    T.('subject') = repmat({subject}, size(T, 1), 1);
    T.('exp_date') = repmat({exp_date}, size(T, 1), 1);
    T.('lab') = repmat(lab, size(T, 1), 1);
    T.('xscl') = repmat(x_scale, size(T, 1), 1);
    T.('yscl') = repmat(y_scale, size(T, 1), 1);
end
