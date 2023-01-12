
clear all;
addpath '../toolbox/nx3k_matlab';

%%//TODO: Ask Danny about these parameters
cond = 10; % condition duration in pulses
hDeg = 20;
vDeg = 15;

% options for data cleaning 
opt.Ofs = 1000;
opt.fs = 120;
ds_rate = opt.fs/opt.Ofs;

opt.peri = round(50*ds_rate); % timepointes to cut around the blinks
opt.thr_movements = 8;
opt.thr_cluster_size = round(50*ds_rate); % 50 in ms
opt.rmSaccades_dist = round(50*ds_rate);
opt.rmSection_length = round(50*ds_rate);

% option for filtering 
opt_flt.fc = 40;
opt_flt.fs = opt.fs;
opt_flt.ord = 5;

% option for velocity detection
opt_mt.fsample = opt.fs;
opt_mt.mindur = 3 ;% round(16*ds_rate)+1;
opt_mt.velthres = 4; % lambda
opt_mt.stdGauss = 2; % std of gaussian filter operating on veloctiy data

%%
base_br_dir = dir('../data/raw/behavior/');
br_dir=base_br_dir(~ismember({base_br_dir.name},{'.','..'}));

%%%%%%%%%%%%%%%% iteration over EXPERIMENTS
for iexp = 1:length(br_dir) 
    
    expname = br_dir(iexp).name; % experiment name 
    runs_dir = dir(fullfile(br_dir(iexp).folder, expname));
    runs_dir=runs_dir(~ismember({runs_dir.name},{'.','..'}));
   
    %%%%%%%%%%%%%%%% iteration over DAYS
    for iday = 1:length(runs_dir) % iteration over days
        
        behavior_path = fullfile(runs_dir(iday).folder, runs_dir(iday).name);
        splt_name = split(runs_dir(iday).name, '_');
        namedate = splt_name{1}; % monkey name and date
        calibration_path = sprintf('../data/raw/calibration/%s', namedate);
        
        date = string(str2double(regexp(namedate,'\d*','match')'));
        monkey = {namedate(1)};
        
        fprintf('running for %s; %s \n', expname, namedate)
        
        matdata_dir = dir(fullfile(behavior_path, 'mat', '*.mat'));
        order_dir = dir(fullfile(behavior_path, 'orderFile', '*.mat'));
        eyedata_dir = dir(fullfile(behavior_path, 'fix', '*.txt'));

        % keep only the valid files
        xlsx_dir = dir(fullfile(behavior_path,  '*.xlsx'));
        acc_table = readtable(fullfile(xlsx_dir.folder,  xlsx_dir.name));
        valid_idx = acc_table.fixPerc >= 89;
        
        order_dir = order_dir(valid_idx);
        matdata_dir = matdata_dir(valid_idx);
        eyedata_dir = eyedata_dir(valid_idx);
        fixPerc = acc_table.fixPerc(valid_idx);

        % sorted array of the order files
        order_files = convertCharsToStrings(natsort({order_dir.name}));

        % extract the calibration point from json file
        calib_path = fullfile(calibration_path, 'eyecalibration.json');       
        lab = 'mri';
        calib_T = fixCalibrationData(calib_path, monkey{1}, date, lab);
        sclFactor = [calib_T.xscl(1), calib_T.yscl(1)];
              
        % name of your conditions 
        if strcmp(br_dir(iexp).name, 'Dyn_mScram')      
            conditions = {'background', ...
                'ori_faces', 'ori_bodies', 'ori_objects', ...
                'scram_faces', 'scram_bodies', 'scram_objects'};
        elseif strcmp(br_dir(iexp).name, 'Dyn_pScram')
            conditions = {'background', ...
                'ori_faces', 'ori_bodies', 'ori_objects', ...
                'pscram_faces', 'pscram_bodies', 'pscram_objects'};
        elseif strcmp(br_dir(iexp).name, 'Dyn_stat')
            conditions = {'background', ...
                'ori_faces', 'ori_bodies', 'ori_objects', ...
                'stat_faces', 'stat_bodies', 'stat_objects'};
        end
              
        data_table_all = [];
        sc_data_table_all = [];
        
        %%%%%%%%%%%%%%%% iteration over RUNS
        for iTrail = 1:length(order_dir)% iterating over the runs
            
            ord_data = load(fullfile(order_dir(iTrail).folder, order_files(iTrail)));
            all_data = load(fullfile(matdata_dir(iTrail).folder, matdata_dir(iTrail).name));
            
            % fixation data 
            fxd  = fileread(fullfile(eyedata_dir(iTrail).folder, eyedata_dir(iTrail).name));
            eval(fxd)
            
            % reformating the names
            fixdata = [];
            for i = 1:length(fixation)
                cond_name_i = fixation(i).id;
                fixdata.(cond_name_i) = fixation(i).windows(2).percent;
            end
            
            % fixing missing files in tnsData
            if isfield(all_data, 'tnsData')
                data = all_data.tnsData;
            else
                data = all_data.data;
            end
            
%              a= utils.getFixwindow(data, sclFactor)

            data_table = [];
            saccade_data_table = [];
            vel_data_table = [];
            
            %%%%%%%%%%%%%%%% iteration over CONDITIONS
            for icondt = 1:length(conditions) %iterating over the conditions
                
                fixEyeData = utils.extrAllEyeData(data, ord_data.onsets{icondt}, sclFactor);
                xdata = fixEyeData(:, 1);
                ydata = fixEyeData(:, 2);
                
                % downsample data              
                xdata_r = resample(xdata, opt.fs, opt.Ofs);
                ydata_r = resample(ydata, opt.fs, opt.Ofs);
                
                % de-meaning the data with the mean value of the background
                % condition
                if icondt == 1
                    mean_x = mean(xdata_r);
                    mean_y = mean(ydata_r);
                end 
                
                xdata = xdata_r-mean_x;
                ydata = ydata_r-mean_y;
                
                % xdata = xdata-mean(xdata); %//TODO discuss this de-meaning 
                % ydata = ydata-mean(ydata);
                
                condt = repmat(conditions(icondt), length(xdata), 1);
                data_table_i = table(xdata, ydata, condt);

                % getting the indices with clean data 
                clean_idx = utils.findCleanSections(xdata, ydata, opt);
                clean_array = utils.findCleanArray(clean_idx);
                if isempty(clean_array)
                    error('Error. \n check the basline position of eyeData \n %s %s', expname, namedate)
                end
                
                % to visualise the cleanded segements 
                plotseg = 0;
                if plotseg
                    utils.plotCleanSegments(xdata, ydata, clean_idx)
                end
                
                % cleanded eye-data in table form
                data_table_cl = data_table_i(clean_array, :);
                
                % survived signal
                survied_perc = length(clean_array)/length(xdata);  
                srv_perc = repmat(survied_perc, height(data_table_cl), 1);
                data_table_cl.('srv_perc') = srv_perc; % putting to the data table
                
                % fixation data per condition               
                if strcmp(expname, 'Dyn_stat') % dirty data names 
                    cond_name = conditions{icondt};
                else
                    cond_name = ord_data.names{icondt};               
                end
                               
                fix_val = fixdata.(cond_name);
                fix_perc = repmat(fix_val, height(data_table_cl), 1);
                data_table_cl.('fixPerc_condt') = fix_perc; % putting to the data table
                
                data_table = [data_table; data_table_cl]; % concating the data for all the conditions in a run
            end

            % de-mean here for each run (may be not required!!)
            data_table.xdata = data_table.xdata - mean(data_table.xdata);
            data_table.ydata = data_table.ydata - mean(data_table.ydata);

            % data_table_val= data_table(((abs(data_table.xdata) < thrs) & (abs(data_table.ydata)) < thrs), :);
            data_table_val = data_table;
            data_table_val.('run') = repmat(string(iTrail), height(data_table_val), 1);
            fixwindowdata = utils.getFixwindow(data, sclFactor);
            data_table_val.('xfixwin') = repmat(fixwindowdata(1),height(data_table_val), 1);
            data_table_val.('yfixwin') = repmat(fixwindowdata(2),height(data_table_val), 1);
            data_table_val.('fixPerc') = repmat(fixPerc(iTrail), height(data_table_val), 1);
                       
            saccade_data_table = utils.getSaccade(data_table_val, opt, opt_flt, opt_mt);
            sc_data_table_all = [sc_data_table_all; saccade_data_table];
      
        end
        
        sc_data_table_all.('monkey') = repmat(monkey, height(sc_data_table_all), 1);
        sc_data_table_all.('date')= repmat(date, height(sc_data_table_all), 1);
        sc_data_table_all.('exp_name') = repmat(expname, height(sc_data_table_all), 1);
      
        % save the table
        writetable(sc_data_table_all, sprintf('../results/data/data_table_%s_%s_%s.csv', expname, monkey{1}, date));     
    end
end
fprintf('finished \n')
        