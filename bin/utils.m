classdef utils
    % static methods for the analysis
    methods(Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% extract data for one onset %%%%%%%%%%%%%
        function XY_inDeg = extrEyeData(data, onset, sclFactor)
            cond_init = 10;
            if onset == 145
                cond = cond_init-5; %//TODO:check this with Danny
            else
                cond = cond_init;
            end

            X1 = tns.extractanalog('Eye X', data.stimuliPerNmrPulse(onset).pulseTime,...
             data.stimuliPerNmrPulse(onset+cond).pulseTime, data);

            Y1 = tns.extractanalog('Eye Y', data.stimuliPerNmrPulse(onset).pulseTime,...
             data.stimuliPerNmrPulse(onset+cond).pulseTime, data);

            XY = [X1 Y1];
            XY_inDeg = XY./sclFactor;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% extract data for all onset (conditions) %%%%%%%%%%%%%
        function XY_inDeg = extrAllEyeData(data, onsets, sclFactor)
            XY_inDeg = [];
            for i = 1:length(onsets)
                iXY_inDeg= utils.extrEyeData(data, onsets(i), sclFactor);
                XY_inDeg = [XY_inDeg; iXY_inDeg];
                % sprintf('Extracted eye data for condition %d', i)
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% remove saccades %%%%%%%%%%%%%
        function sections = findCleanSections(xdata, ydata,opt)
            % outputes the range of indices without the noise clusters.
            % also removes the clean segments shorter than the cluster threshold
            % length.

            % peri = time points to remove before and after of a cluster of noises.
            % thr_val = threshold to determinve the noise in eye movements
            % thr_cluster_size = minmum cluster size to keep

            thr_val = opt.thr_movements;
            peri = opt.peri;
            thr_cluster = opt.thr_cluster_size;

            x_indices = find(abs(xdata)>thr_val);
            y_indices = find(abs(ydata)>thr_val);

            c_all = unique(sort([x_indices; y_indices]));
            if isempty(c_all)
                sections.start = 1;
                sections.end = length(xdata);
            else
                diff_c_all = diff(c_all);
                new_idx = diff_c_all < thr_cluster;
                new_new_idx = find(~new_idx);
                st_idx = c_all(new_new_idx)+peri;
                end_idx = c_all(new_new_idx+1)-peri;
                strt_idx = [1; st_idx; c_all(end)+peri];
                ends_idx = [c_all(1)-peri;end_idx; length(xdata)];

                % find and remove the clean segments less than the cluster length 
                rm_idx = ends_idx -  strt_idx < opt.rmSection_length;
                strt_idx(rm_idx) = [];
                ends_idx(rm_idx) = [];

                sections.start= strt_idx;
                sections.end = ends_idx;
            end
        end  



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% indices of the cleaned signal %%%%%%%%%%%%%
        function cln_array = findCleanArray(clean_idx)
            cln_array = [];
            for i = 1:length(clean_idx.start)
                cln_array = [cln_array, clean_idx.start(i):clean_idx.end(i)];
            end
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% plot the cleaned signal %%%%%%%%%%%%%
       function plotCleanSegments(xdata, ydata, clean_idx)
            figure('position', [50, 90, 1500, 300]);
            plot(xdata)
            hold on ;
            plot(ydata);
            hold on;
            
            h = arrayfun(@(a)xline(a),clean_idx.start);
            h = arrayfun(@(a)xline(a),clean_idx.end);
%             
%             hold on;
%             clean_array = utils.findCleanArray(clean_idx);
%             plot(xdata(clean_array));
%             plot(ydata(clean_array));
%             hold off;        
       end



       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Butterworth filtering the data %%%%%%%%%%%%%
        function flt_dat = BWfilter(dat_cln, opt)
            
            cln_xdata = dat_cln(:, 1);
            cln_ydata = dat_cln(:, 2);
            
            % removing mean
            cln_xdata = cln_xdata - mean(cln_xdata);
            cln_ydata = cln_ydata - mean(cln_ydata);
   
            [b, a] = butter(opt.ord,opt.fc/(opt.fs/2));
            xdata_flt = filter(b, a, cln_xdata);
            ydata_flt = filter(b, a, cln_ydata);
            
            flt_dat = [xdata_flt, ydata_flt];
        end


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% get the fixation window form data %%%%%%%%%%%%% 
        function ab = getFixwindow(data, sclFactor)
            ab= [];
            data_length = length(data.EyeWindows);
            for i = 1:data_length
                a = data.EyeWindows(i).windows;
                x = a(3)-a(1);
                y = a(4)-a(2);
                ab(i, :) = [x/sclFactor(1), y/sclFactor(2)];
%                 ab(i, :)= [x, y];
%                 mean_fixwindow = median(ab);
            end
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% get the veloctiy from Victor's method %%%%%%%%%%%%%
        function vel = getVel(dat)
            % function to calculate velocity from position data
            % based on victor's code
            % input: dat = [x; y] position data
            % output: vel = [vx; vy] velocity data
    
            xdata = dat(1, :);
            ydata = dat(2, :);

            ydata_new = [0 ydata 0 0];
            xdata_new = [0 xdata 0 0];

            for k = 2:length(xdata_new)-2
                vxdata(:, k-1) = (sum(xdata_new([k+2, k+1])) - sum(xdata_new([k, k-1])))*(120/4);
                vydata(:, k-1) = (sum(ydata_new([k+2, k+1])) - sum(ydata_new([k, k-1])))*(120/4);
            end

            vel(1, :) = vxdata;
            vel(2,:) = vydata;
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% get the movement data %%%%%%%%%%%%%
        function [movement, vlt] = getMovements(filt_dat, opt)
            % adopted from the function of the same name in the FieldTrip toolbox
            % depndencies: getVel.m
            % INPUTS:

            dat = filt_dat';   
            method = 'victor'; % or 
        %     method = 'Engbert';
            
            if method == 'victor'
                vel = utils.getVel(dat);
            elseif method == 'engbert'    
                time = (1:length(dat))/opt.fsample;
                kernel = [1 1 0 -1 -1].*(opt.fsample/6); % original kernel (Engbert et al., 2003)
                % kernel = [1 0 -1].*(opt.fsample/3); % adopted for 3 samples
                vel =  convn(dat,   kernel,   'same');
            end
            
            %% microsaccade detection
            % compute velocity thresholds as in Engbert et al (2003) Vis Res, eqn. (2)
            medianstd = sqrt( median(vel.^2,2) - (median(vel,2)).^2 );

            % Engbert et al (2003) Vis Res, eqn. (3)
            radius =opt.velthres*medianstd;
            ndatsample = size(dat,2);

            % compute test criterion: ellipse equation
            test = sum((vel./radius(:,ones(1,ndatsample))).^2,1);
            sacsmp = find(test>1);% microsaccade's indexing

            %% determine microsaccades per trial
            % first find eye movements of n-consecutive time points
            j = find(diff(sacsmp)==1);
            j1 = [j; j+1];
            com = intersect(j,j+1);
            cut = ~ismember(j1,com);
            sacidx = reshape(j1(cut),2,[]);

            movement = [];
            peak_velComp = [];
            peak_vel = [];

            for k=1:size(sacidx,2)
            duration = sacidx(1,k):sacidx(2,k);
                if size(duration,2) >= opt.mindur 
                % finding peak velocity by Pitagoras
                begtrl = sacsmp(duration(1,1));
                endtrl = sacsmp(duration(1,end));

                [peakvel, smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
                veltrl = sacsmp(duration(1,smptrl));% peak velocity microsaccade sample -> important for spike conversion

                trlsmp = 1:ndatsample;
                begsample = trlsmp(1, begtrl); % begining microsaccade sample
                endsample = trlsmp(1, endtrl); % end microsaccade sample
                velsample = trlsmp(1, veltrl); % velocity peak microsaccade sample
                movement(end+1,:) = [begsample endsample velsample];
                end
            end   
            vlt = vel'; 
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% wrapper function to get the saccade data %%%%%%%%%%%%%
        function data_table_cl =   getSaccade(data_table_cl, opt, opt_flt, opt_mt)
            % funtction to get saccade data from the position data
            % and to store them in a data table 

            % clean data after removeing the blinks
            dat_cln = [data_table_cl.xdata, data_table_cl.ydata];

            % filtering position data
            flt_data= utils.BWfilter(dat_cln, opt_flt);

            % getting the movement information 
            [mv, vl] = utils.getMovements(flt_data, opt_mt);

            % saccade data
            Xvel = vl(:, 1);
            Yvel = vl(:, 2);
            sacc = zeros(length(Xvel), 1);
            dist = zeros(length(Xvel), 1);
            if ~isempty(mv)                 
                % remove saccades that are close to each other
                % follows the victor's paper
                rm_idx = find(diff(mv(:, 3))<opt.rmSaccades_dist);
                mv(rm_idx+1, :) = [];
                rm_idx2 = find(diff(mv(:, 3))<opt.rmSaccades_dist);
                if ~isempty(rm_idx2)
                    error('found saccades in the proximity of %d ms \n', opt.rmSaccades_dist)
                end

                % tagging saccades and distances 
                vel_idx = mv(:, 3);                   
                sacc(vel_idx) = 1;       
                distance_ini =  sqrt(sum(flt_data.^2, 2));
                dist(vel_idx) = abs(distance_ini(mv(:, 1)) - distance_ini(mv(:, 2))); % distance in origingal data for saccades;
            end

            data_table_cl.('Xvel') = Xvel;
            data_table_cl.('Yvel') = Yvel;
            data_table_cl.('sacc') = sacc;
            data_table_cl.('dist') = dist;
            data_table_cl.('xvel_flt') = utils.gaussFilt(Xvel, opt_mt);
            data_table_cl.('yvel_flt') = utils.gaussFilt(Yvel, opt_mt);
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Gaussina filter kernal %%%%%%%%%%%%%
        function fval = gaussFilt(val, opt_mt)
            FilterSize = 13;
            std = opt_mt.stdGauss;
            gfilter_1d=fspecial('gaussian',[1 FilterSize], std);
        %     figure; plot(gfilter_1d)
            fval=imfilter(val,gfilter_1d);
        end

        
            
        %%% static clss ends here %%%          
    end
end

