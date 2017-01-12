classdef AWDataSegment
    properties
        time
        x
        y
        z
        gx
        gy
        gz
        class
    end
    methods
        function obj = AWDataSegment(time,x,y,z,gx,gy,gz,class)
            obj.time = time;
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.gx = gx;
            obj.gy = gy;
            obj.gz = gz;
            obj.class = class;
            obj = obj.removeErrors;
            obj = obj.correctTimeSeconds;
        end
        % Remove info at given indeces
        function obj = removeI(obj, i)
            time = obj.time;
            x = obj.x;
            y = obj.y;
            z = obj.z;
            gx = obj.gx;
            gy = obj.gy;
            gz = obj.gz;
            time(i) = [];
            x(i) = [];
            y(i) = [];
            z(i) = [];
            gx(i) = [];
            gy(i) = [];
            gz(i) = [];
            obj = AWDataSegment(time,x,y,z,gx,gy,gz,obj.class);
        end
        function obj = removeErrors(obj)
            for i = 2:length(obj.time)
                if obj.time(i-1) > obj.time(i)
                    obj=obj.removeI(i-1);
                    break;
                end
            end
        end         
        function obj = removeOutliers(obj, percentRemove)
            orgNumber = length(obj.time);
            numRemove = ceil(orgNumber * (percentRemove / 100));
            % Remove in X direction of accel.
            [xsort, i] = sort(obj.x);
            obj = obj.removeI([i(1:numRemove); i(length(i)-numRemove+1:length(i))]);
            % Remove in Y direction of accel.
            [xsort, i] = sort(obj.y);
            obj = obj.removeI([i(1:numRemove); i(length(i)-numRemove+1:length(i))]);
            % Remove in Z direction of accel.
            [xsort, i] = sort(obj.z);
            obj = obj.removeI([i(1:numRemove); i(length(i)-numRemove+1:length(i))]);
            % Remove in X of gyro
            [xsort, i] = sort(obj.gx);
            obj = obj.removeI([i(1:numRemove); i(length(i)-numRemove+1:length(i))]);
            % Remove in Y of gyro
            [xsort, i] = sort(obj.gy);
            obj = obj.removeI([i(1:numRemove); i(length(i)-numRemove+1:length(i))]);
            % Remove in Z of gyro
            [xsort, i] = sort(obj.gz);
            obj = obj.removeI([i(1:numRemove); i(length(i)-numRemove+1:length(i))]);
        end
        function gcm = gravityCorrectedMagnitudeAccel(obj)
            gcm = zeros(length(obj.time), 0);
            sum = 0;
            for i = 1:length(obj.time)
                gcm(i) = sqrt((obj.x(i)).^2+(obj.y(i)).^2+(obj.z(i)).^2);
                sum = sum + gcm(i);
            end
            term2 = sum / length(obj.time);
            for i = 1:length(obj.time)
                gcm(i) = gcm(i) - term2;
            end
        end
        function gcm = gcmA(obj)
            gcm = obj.gravityCorrectedMagnitudeAccel;
        end
        function gcm = gravityCorrectedMagnitudeGyro(obj)
            gcm = zeros(length(obj.time), 0);
            sum = 0;
            for i = 1:length(obj.time)
                gcm(i) = sqrt((obj.gx(i))^2+(obj.gy(i))^2+(obj.gz(i))^2);
                sum = sum + gcm(i);
            end
            term2 = sum / length(obj.time);
            for i = 1:length(obj.time)
                gcm(i) = gcm(i) - term2;
            end
        end
        function gcm = gcmG(obj)
            gcm = obj.gravityCorrectedMagnitudeGyro;
        end
        function obj = correctTimeSeconds(obj)
            orgTime = obj.time(1);
            if orgTime == 0
                return
            end
            newTimes = obj.time;
            for i = 1:length(newTimes)
                newTimes(i) = (newTimes(i) - orgTime) / 1000;
            end
            obj = AWDataSegment(newTimes,obj.x,obj.y,obj.z,obj.gx,obj.gy,obj.gz,obj.class);
        end
        function obj = correctTimeSeg(obj)
            orgTime = obj.time(1);
            if orgTime == 0
                return
            end
            newTimes = obj.time;
            for i = 1:length(newTimes)
                newTimes(i) = (newTimes(i) - orgTime);
            end
            obj = AWDataSegment(newTimes,obj.x,obj.y,obj.z,obj.gx,obj.gy,obj.gz,obj.class);
        end
        function steps = getNumSteps(obj)
            pks = findpeaks(obj.gcmA);
            gcmA = obj.gcmA;
            average = mean(gcmA(:));
            stdev = std(gcmA);
            max = average + stdev;
            min = average - stdev;
            steps = 0;
            for i = 1:length(pks)
                if pks(i) > max | pks(i) < min
                    steps = steps + 1;
                end
            end
        end
        function loc = getStepLocations(obj)
            pks = findpeaks(obj.gcmA);
            gcmA = obj.gcmA;
            average = mean(gcmA(:));
            stdev = std(gcmA);
            max = average + stdev;
            min = average - stdev;
            loc = [];
            for i = 1:length(pks)
                if pks(i) > max | pks(i) < min
                    loc = [loc i];
                end
            end
        end
        function cadence = getCadence(obj)
            numSteps = obj.getNumSteps;
            cTime = obj.correctTimeSeg.time;
            allTime = cTime(length(cTime));
            cadence = (60 / allTime) * numSteps;
        end
        function [skew, kurt] = getSkewAndKurt(obj)
            gcmA = obj.gcmA;
            skew = skewness(gcmA);
            kurt = kurtosis(gcmA);
        end
        function area = getXZSwayArea(obj)
            ellipse = fit_ellipse(obj.gx, obj.gz);
            area = ellipse.a * ellipse.b * pi;
        end
        function area = getXYSwayArea(obj)
            ellipse = fit_ellipse(obj.gx, obj.gy);
            area = ellipse.a * ellipse.b * pi;
        end
        function area = getYZSwayArea(obj)
            ellipse = fit_ellipse(obj.gy, obj.gz);
            area = ellipse.a * ellipse.b * pi;
        end
        function volume = getSwayVolume(obj)
            volume = calculateVolume(obj.gx,obj.gy,obj.gz);
        end
        function showAccelerationData(obj)
%             seg = obj.correctTimeSeconds;
%             plot(seg.time,seg.gcmA);
            plot(obj.time,obj.gcmA);
        end
        function showGyroscopeData(obj)
            seg = obj.correctTimeSeg;
            plot(seg.time,seg.gcmG);
        end
        function sl = getStepLength(obj)
            % Calibrated to 5 seconds of walking in Olin Hall bottom floor
            % IN METERS
            segDistance = 5;
            steps = obj.getNumSteps;
            sl = segDistance / steps;
        end
        function gv = getGaitVelocity(obj)
            stepLength = obj.getStepLength;
            cadence = obj.getCadence;
            gv = cadence / stepLength;
        end
        function rsl = getResidualStepLength(obj)
            s = obj.getNumSteps;
            asl = obj.getStepLength;
            v = obj.getGaitVelocity;
            l = obj.getStepLocations;

            timeOfStep = [obj.time(l(1))];
            stepLength = v*(timeOfStep(1) - 0);
            stepLengthArray = [stepLength];

            for i=2:s
                timeOfStep = [timeOfStep obj.time(l(i))];
                stepLength = v*(timeOfStep(i) - timeOfStep(i-1));
                stepLengthArray = [stepLengthArray stepLength];

            end

            residual = asl - stepLengthArray(1);

            for j=2:s
                difference = asl - stepLengthArray(j);
                residual = [residual difference];
            end

            rsl = sum(residual);
        end
        function r = getRatio(obj)
            sp = spectralPeaks(obj.x,obj.y,obj.z);
            sp = sp(1:ceil(end/2));

            numSP = length(sp);

            totalSP = sum(sp);

            high = [];
            low = [];

            for i=1:numSP
                if(sp(i) < (mean(sp)))
                    low = [low sp(i)];
                elseif (sp(i) > (mean(sp)))
                    high = [high sp(i)];
                end
            end

            r = sum(low)./sum(high);
        end
        function obj = takeMovingAverage(obj, n);
            obj.time = movmean(obj.time,n);
            obj.x = movmean(obj.x,n);
            obj.y = movmean(obj.y,n);
            obj.z = movmean(obj.z,n);
            obj.gx = movmean(obj.gx,n);
            obj.gy = movmean(obj.gy,n);
            obj.gz = movmean(obj.gz,n);
        end
        function rst = getResidualStepTime(obj)
            s = obj.getNumSteps;
            l = obj.getStepLocations;
            timeOfStep = [obj.time(l(1))];
            stepTime = timeOfStep(1) - 0;
            stepTimeArray = [stepTime];

            for i=2:s
                timeOfStep = [timeOfStep obj.time(l(i))];
                stepTime = timeOfStep(i) - timeOfStep(i-1);
                stepTimeArray = [stepTimeArray stepTime];

            end

            totalTime = sum(stepTimeArray);

            as = s/totalTime;

            residual = as - stepTimeArray(1);

            for j=2:s
                difference = as - stepTimeArray(j);
                residual = [residual difference];
            end

            rst = sum(residual);
        end
        function bp = getBandpower(obj)
            mag = ((obj.x).^2+(obj.y).^2+(obj.z).^2);
            magNoG = mag - mean(mag);
            bp = bandpower(magNoG);
        end
        function sn = getSignalNoiseRatio(obj)
            mag = ((obj.x).^2+(obj.y).^2+(obj.z).^2);
            magNoG = mag - mean(mag);
            sn = snr(magNoG);
        end
        function td = getTotalHarmonicDistortion(obj)
            mag = ((obj.x).^2+(obj.y).^2+(obj.z).^2);
            magNoG = mag - mean(mag);
            td = thd(magNoG);
        end
        function fullFile = outputArff(obj, height, weight, age, gender, pants)
            [fileName, pathName] = uiputfile({'.arff'});
            fullFile = fullfile(pathName,fileName);
            fileID = fopen(fullFile,'w');
            fprintf(fileID, '@relation alcowatch_features\n\n');
            fprintf(fileID, '@attribute numSteps numeric\n');
            fprintf(fileID, '@attribute cadence numeric\n');
            fprintf(fileID, '@attribute skew numeric\n');
            fprintf(fileID, '@attribute kurtosis numeric\n');
            fprintf(fileID, '@attribute gaitVelocity numeric\n');
            fprintf(fileID, '@attribute residualStepLength numeric\n');
            fprintf(fileID, '@attribute ratio numeric\n');
            fprintf(fileID, '@attribute residualStepTime numeric\n');
            fprintf(fileID, '@attribute bandpower numeric\n');
            fprintf(fileID, '@attribute signalNoiseRatio numeric\n');
            fprintf(fileID, '@attribute totalHarmonicDistortion numeric\n');
            fprintf(fileID, '@attribute xzSwayArea numeric\n');
            fprintf(fileID, '@attribute xySwayArea numeric\n');
            fprintf(fileID, '@attribute yzSwayArea numeric\n');
            fprintf(fileID, '@attribute swayVolume numeric\n');
            fprintf(fileID, '@attribute height numeric\n');
            fprintf(fileID, '@attribute weight numeric\n');
            fprintf(fileID, '@attribute age numeric\n');
            fprintf(fileID, '@attribute gender {0,1,2}\n');
            fprintf(fileID, '@attribute pants {0,1,2,3}\n');
            fprintf(fileID, '@attribute class {no_goggles, green_goggles, black_goggles, red_goggles, orange_goggles}\n');
            fprintf(fileID, '\n');
            fprintf(fileID, '@data\n');
            fclose(fileID);
            obj.writeDataToArff(fullfile(pathName,fileName), height, weight, age, gender, pants);
        end
        function writeDataToArff(obj, fullFile, height, weight, age, gender, pants)
            fileID = fopen(fullFile, 'a');
            [skew, kurt] = obj.getSkewAndKurt;
            fprintf(fileID,'%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',[obj.getNumSteps, obj.getCadence, skew, kurt, obj.getGaitVelocity, obj.getResidualStepLength,obj.getRatio,obj.getResidualStepTime,obj.getBandpower, obj.getSignalNoiseRatio, obj.getTotalHarmonicDistortion, obj.getXZSwayArea, obj.getXYSwayArea, obj.getYZSwayArea]);
            sVolume = obj.getSwayVolume;
            if sVolume > 0
                fprintf(fileID,',%f', sVolume);
            else 
                fprintf(fileID, ',?');
            end
            fprintf(fileID,',%d', height);
            fprintf(fileID,',%f', weight);
            fprintf(fileID,',%d', age);
            fprintf(fileID,',%d', gender);
            fprintf(fileID,',%d', pants);
            fprintf(fileID,',%s\n', cell2mat(obj.class));
            fclose(fileID);
        end
    end
end
