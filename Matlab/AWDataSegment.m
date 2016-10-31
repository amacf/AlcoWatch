classdef AWDataSegment
    properties
        time
        x
        y
        z
        gx
        gy
        gz
        gcmA
        gcmG
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
            obj.gcmA = obj.gravityCorrectedMagnitudeAccel;
            obj.gcmG = obj.gravityCorrectedMagnitudeGyro;
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
        function steps = getNumSteps(obj)
            pks = findpeaks(obj.gcmA);
            average = mean(obj.gcmA(:));
            stdev = std(obj.gcmA);
            max = average + stdev;
            min = average - stdev;
            steps = 0;
            for i = 1:length(pks)
                if pks(i) > max | pks(i) < min
                    steps = steps + 1;
                end
            end
        end
        function cadence = getCadence(obj)
            numSteps = obj.getNumSteps;
            cTime = obj.correctTimeSeconds.time;
            allTime = cTime(length(cTime));
            cadence = (60 / allTime) * numSteps;
        end
        function [skew, kurt] = getSkewAndKurt(obj)
            skew = skewness(obj.gcmA);
            kurt = kurtosis(obj.gcmA);
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
            [ellipsoidCenter, ellipsoidRadii] = ellipsoid_fit([obj.gx obj.gy obj.gz]);
            volume = abs((4/3) * pi * ellipsoidRadii(1) * ellipsoidRadii(2) * ellipsoidRadii(3));
        end
        function showAccelerationData(obj)
%             seg = obj.correctTimeSeconds;
%             plot(seg.time,seg.gcmA);
            data = zeros(length(obj.time), 0);
            for i = 1:length(obj.time)
                data(i) = sqrt((obj.x(i))^2+(obj.y(i))^2+(obj.z(i))^2);
            end
            plot(obj.time,obj.gcmA);
        end
        function showGyroscopeData(obj)
            seg = obj.correctTimeSeconds;
            plot(seg.time,seg.gcmG);
        end
        function fullFile = outputArff(obj)
            [fileName, pathName] = uiputfile({'.arff'});
            fullFile = fullfile(pathName,fileName);
            fileID = fopen(fullFile,'w');
            fprintf(fileID, '@relation alcowatch_features\n\n');
            fprintf(fileID, '@attribute numSteps numeric\n');
            fprintf(fileID, '@attribute cadence numeric\n');
            fprintf(fileID, '@attribute skew numeric\n');
            fprintf(fileID, '@attribute kurtosis numeric\n');
            fprintf(fileID, '@attribute xzSwayArea numeric\n');
            fprintf(fileID, '@attribute xySwayArea numeric\n');
            fprintf(fileID, '@attribute yzSwayArea numeric\n');
            fprintf(fileID, '@attribute swayVolume numeric\n');
            fprintf(fileID, '@attribute class {no_goggles, green_goggles, black_goggles, red_goggles, orange_goggles}\n');
            fprintf(fileID, '\n');
            fprintf(fileID, '@data\n');
            fclose(fileID);
            obj.writeDataToArff(fullfile(pathName,fileName));
        end
        function writeDataToArff(obj, fullFile)
            fileID = fopen(fullFile, 'a');
            [skew, kurt] = obj.getSkewAndKurt
            fprintf(fileID,'%d,%f,%f,%f,%f,%f,%f,%f',[obj.getNumSteps, obj.getCadence, skew, kurt, obj.getXZSwayArea, obj.getXYSwayArea, obj.getYZSwayArea, obj.getSwayVolume]);
            fprintf(fileID,',%s\n', cell2mat(obj.class));
            fclose(fileID);
        end
    end
end
