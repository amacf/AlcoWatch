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
        function obj = AWDataSegment(time,x,y,z,gx,gy,gz,class,skipTimeCorrection)
            if nargin == 8
                skipTimeCorrection = false;
            end
            obj.time = time;
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.gx = gx;
            obj.gy = gy;
            obj.gz = gz;
            obj.class = class;
            obj = obj.removeErrors;
            if ~skipTimeCorrection
                obj = obj.correctTimeSeconds;
            end
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
        % Remove pieces of raw data that are out of order
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
        % Compute the gravityCorrectedMagnitude of the gyroscope
        % Returns a list with the magnitude of gyroscope, combined in each
        % direction
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
        function [fSegs,bSegs] = getForwardBackSegments(obj)
            movAvg = obj;
            gcmA = movAvg.gcmA;
            [maxpks,maxpksi] = findpeaks(gcmA);
            [minpks,minpksi] = findpeaks(-1*gcmA);
            minpks = -1 * minpks;
            % Maxpks, minpks gauranteed to be same length
            average = mean(gcmA(:));
            stdev = std(gcmA);
            max = average + stdev;
            min = average - stdev;
            prunedMaxi = [];
            prunedMini = [];
            for i = 1:length(maxpksi)
                if maxpks(i) > max
                    prunedMaxi = [prunedMaxi maxpksi(i)];
                end
            end
            for i = 1:length(minpksi)
                if minpks(i) < min
                    prunedMini = [prunedMini minpksi(i)];
                end
            end
            curIndex = 0;
            if prunedMaxi(1) < prunedMini(1)
                chooseMax = true;
                curIndex = prunedMaxi(1);
            else
                chooseMax = false;
                curIndex = prunedMini(1);
            end
            fSegs = [];
            bSegs = [];
            while curIndex > 0
                if chooseMax
                    nextIndex = nextSmallest(prunedMini,curIndex);
                    if nextIndex < 0
                        curIndex = nextIndex;
                        break;
                    end
                    bSegs = [bSegs AWDataSegment(obj.time(curIndex:nextIndex),obj.x(curIndex:nextIndex), obj.y(curIndex:nextIndex),obj.z(curIndex:nextIndex),obj.gx(curIndex:nextIndex),obj.gy(curIndex:nextIndex),obj.gz(curIndex:nextIndex), obj.class,true)];
                    curIndex = nextIndex;
                    chooseMax = false;
                else
                    nextIndex = nextSmallest(prunedMaxi,curIndex);
                    if nextIndex < 0
                        curIndex = nextIndex;
                        break;
                    end
                    fSegs = [fSegs AWDataSegment(obj.time(curIndex:nextIndex),obj.x(curIndex:nextIndex),obj.y(curIndex:nextIndex),obj.z(curIndex:nextIndex),obj.gx(curIndex:nextIndex),obj.gy(curIndex:nextIndex),obj.gz(curIndex:nextIndex),obj.class,true)];
                    curIndex = nextIndex;
                    chooseMax = true;
                end
            end
        end
        function [fSegs,bSegs] = getForwardBackSegmentsGyro(obj)
            movAvg = obj;
            gcmG = movAvg.gcmG;
            [maxpks,maxpksi] = findpeaks(gcmG);
            [minpks,minpksi] = findpeaks(-1*gcmG);
            minpks = -1 * minpks;
            % Maxpks, minpks gauranteed to be same length
            average = mean(gcmG(:));
            stdev = std(gcmG);
            max = average + stdev;
            min = average - stdev;
            prunedMaxi = [];
            prunedMini = [];
            for i = 1:length(maxpksi)
                if maxpks(i) > max
                    prunedMaxi = [prunedMaxi maxpksi(i)];
                end
            end
            for i = 1:length(minpksi)
                if minpks(i) < min
                    prunedMini = [prunedMini minpksi(i)];
                end
            end
            fSegs = [];
            bSegs = [];
            if isempty(prunedMini) | isempty(prunedMaxi)
                return;
            end
            curIndex = 0;
            if prunedMaxi(1) < prunedMini(1)
                chooseMax = true;
                curIndex = prunedMaxi(1);
            else
                chooseMax = false;
                curIndex = prunedMini(1);
            end
            while curIndex > 0
                if chooseMax
                    nextIndex = nextSmallest(prunedMini,curIndex);
                    if nextIndex < 0
                        curIndex = nextIndex;
                        break;
                    end
                    bSegs = [bSegs AWDataSegment(obj.time(curIndex:nextIndex),obj.x(curIndex:nextIndex), obj.y(curIndex:nextIndex),obj.z(curIndex:nextIndex),obj.gx(curIndex:nextIndex),obj.gy(curIndex:nextIndex),obj.gz(curIndex:nextIndex), obj.class,true)];
                    curIndex = nextIndex;
                    chooseMax = false;
                else
                    nextIndex = nextSmallest(prunedMaxi,curIndex);
                    if nextIndex < 0
                        curIndex = nextIndex;
                        break;
                    end
                    fSegs = [fSegs AWDataSegment(obj.time(curIndex:nextIndex),obj.x(curIndex:nextIndex),obj.y(curIndex:nextIndex),obj.z(curIndex:nextIndex),obj.gx(curIndex:nextIndex),obj.gy(curIndex:nextIndex),obj.gz(curIndex:nextIndex),obj.class,true)];
                    curIndex = nextIndex;
                    chooseMax = true;
                end
            end
        end
        function [fmed, fvar, bmed, bvar] = getRollVelocity(obj)
            % median and variance angular velocity about the angle of the 
            % arm for the forward and backward stages
            % X axis points down arm towards fingers when on left wrist
            [fsegs, bsegs] = obj.getForwardBackSegmentsGyro;
            frolls = [];
            brolls = [];
            for i = 1:length(fsegs)
                frolls = [frolls mean(fsegs(i).gx)];
            end
            for i = 1:length(bsegs)
                brolls = [brolls mean(bsegs(i).gx)];
            end
            fmed = median(frolls);
            bmed = median(brolls);
            fvar = var(frolls);
            bvar = var(brolls);
        end
        function [fmed, fvar, bmed, bvar] = getPitchVelocity(obj)
            % median and variance angular velocity about the angle of the 
            % arm relative to the body
            % Y axis point in direction of extended thumb
            % This is perhaps the most interesting because you extend
            % your arms outward (away from your body) to steady yourself
            [fsegs, bsegs] = obj.getForwardBackSegmentsGyro;
            frolls = [];
            brolls = [];
            for i = 1:length(fsegs)
                frolls = [frolls mean(fsegs(i).gy)];
            end
            for i = 1:length(bsegs)
                brolls = [brolls mean(bsegs(i).gy)];
            end
            fmed = median(frolls);
            bmed = median(brolls);
            fvar = var(frolls);
            bvar = var(brolls);
        end
        function [fmed, fvar, bmed, bvar] = getYawVelocity(obj)
            % median and variance angular velocity 
            % Z axis points out of hand towards body
            [fsegs, bsegs] = obj.getForwardBackSegmentsGyro;
            frolls = [];
            brolls = [];
            for i = 1:length(fsegs)
                frolls = [frolls mean(fsegs(i).gz)];
            end
            for i = 1:length(bsegs)
                brolls = [brolls mean(bsegs(i).gz)];
            end
            fmed = median(frolls);
            bmed = median(brolls);
            fvar = var(frolls);
            bvar = var(brolls);
        end
        function [f,b] = getRoll(obj)
            % mean angular change about the axis of the arm in forward and backward stages
            % x axis is axis parallel to the arm
            [fsegs,bsegs] = obj.getForwardBackSegmentsGyro;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = mean(fsegs(i).gx);
                roll = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + roll;
            end
            for i = 1:length(bsegs)
                vel = mean(bsegs(i).gx);
                roll = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + roll;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function [f,b] = getPitch(obj)
            % mean angular change about the y axis of the arm in forward and backward stages
            [fsegs,bsegs] = obj.getForwardBackSegmentsGyro;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = mean(fsegs(i).gy);
                roll = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + roll;
            end
            for i = 1:length(bsegs)
                vel = mean(bsegs(i).gy);
                roll = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + roll;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function [f,b] = getYaw(obj)
            % mean angular change about z axis in forward and backward stages
            [fsegs,bsegs] = obj.getForwardBackSegmentsGyro;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = mean(fsegs(i).gz);
                roll = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + roll;
            end
            for i = 1:length(bsegs)
                vel = mean(bsegs(i).gz);
                roll = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + roll;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function [fm,fv,bm,bv] = getZSpeed(obj)
            % fm = mean in speed in z direction of forward swing
            % fv = variance in speeds in z direction of forward swing
            % bm = mean in speed in z direction of backward swing
            % bv = variance in speeds in z direction of backward swing
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fvels = [];
            bvels = [];
            for i = 1:length(fsegs)
                fvels = [fvels trapz(fsegs(i).time,fsegs(i).z)];
            end
            for i = 1:length(bsegs)
                bvels = [bvels trapz(bsegs(i).time,bsegs(i).z)];
            end
            fm = mean(fvels);
            bm = mean(bvels);
            fv = var(fvels);
            bv = var(bvels);
        end
        function [fm,fv,bm,bv] = getXSpeed(obj)
            % fm = mean in speed in x direction of forward swing
            % fv = variance in speeds in x direction of forward swing
            % bm = mean in speed in x direction of backward swing
            % bv = variance in speeds in x direction of backward swing
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fvels = [];
            bvels = [];
            for i = 1:length(fsegs)
                fvels = [fvels trapz(fsegs(i).time,fsegs(i).x)];
            end
            for i = 1:length(bsegs)
                bvels = [bvels trapz(bsegs(i).time,bsegs(i).x)];
            end
            fm = mean(fvels);
            bm = mean(bvels);
            fv = var(fvels);
            bv = var(bvels);
        end
        function [fm,fv,bm,bv] = getYSpeed(obj)
            % fm = mean in speed in y direction of forward swing
            % fv = variance in speeds in y direction of forward swing
            % bm = mean in speed in y direction of backward swing
            % bv = variance in speeds in y direction of backward swing
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fvels = [];
            bvels = [];
            for i = 1:length(fsegs)
                fvels = [fvels trapz(fsegs(i).time,fsegs(i).y)];
            end
            for i = 1:length(bsegs)
                bvels = [bvels trapz(bsegs(i).time,bsegs(i).y)];
            end
            fm = mean(fvels);
            bm = mean(bvels);
            fv = var(fvels);
            bv = var(bvels);
        end
        function [f,b] = getZDist(obj)
            % Average net z displacement of up and down segments
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = trapz(fsegs(i).time,fsegs(i).z);
                dist = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + dist;
            end
            for i = 1:length(bsegs)
                vel = trapz(bsegs(i).time,bsegs(i).z);
                dist = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + dist;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function [f,b] = getXDist(obj)
            % Average net x displacement of up and down segments
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = trapz(fsegs(i).time,fsegs(i).x);
                dist = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + dist;
            end
            for i = 1:length(bsegs)
                vel = trapz(bsegs(i).time,bsegs(i).x);
                dist = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + dist;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function [f,b] = getDist(obj)
            % Average net displacement of up and down segments
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = trapz(fsegs(i).time,fsegs(i).gcmA);
                dist = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + dist;
            end
            for i = 1:length(bsegs)
                vel = trapz(bsegs(i).time,bsegs(i).gcmA);
                dist = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + dist;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function [f,b] = getYDist(obj);
            % Average net y displacement of up and down segments
            [fsegs,bsegs] = obj.getForwardBackSegments;
            fsum = 0;
            bsum = 0;
            for i = 1:length(fsegs)
                vel = trapz(fsegs(i).time,fsegs(i).y);
                dist = vel*(fsegs(i).time(length(fsegs(i).time))-fsegs(i).time(1));
                fsum = fsum + dist;
            end
            for i = 1:length(bsegs)
                vel = trapz(bsegs(i).time,bsegs(i).y);
                dist = vel*(bsegs(i).time(length(bsegs(i).time))-bsegs(i).time(1));
                bsum = bsum + dist;
            end
            f = fsum/length(fsegs);
            b = bsum/length(bsegs);
        end
        function steps = getNumSteps(obj)
            gcmA = obj.gcmA;
            pks = findpeaks(gcmA);
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
        function obj = takeMovingAverage(obj, n)
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
            fprintf(fileID, '@attribute distXForward numeric\n');
            fprintf(fileID, '@attribute distXBackward numeric\n');
            fprintf(fileID, '@attribute distYForward numeric\n');
            fprintf(fileID, '@attribute distYBackward numeric\n');
            fprintf(fileID, '@attribute distZForward numeric\n');
            fprintf(fileID, '@attribute distZBackward numeric\n');
            fprintf(fileID, '@attribute distForward numeric\n');
            fprintf(fileID, '@attribute distBackward numeric\n');
            fprintf(fileID, '@attribute rollVelMedianForward numeric\n');
            fprintf(fileID, '@attribute rollVelVarianceForward numeric\n');
            fprintf(fileID, '@attribute rollVelMedianBackward numeric\n');
            fprintf(fileID, '@attribute rollVelVarianceBackward numeric\n');
            fprintf(fileID, '@attribute pitchVelMedianForward numeric\n');
            fprintf(fileID, '@attribute pitchVelVarianceForward numeric\n');
            fprintf(fileID, '@attribute pitchVelMedianBackward numeric\n');
            fprintf(fileID, '@attribute pitchVelVarianceBackward numeric\n');
            fprintf(fileID, '@attribute yawVelMedianForward numeric\n');
            fprintf(fileID, '@attribute yawVelVarianceForward numeric\n');
            fprintf(fileID, '@attribute yawVelMedianBackward numeric\n');
            fprintf(fileID, '@attribute yawVelVarianceBackward numeric\n');
            fprintf(fileID, '@attribute rollForward numeric\n');
            fprintf(fileID, '@attribute rollBackward numeric\n');
            fprintf(fileID, '@attribute pitchForward numeric\n');
            fprintf(fileID, '@attribute pitchBackward numeric\n');
            fprintf(fileID, '@attribute yawForward numeric\n');
            fprintf(fileID, '@attribute yawBackward numeric\n');
            fprintf(fileID, '@attribute xVelMedianForward numeric\n');
            fprintf(fileID, '@attribute xVelVarianceForward numeric\n');
            fprintf(fileID, '@attribute xVelMedianBackward numeric\n');
            fprintf(fileID, '@attribute xVelVarianceBackward numeric\n');
            fprintf(fileID, '@attribute yVelMedianForward numeric\n');
            fprintf(fileID, '@attribute yVelVarianceForward numeric\n');
            fprintf(fileID, '@attribute yVelMedianBackward numeric\n');
            fprintf(fileID, '@attribute yVelVarianceBackward numeric\n');
            fprintf(fileID, '@attribute zVelMedianForward numeric\n');
            fprintf(fileID, '@attribute zVelVarianceForward numeric\n');
            fprintf(fileID, '@attribute zVelMedianBackward numeric\n');
            fprintf(fileID, '@attribute zVelVarianceBackward numeric\n');
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
            [distxf, distxb] = obj.getXDist;
            [distyf, distyb] = obj.getYDist;
            [distzf, distzb] = obj.getZDist;
            [distf,distb] = obj.getDist;
            [rvfmed, rvfvar, rvbmed, rvbvar] = obj.getRollVelocity;
            [pvfmed, pvfvar, pvbmed, pvbvar] = obj.getPitchVelocity;
            [yvfmed, yvfvar, yvbmed, yvbvar] = obj.getYawVelocity;
            [rollf,rollb] = obj.getRoll;
            [pitchf,pitchb] = obj.getPitch;
            [yawf,yawb] = obj.getYaw;
            [xspeedfm,xspeedfv,xspeedbm,xspeedbv] = obj.getXSpeed;
            [yspeedfm,yspeedfv,yspeedbm,yspeedbv] = obj.getYSpeed;
            [zspeedfm,zspeedfv,zspeedbm,zspeedbv] = obj.getZSpeed;
            fprintf(fileID,'%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',[obj.getNumSteps, obj.getCadence, skew, kurt, obj.getGaitVelocity, obj.getResidualStepLength,obj.getRatio,obj.getResidualStepTime,obj.getBandpower, obj.getSignalNoiseRatio, obj.getTotalHarmonicDistortion, obj.getXZSwayArea, obj.getXYSwayArea, obj.getYZSwayArea]);
            sVolume = obj.getSwayVolume;
            if sVolume > 0
                fprintf(fileID,',%f', sVolume);
            else 
                fprintf(fileID, ',?');
            end
            fprintf(fileID,',%f', distxf);
            fprintf(fileID,',%f', distxb);
            fprintf(fileID,',%f', distyf);
            fprintf(fileID,',%f', distyb);
            fprintf(fileID,',%f', distzf);
            fprintf(fileID,',%f', distzb);
            fprintf(fileID,',%f', distf);
            fprintf(fileID,',%f', distb);
            fprintf(fileID,',%f',rvfmed);
            fprintf(fileID,',%f',rvfvar);
            fprintf(fileID,',%f',rvbmed);
            fprintf(fileID,',%f',rvbvar);
            fprintf(fileID,',%f',pvfmed);
            fprintf(fileID,',%f',pvfvar);
            fprintf(fileID,',%f',pvbmed);
            fprintf(fileID,',%f',pvbvar);
            fprintf(fileID,',%f',yvfmed);
            fprintf(fileID,',%f',yvfvar);
            fprintf(fileID,',%f',yvbmed);
            fprintf(fileID,',%f',yvbvar);
            fprintf(fileID,',%f',rollf);
            fprintf(fileID,',%f',rollb);
            fprintf(fileID,',%f',pitchf);
            fprintf(fileID,',%f',pitchb);
            fprintf(fileID,',%f',yawf);
            fprintf(fileID,',%f',yawb);
            fprintf(fileID,',%f',xspeedfm);
            fprintf(fileID,',%f',xspeedfv);
            fprintf(fileID,',%f',xspeedbm);
            fprintf(fileID,',%f',xspeedbv);
            fprintf(fileID,',%f',yspeedfm);
            fprintf(fileID,',%f',yspeedfv);
            fprintf(fileID,',%f',yspeedbm);
            fprintf(fileID,',%f',yspeedbv);
            fprintf(fileID,',%f',zspeedfm);
            fprintf(fileID,',%f',zspeedfv);
            fprintf(fileID,',%f',zspeedbm);
            fprintf(fileID,',%f',zspeedbv);
            fprintf(fileID,',%d', height);
            fprintf(fileID,',%f', weight);
            fprintf(fileID,',%d', age);
            fprintf(fileID,',%d', gender);
            fprintf(fileID,',%d', pants);
            fprintf(fileID,',%s\n', cell2mat(obj.class));
            fclose(fileID);
        end
        % Calculate selected features and return as string
        function outString = getSelectedFeatures(obj)
            outString = '';
            outString = strcat(outString,sprintf('%d,',obj.getNumSteps));
            outString = strcat(outString,sprintf('%f,',obj.getBandpower));
            outString = strcat(outString,sprintf('%f,',obj.getTotalHarmonicDistortion));
            outString = strcat(outString,sprintf('%f,',obj.getBandpower));
            outString = strcat(outString,sprintf('%f,',obj.getXZSwayArea));
            outString = strcat(outString,sprintf('%f,',obj.getXYSwayArea));
            outString = strcat(outString,sprintf('%f,',obj.getYZSwayArea));
            [distyf, distyb] = obj.getYDist;
            outString = strcat(outString,sprintf('%f,',distyf));
            [rvfmed, rvfvar, rvbmed, rvbvar] = obj.getRollVelocity;
            outString = strcat(outString,sprintf('%f,',rvfvar));
            outString = strcat(outString,sprintf('%f,',rvbmed));
            outString = strcat(outString,sprintf('%f,',rvbvar));
            [pvfmed, pvfvar, pvbmed, pvbvar] = obj.getPitchVelocity;
            outString = strcat(outString,sprintf('%f,',pvfmed));
            outString = strcat(outString,sprintf('%f,',pvfvar));
            outString = strcat(outString,sprintf('%f,',pvbmed));
            outString = strcat(outString,sprintf('%f,',pvbvar));
            [yvfmed, yvfvar, yvbmed, yvbvar] = obj.getYawVelocity;
            outString = strcat(outString,sprintf('%f,',yvfmed));
            outString = strcat(outString,sprintf('%f,',yvfvar));
            outString = strcat(outString,sprintf('%f,',yvbmed));
            outString = strcat(outString,sprintf('%f,',yvbvar));
            [pitchf,pitchb] = obj.getPitch;
            outString = strcat(outString,sprintf('%f,',pitchf));
            outString = strcat(outString,sprintf('%f,',pitchb));
            [yawf,yawb] = obj.getYaw;
            outString = strcat(outString,sprintf('%f,',yawf));
            [xspeedfm,xspeedfv,xspeedbm,xspeedbv] = obj.getXSpeed;
            outString = strcat(outString,sprintf('%f,',xspeedbm));
            [yspeedfm,yspeedfv,yspeedbm,yspeedbv] = obj.getYSpeed;
            outString = strcat(outString,sprintf('%f,',yspeedfm));
            outString = strcat(outString,sprintf('%f,',yspeedfv));
            outString = strcat(outString,sprintf('%f',yspeedbm));
        end
    end
end
