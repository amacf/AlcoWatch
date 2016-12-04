classdef AWDataFile
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
        function obj = AWDataFile(time,x,y,z,gx,gy,gz,class)
            obj.time = time;
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.gx = gx;
            obj.gy = gy;
            obj.gz = gz;
            obj.class = class;
        end 
        function segs = segmentsWithSize(obj, segSize)
            numSegs = floor(length(obj.time)/segSize)
            segs = AWDataSegment.empty(0, numSegs)
            for i = 1:numSegs
                segStart = (i-1)*segSize + 1
                segEnd = i * segSize
                segs(i) = AWDataSegment(obj.time(segStart:segEnd), obj.x(segStart:segEnd), obj.y(segStart:segEnd), obj.z(segStart:segEnd), obj.gx(segStart:segEnd), obj.gy(segStart:segEnd), obj.gz(segStart:segEnd), obj.class(1));
            end
        end
        function showAccelerationData(obj)
            l = length(obj.time);
            seg = obj.segmentsWithSize(l);
            seg = seg.correctTimeSeconds;
            plot(seg.time,seg.gcmA);
        end
        function showGyroscopeData(obj)
            l = length(obj.time);
            seg = obj.segmentsWithSize(l);
            seg = seg.correctTimeSeconds;
            plot(seg.time,seg.gcmG);
        end
    end
    methods (Static)
        function obj = AWDataFileFromFile(filename)
            if(nargin==0)
                [name, path] = uigetfile({'*.arff'});
                filename = fullfile(path,name);
            end
            delimiter = ',';
            startRow = 10;
            timeformatSpec = '%s%[^\n\r]';
            otherformatSpec = '%*s%s%s%s%s%s%s%s%[^\n\r]';
            fileID = fopen(filename, 'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
            timedataArray = textscan(fileID, timeformatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            fclose(fileID);
            fileID = fopen(filename, 'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
            otherdataArray = textscan(fileID, otherformatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            fclose(fileID);
            otherraw = repmat({''},length(otherdataArray{1}),length(otherdataArray)-1);
            for col=1:length(otherdataArray)-1
                otherraw(1:length(otherdataArray{col}),col) = otherdataArray{col};
            end
            othernumericData = NaN(size(otherdataArray{1},1),size(otherdataArray,2));
            for col=[1,2,3,4,5,6]
                % Converts strings in the input cell array to numbers. Replaced non-numeric
                % strings with NaN.
                otherrawData = otherdataArray{col};
                for row=1:size(otherrawData, 1);
                    % Create a regular expression to detect and remove non-numeric prefixes and
                    % suffixes.
                    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                    try
                        result = regexp(otherrawData{row}, regexstr, 'names');
                        numbers = result.numbers;

                        % Detected commas in non-thousand locations.
                        invalidThousandsSeparator = false;
                        if any(numbers==',');
                            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                            if isempty(regexp(numbers, thousandsRegExp, 'once'));
                                numbers = NaN;
                                invalidThousandsSeparator = true;
                            end
                        end
                        % Convert numeric strings to numbers.
                        if ~invalidThousandsSeparator;
                            numbers = textscan(strrep(numbers, ',', ''), '%f');
                            othernumericData(row, col) = numbers{1};
                            otherraw{row, col} = numbers{1};
                        end
                    catch me
                    end
                end
            end
            timeraw = repmat({''},length(timedataArray{1}),length(timedataArray)-1);
            for col=1:length(timedataArray)-1
                timeraw(1:length(timedataArray{col}),col) = timedataArray{col};
            end
            timenumericData = NaN(size(timedataArray{1},1),size(timedataArray,2));
            timerawData = timedataArray{1};
            for row=1:size(timerawData, 1);
                % Create a regular expression to detect and remove non-numeric prefixes and
                % suffixes.
                regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                try
                    result = regexp(rawData{row}, regexstr, 'names');
                    numbers = result.numbers;

                    % Detected commas in non-thousand locations.
                    invalidThousandsSeparator = false;
                    if any(numbers==',');
                        thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                        if isempty(regexp(numbers, thousandsRegExp, 'once'));
                            numbers = NaN;
                            invalidThousandsSeparator = true;
                        end
                    end
                    % Convert numeric strings to numbers.
                    if ~invalidThousandsSeparator;
                        numbers = textscan(strrep(numbers, ',', ''), '%f');
                        timenumericData(row, 1) = numbers{1};
                        timeraw{row, 1} = numbers{1};
                    end
                catch me
                end
            end
            time = str2double(timeraw(2:end));
            x = cell2mat(otherraw(2:end, 1));
            y = cell2mat(otherraw(2:end, 2));
            z = cell2mat(otherraw(2:end, 3));
            gx = cell2mat(otherraw(2:end, 4));
            gy = cell2mat(otherraw(2:end, 5));
            gz = cell2mat(otherraw(2:end, 6));
            class = otherraw(2:end, 7);
            obj = AWDataFile(time,x,y,z,gx,gy,gz, class);
        end 
    end
end