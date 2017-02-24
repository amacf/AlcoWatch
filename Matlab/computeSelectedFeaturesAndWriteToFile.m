function computeSelectedFeaturesAndWriteToFile()
try
    logID = fopen('matlabLog.txt','w');
    fprintf(logID,'Started\n');
    fileID = fopen('input.txt','r');
    inputJSON = fscanf(fileID,'%s');
    fclose(fileID);
    fprintf(logID,'Got input...looks like:\n');
    fprintf(logID,inputJSON);
    inputJSON = JSON.parse(inputJSON);
    seg = AWDataSegment(getDoubleList(inputJSON.accelerometer.t), ...
                        getDoubleList(inputJSON.accelerometer.x), ...
                        getDoubleList(inputJSON.accelerometer.y), ...
                        getDoubleList(inputJSON.accelerometer.z), ...
                        getDoubleList(inputJSON.gyroscope.x), ...
                        getDoubleList(inputJSON.gyroscope.y), ...
                        getDoubleList(inputJSON.gyroscope.z), ...
                        'unknown');
    fprintf(logID,'\nMade into seg\n');
    seg = seg.removeOutliers(1).correctTimeSeconds;
    fprintf(logID,'Removed Outliers\n');
    seg = seg.takeMovingAverage(10);
    fprintf(logID,'Performed Smoothing\n');
    outString = seg.getSelectedFeatures;
    fprintf(logID,'Calculated Features\n');
    fileID = fopen('resultingGeneratedFeatures.csv','w');
    fprintf(fileID,outString);
    fprintf(logID,'Printed to CSV\n');
    fclose(logID);
    fclose(fileID);
catch Exception
    fprintf(logID,getReport(Exception));
end
end


function dlist = getDoubleList(str)
str = str(2:length(str)-1);
dlist = str2double(strsplit(str,','));
end
