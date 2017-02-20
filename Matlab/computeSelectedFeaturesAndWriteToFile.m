function computeSelectedFeaturesAndWriteToFile()
fileID = fopen('input.txt','r');
inputJSON = JSON.parse(fscanf(fileID,'%s'));
seg = AWDataSegment(getDoubleList(inputJSON.accelerometer.t), ...
                    getDoubleList(inputJSON.accelerometer.x), ...
                    getDoubleList(inputJSON.accelerometer.y), ...
                    getDoubleList(inputJSON.accelerometer.z), ...
                    getDoubleList(inputJSON.gyroscope.x), ...
                    getDoubleList(inputJSON.gyroscope.y), ...
                    getDoubleList(inputJSON.gyroscope.z), ...
                    'unknown');
outString = seg.getSelectedFeatures;
fileID = fopen('resultingGeneratedFeatures.csv','w');
fprintf(fileID,outString);
end


function dlist = getDoubleList(str)
str = str(2:length(str)-1);
dlist = str2double(strsplit(str,','));
end