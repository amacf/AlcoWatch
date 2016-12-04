function m = processAllFiles(filesLocation)
    [name, path] = uigetfile({'*.csv'});
    filename = fullfile(path,name);
    files = dir(filesLocation);
    file1 = files(3).name;
    h = waitbar(0,'Processing Files...');
    fileToStart = AWDataFile.AWDataFileFromFile(fullfile(filesLocation,file1));
    segs = fileToStart.segmentsWithSize(500);
    for i = 1:length(segs)
        segs(i) = segs(i).removeOutliers(1).correctTimeSeconds;
    end
    [height,weight,age,gender,pants] = getInfoForID(filename,file1);
    newFile = segs(1).outputArff(height, weight, age, gender, pants);
    for i = 2:length(segs)
        segs(i).writeDataToArff(newFile, height, weight, age, gender, pants);
    end
    waitbar(1/(length(files)-2));
    for j = 4:length(files)
        thisFile = files(j).name;
        fileToAppend = AWDataFile.AWDataFileFromFile(fullfile(filesLocation,thisFile));
        [height,weight,age,gender,pants] = getInfoForID(filename,thisFile);
        segs = fileToAppend.segmentsWithSize(500);
        for i = 1:length(segs)
            segs(i) = segs(i).removeOutliers(1).correctTimeSeconds;
            segs(i).writeDataToArff(newFile, height, weight, age, gender, pants);
        end
        waitbar((j-2)/(length(files)-2));
    end
    waitbar(1);
end