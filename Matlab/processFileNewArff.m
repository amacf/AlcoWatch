function processFileNewArff(height, weight, age, gender, pants)
    fileToAppend = AWDataFile.AWDataFileFromFile;
    segs = fileToAppend.segmentsWithSize(500);
    for i = 1:length(segs)
        segs(i) = segs(i).removeOutliers(1).correctTimeSeconds;
    end
    newFile = segs(1).outputArff(height, weight, age, gender, pants);
    for i = 2:length(segs)
        segs(i).writeDataToArff(newFile, height, weight, age, gender, pants);
    end
end