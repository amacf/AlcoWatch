function processFileNewArff()
    fileToAppend = AWDataFile.AWDataFileFromFile;
    segs = fileToAppend.segmentsWithSize(400);
    for i = 1:length(segs)
        segs(i) = segs(i).removeOutliers(1).correctTimeSeconds;
    end
    newFile = segs(1).outputArff;
    for i = 2:length(segs)
        segs(i).writeDataToArff(newFile);
    end
end