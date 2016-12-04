function processFileAppendExisting(existingArffFile, height, weight, age, gender, pants)
    fileToAppend = AWDataFile.AWDataFileFromFile;
    segs = fileToAppend.segmentsWithSize(500);
    for i = 1:length(segs)
        segs(i) = segs(i).removeOutliers(1);
        segs(i).writeDataToArff(existingArffFile, height, weight, age, gender, pants);
    end
end