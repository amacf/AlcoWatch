function processFileAppendExisting(existingArffFile)
    fileToAppend = AWDataFile.AWDataFileFromFile;
    segs = fileToAppend.segmentsWithSize(400);
    for i = 1:length(segs)
        segs(i) = segs(i).removeOutliers(1);
        segs(i).writeDataToArff(existingArffFile);
    end
end