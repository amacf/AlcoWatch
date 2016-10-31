function showFileGCM()
    file = AWDataFile.AWDataFileFromFile;
    segs = file.segmentsWithSize(length(file.time));
    segs = segs.correctTimeSeconds;
    plot(segs.time,segs.gcmG);
end