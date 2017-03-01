function computeNormalization()
[name, path] = uigetfile({'*.arff'});
filename = fullfile(path,name);
fileID = fopen(filename, 'r');

outFilename = strcat('normalized_',name);
outFileID = fopen(outFilename, 'w');

delimiter = ',';
startRow = 62;

%Skip first 64 rows
header = textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% Print skipped part to output file
for i = 1:startRow - 1
    fprintf(outFileID,'%s\n', header{1}{i});
end
%Read data rows
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
instances = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);

numInstances = length(instances{1});
numFeaturesNormalize = 53;

i = 1;
while(i<numInstances)
    featureSums = zeros(59,1);
    numZeroInstances = 0;
    while(strcmp(instances{59}(i),'no_goggles'))
        numZeroInstances = numZeroInstances + 1;
        for j=1:numFeaturesNormalize
            thisFeature = instances{j}{i};
            if strcmp(thisFeature,'?')
                continue
            end
            thisFeature = str2double(thisFeature);
            featureSums(j) = featureSums(j) + thisFeature;
        end
        i = i + 1;
    end
    
    featureSums = featureSums ./ numZeroInstances;
    
    for j = 1:numZeroInstances
        currentInstanceI = i - numZeroInstances - 1 + j;
        for k = 1:numFeaturesNormalize
            nfeature = str2double(instances{k}{currentInstanceI}) / featureSums(k);
            if isnan(nfeature)
                fprintf(outFileID,'?,');
            else
                fprintf(outFileID,'%f,',nfeature);
            end
        end
        fprintf(outFileID,'%s,',instances{54}{currentInstanceI});
        fprintf(outFileID,'%s,',instances{55}{currentInstanceI});
        fprintf(outFileID,'%s,',instances{56}{currentInstanceI});
        fprintf(outFileID,'%s,',instances{57}{currentInstanceI});
        fprintf(outFileID,'%s,',instances{58}{currentInstanceI});
        fprintf(outFileID,'no_goggles\n');
    end
    
    while(i<=numInstances && (~strcmp(instances{59}(i),'no_goggles')))
        for k = 1:numFeaturesNormalize
            nfeature = str2double(instances{k}{i}) / featureSums(k);
            if isnan(nfeature)
                fprintf(outFileID,'?,');
            else
                fprintf(outFileID,'%f,',nfeature);
            end
        end
        fprintf(outFileID,'%s,',instances{54}{i});
        fprintf(outFileID,'%s,',instances{55}{i});
        fprintf(outFileID,'%s,',instances{56}{i});
        fprintf(outFileID,'%s,',instances{57}{i});
        fprintf(outFileID,'%s,',instances{58}{i});
        fprintf(outFileID,'%s\n',instances{59}{i});
        i = i + 1;
    end
end

end
