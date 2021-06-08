clear; clc

% Path
addpath('/Users/gisele.miranda/Box Sync/BIIF/Projects/AnnelieTjernlund2019-1/Gabriela/neuritenessscriptsboar/lib')
rad = 4; % radius of median filter

% selecting working directory
selpath = uigetdir;
d = dir(selpath);
nameFolds = {d.name}';
nameFolds(ismember(nameFolds,{'.','..','.DS_Store','analysis','summary'})) = [];

for countImg = 1:length(nameFolds) % for each image folder
    countImg
    imageFolder = strcat(selpath,'/',nameFolds{countImg});
    
    if ~contains(nameFolds{countImg},'analysis') && ~contains(nameFolds{countImg},'overlay') && ~contains(nameFolds{countImg},'summary')
        nameFolds{countImg}
        
        channels = dir(imageFolder);
        channelNames = {channels.name}';
        channelNames(ismember(channelNames,{'.','..','.DS_Store','newTiffImages'})) = [];

        for countChannel = 1:length(channelNames) % for each channel
            if endsWith(channelNames{countChannel}, '.tif') && ~contains(channelNames{countChannel},'DAPI') && ~contains(channelNames{countChannel},'Ne_')
                channel = strcat(selpath,'/',nameFolds{countImg},'/',channelNames{countChannel});
                im = imread(channel); % load images
                im = medfilt2(im,[rad rad]); % apply median filter
    
                sigma = 4;
                [imf,~,~] = NeuriteneesFilter2D(im,sigma);
                aux = split(channelNames{countChannel}, '_');
                aux = aux(2:end);
                aux = join(aux,'_');
    
                outName = strcat(selpath,'/',nameFolds{countImg},'/Ne_',nameFolds{countImg},'_',aux);
                outName = [outName{:}];
                imwrite(imf, outName);
            end
        end
    end

end
