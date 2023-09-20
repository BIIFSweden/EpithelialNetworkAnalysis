/************** parameters **************/
useNucleiChannel = false;
/****************************************/

// source: https://github.com/miura/IJ_BCautoMacro/tree/master
function autoAdjust() {
	auto_thres = 5000;
	getRawStatistics(pixcount);
	limit = pixcount/10;
	threshold = pixcount/auto_thres;
	nBins = 256;
	getHistogram(values, histA, nBins);
	i = -1;
	found = false;
	do {
		counts = histA[++i];
		if (counts > limit) counts = 0;
		found = counts > threshold;
	} while ((!found) && (i < histA.length-1))
	hmin = values[i];
	
	i = histA.length;
	do {
		counts = histA[--i];
		if (counts > limit) counts = 0;
		found = counts > threshold;
	} while ((!found) && (i > 0))
	hmax = values[i];
	
	setMinAndMax(hmin, hmax);
	run("Apply LUT");
}

// iterate over each file of the selected folder
path = getDirectory("Choose a Directory");
pathMax = "";
print(path);
aux = split(path, "/");

for(i=0; i<aux.length-1; i++) {
	pathMax = pathMax + "/" + aux[i];
}

pathMax = pathMax + "/" + aux[aux.length-1] + "_MAX_projection/";
print(pathMax);

if(!File.isDirectory(pathMax))
	File.makeDirectory(pathMax);
	
dir = getFileList(path);

for(cont=0; cont<dir.length; cont++) {
	file = dir[cont];
	images = getFileList(path+file);
	print(path+file);
	process = false;
	
	if(!startsWith(file, "overlay_") && !startsWith(file, "summary_")) { // exclude the overlay folder from the analysis if it is already created
		newImgPath = path + file + "newTiffImages/";
		File.makeDirectory(newImgPath);

		// open all channels
		for(cont2=0; cont2<images.length; cont2++) {
			image = images[cont2];

			if(matches(image, ".*.tif.*") && !startsWith(image, "Ne_") && !startsWith(image, "output_") && !startsWith(image, "newTiff")) {
				aux = split(image,".");
	
				imgMarkerName = split(image,"_");
				temp = split(file,"/");
	
				if(!useNucleiChannel) {
					if(!matches(image, ".*DAPI.*")) { // check if file is an image
						run("Bio-Formats Importer", "open=[" + path + file + image + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
						saveAs("Tiff", newImgPath + temp[0] + "_" + imgMarkerName[1] + "_" + imgMarkerName[2]);
						autoAdjust();
						process = true;
					} 
				} else{
					run("Bio-Formats Importer", "open=[" + path + file + image + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
					saveAs("Tiff", newImgPath + temp[0] + "_" + imgMarkerName[1] + "_" + imgMarkerName[2]);
					autoAdjust();
					process = true;
				}
			}
		}
	
		// if channels are validated then continue processing
		if(process) {
			annotation = split(file, "/");
			outPath = pathMax + annotation[0];
			print(outPath);
			
			run("Images to Stack", "name=Stack title=[] use");
			run("Z Project...", "projection=[Max Intensity]"); // Max Intensity projection
			
			rename("image");
			saveAs("Tiff", outPath);
		}

		run("Close All");
	}
}