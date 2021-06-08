
path = getDirectory("Choose a Directory");

/************** parameters **************/
radiusDilation = 4;
radiusOpening = 10;
radiusMedianFilter = 6;
radiusMeanFilter = 4;
threshold = "Huang";
useNucleiChannel = false; // false
/****************************************/

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
	//print(hmin, hmax);
	run("Apply LUT");
}

// file that contains the summary of the results
if(useNucleiChannel)
	summaryPath = path + "summary_radDil-" + radiusDilation + "_radOpen-" + radiusOpening + "_radMedian-" + radiusMedianFilter + "_radMean-" + radiusMeanFilter + "_thr-" + threshold + "_withNuclei.csv";
else
	summaryPath = path + "summary_radDil-" + radiusDilation + "_radOpen-" + radiusOpening + "_radMedian-" + radiusMedianFilter + "_radMean-" + radiusMeanFilter + "_thr-" + threshold + ".csv";

if(File.exists(summaryPath))
	File.delete(summaryPath);

// writing the header of the summary file
summaryFile = File.open(summaryPath);
print(summaryFile, "image height_basal height_apical area_basal area_apical std_basal std_apical min_basal min_apical max_basal max_apical \n");

// iterate over each file of the selected folder
dir = getFileList(path);
overlayPath = path + "overlay_radDil-" + radiusDilation + "_radOpen-" + radiusOpening + "_radMedian-" + radiusMedianFilter + "_radMean-" + radiusMeanFilter + "_thr-" + threshold + "/";
File.makeDirectory(overlayPath);
for(cont=0; cont<dir.length; cont++) {
	file = dir[cont];
	images = getFileList(path+file);
	print(path+file);
	process = false;

	print(file);
	if(!startsWith(file, "overlay_") && !startsWith(file, "summary_")) { // exclude the overlay folder from the analysis if it is already created
		newImgPath = path + file + "newTiffImages/";
		File.makeDirectory(newImgPath);

		// open all channels
		for(cont2=0; cont2<images.length; cont2++) {
			image = images[cont2];			
			aux = split(image,".");

			imgMarkerName = split(image,"_");
			temp = split(file,"/");

			if(!useNucleiChannel) {
				if(!matches(image, ".*DAPI.*") && matches(image, ".*.tif.*")) { // check if file is an image
				//if(matches(image, ".*.tif.*")) { // check if file is an image
					run("Bio-Formats Importer", "open=[" + path + file + image + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
					saveAs("Tiff", newImgPath + temp[0] + "_" + imgMarkerName[1] + "_" + imgMarkerName[2]);
					autoAdjust();
					process = true;
				} 
			} else{
				if(matches(image, ".*.tif.*")) { // check if file is an image
					run("Bio-Formats Importer", "open=[" + path + file + image + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
					saveAs("Tiff", newImgPath + temp[0] + "_" + imgMarkerName[1] + "_" + imgMarkerName[2]);
					autoAdjust();
					process = true;
				} 
			}
			
			
		}
	
		// if channels are validated then continue processing
		if(process) {
			if(useNucleiChannel)
				outPath = path + file + "analysis_radDil-" + radiusDilation + "_radOpen-" + radiusOpening + "_radMedian-" + radiusMedianFilter + "_radMean-" + radiusMeanFilter + "_thr-" + threshold + "_withNuclei/";
			else
				outPath = path + file + "analysis_radDil-" + radiusDilation + "_radOpen-" + radiusOpening + "_radMedian-" + radiusMedianFilter + "_radMean-" + radiusMeanFilter + "_thr-" + threshold + "/";	
			File.makeDirectory(outPath);
			
			run("Images to Stack", "name=Stack title=[] use");
			run("Z Project...", "projection=[Max Intensity]"); // Max Intensity projection
			
			rename("image");
			run("Duplicate...", "signal");
			rename("signal");
			run("Duplicate...", "outerParts");
			rename("outerParts");
			run("Duplicate...", "firstThres");
			rename("firstThres");
			close("Stack");

			// first threshold: select all pixels above 0
			selectWindow("firstThres");
			setThreshold(1, 65535);
			run("Convert to Mask");
			run("Fill Holes");
			run("Create Selection");

			selectWindow("image");
			w = getWidth();
			h = getHeight();

			run("Enhance Contrast", "saturated=0.35");
			run("Apply LUT");
			run("Median...", "radius="+radiusMedianFilter);

			// morphological operations
			run("Morphological Filters", "operation=Dilation element=Square radius="+radiusDilation);
			selectWindow("image-Dilation");
			run("Morphological Filters", "operation=Opening element=Square radius="+radiusOpening);
			selectWindow("image-Dilation-Opening");
			run("Enhance Contrast", "saturated=0.35");
			run("Apply LUT");
			run("Mean...", "radius="+radiusMeanFilter);
			
			selectWindow("image-Dilation-Opening");
			run("Restore Selection");

			// apply threshold
			setAutoThreshold(threshold+" dark");
			close("image");
			selectWindow("image-Dilation-Opening");
			rename("image");
			
			setOption("BlackBackground", true);
			run("Convert to Mask");
			run("Fill Holes");
			run("Analyze Particles...", "size=0-Infinity show=[Masks] clear record add");

			// find the biggest particle in the binary mask corresponding to the epithelium
			roiManager("Measure");
			biggestROI = 0; // area of the first biggest ROI
			biggestInd = 0; // index of the first biggest ROI

			for (i=0; i<nResults; i++) {
				area_i = getResult("Area", i);
				if(area_i > biggestROI) {
					biggestROI = area_i;
					biggestInd = i;
				}
			}

			roiManager("Select", biggestInd);
			run("Create Mask");
			rename("epithelium");

			selectWindow("ROI Manager"); 
			run("Close");

			selectWindow("Results"); 
			run("Close");
		
			selectWindow("epithelium");
			save(outPath+"mask-epithelium.tif");
		
			selectWindow("outerParts");
			setAutoThreshold("Default dark");
			setThreshold(211, 65535);
			run("Convert to Mask");

			// subtract the epithelium mask from the original image, the intersections will the apical and basal layers
			imageCalculator("Subtract create", "outerParts","epithelium");
			selectWindow("Result of outerParts");
		
			run("Analyze Particles...", "size=0-Infinity show=Masks clear add");
			selectWindow("Mask of Result of outerParts");
			run("Convert to Mask");

			countRois = roiManager("count");

			// create the contours of the apical and basal layers
			if(countRois > 1) {
				// finding the 2 biggest ROIs
				roiManager("Measure");
				firstROI = 0; // area of the first biggest ROI
				firstInd = 0; // index of the first biggest ROI
				secondROI = 0; // area of the second biggest ROI
				secondInd = 0; // index of the second biggest ROI

				for (i=0; i<nResults; i++) {
					area_i = getResult("Area", i);
					if(area_i > firstROI) {
						secondROI = firstROI;
						secondInd = firstInd;
						firstROI = area_i;
						firstInd = i;
					} else if(area_i > secondROI) {
						secondROI = area_i;
						secondInd = i;
					}
				}
				
				roiManager("Select", firstInd);
				run("Create Mask");
				run("Fill Holes");
				rename("mask1");
				selectWindow("Mask of Result of outerParts");
			
				roiManager("Select", secondInd);
				run("Create Mask");
				run("Fill Holes");
				rename("mask2");
			
				close("image");
				close("outerParts");
				close("Result of outerParts");
				close("Mask of Result of outerParts");
			
				selectWindow("ROI Manager"); 
				run("Close");

				selectWindow("Results"); 
				run("Close");
			
				selectWindow("epithelium");
				run("Dilate");
				    	
				imageCalculator("AND create", "epithelium","mask1");
				selectWindow("Result of epithelium");
				rename("contour1");
				save(outPath+"contour1.tif");
				selectWindow("epithelium");
				imageCalculator("AND create", "epithelium","mask2");
				selectWindow("Result of epithelium");
				rename("contour2");
				save(outPath+"contour2.tif");
			
				//selectWindow("contour1");
				newImage("blank", "8-bit white", w, h, 1);
				imageCalculator("Subtract create", "blank", "contour1");
				selectWindow("Result of blank");
				run("Exact Euclidean Distance Transform (3D)");
				selectWindow("EDT");
				rename("EDT-contour1");
				selectWindow("contour2");
				run("Set Measurements...", "area mean standard min perimeter redirect=None decimal=4");
				run("Analyze Particles...", "clear add");
				selectWindow("EDT-contour1");
				run("Conversions...", "weighted");
				setOption("ScaleConversions", false);
				run("16-bit");
				run("From ROI Manager");
				save(outPath+"EDT-contour1.tif");
				//run("Measure");
				roiManager("Measure");
				close("Result of blank");
				meanContour1 = getResult("Mean", 0);
				areaContour1 = getResult("Area", 0);
				stdContour1 = getResult("StdDev", 0);
				minContour1 = getResult("Min", 0);
				maxContour1 = getResult("Max", 0);

				selectWindow("Results"); 
				run("Close");
						
				imageCalculator("Subtract create", "blank", "contour2");
				selectWindow("Result of blank");
				run("Exact Euclidean Distance Transform (3D)");
				selectWindow("EDT");
				rename("EDT-contour2");
				selectWindow("contour1");
				run("Analyze Particles...", "clear add");
				selectWindow("EDT-contour2");
				run("Conversions...", "weighted");
				setOption("ScaleConversions", false);
				run("16-bit");
				run("From ROI Manager");
				save(outPath+"EDT-contour2.tif");
				//run("Measure");
				roiManager("Measure");
				close("Result of blank");
				meanContour2 = getResult("Mean", 0);
				areaContour2 = getResult("Area", 0);
				stdContour2 = getResult("StdDev", 0);
				minContour2 = getResult("Min", 0);
				maxContour2 = getResult("Max", 0);

				selectWindow("Results"); 
				run("Close");
				selectWindow("ROI Manager"); 
				run("Close");

				// save the measures in the summary file
				selectWindow("signal");
				img = replace(file, " ", "_");
				aux = split(file,"/");
				if(areaContour1 > areaContour2) {
					print(summaryFile, img + " " + meanContour1 + " " + meanContour2 + " " + 
												   areaContour1 + " " + areaContour2 + " " + 
												   stdContour1 + " " + stdContour2 + " " + 
												   minContour1 + " " + minContour2 + " " + 
												   maxContour1 + " " + maxContour2 + "\n");
					// save overlay
					selectWindow("EDT-contour1");
					run("To ROI Manager");
					roiManager("Set Color", "red");
					roiManager("Set Line Width", 5);
					selectWindow("signal");
					run("From ROI Manager");
					selectWindow("ROI Manager"); 
					run("Close");

					selectWindow("EDT-contour2");
					run("To ROI Manager");
					roiManager("Set Color", "yellow");
					roiManager("Set Line Width", 5);
					selectWindow("signal");
					run("From ROI Manager");

					save(outPath + "_overlay.tif");
					save(overlayPath + aux[0] + "_overlay.jpeg");
					
					selectWindow("ROI Manager"); 
					run("Close");
				} else {
					print(summaryFile, img + " " + meanContour2 + " " + meanContour1 + " " + 
												   areaContour2 + " " + areaContour1 + " " + 
												   stdContour2 + " " + stdContour1 + " " + 
												   minContour2 + " " + minContour1 + " " + 
												   maxContour2 + " " + maxContour1 + "\n");
				    // saving overlay
					selectWindow("EDT-contour1");
					run("To ROI Manager");
					roiManager("Set Color", "yellow");
					roiManager("Set Line Width", 5);
					selectWindow("signal");
					run("From ROI Manager");
					selectWindow("ROI Manager"); 
					run("Close");

					selectWindow("EDT-contour2");
					run("To ROI Manager");
					roiManager("Set Color", "red");
					roiManager("Set Line Width", 5);
					selectWindow("signal");
					run("From ROI Manager");
					
					save(outPath + "_overlay.tif");
					save(overlayPath + aux[0] + "_overlay.jpeg");
					
					selectWindow("ROI Manager"); 
					run("Close");
				}
				run("Close All");
				
			} else {
				run("Close All");
				selectWindow("ROI Manager"); 
				run("Close");
			}
		}
	}
}

selectWindow("Log"); 
run("Close");

File.close(summaryFile);