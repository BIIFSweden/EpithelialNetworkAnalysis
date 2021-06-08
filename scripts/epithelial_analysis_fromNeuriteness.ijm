
/*************** parameters ******************************/

type = "fromNeuriteness";
path = "/Users/gisele.miranda/Desktop/CIZ pics/sampleSNIs/";
maxProjPath = "/Users/gisele.miranda/Desktop/CIZ pics/max_projection/";
channel = "FITC";
threshold_method = "Otsu";
distance_threshold = 15;
correc_factor = 1.6;
area_threshold = 15000;
radiusDilation = 5;
nucleiSize = 1000;
pixPerMic = 3.07693;
useScale = true; // or 'false'

/*********************************************************/

function createLayer(ind) {
	roiManager("Select", ind);
	rName = Roi.getName();
	run("Create Mask");	
	if(rName == "A") rename("apical");
	else rename("basal");
}

function calculateDistance(refLayer, targetLayer, measures, areaValue) {
	selectWindow(refLayer);
	run("Invert");
	run("Exact Euclidean Distance Transform (3D)");
	rename("EDT-"+refLayer);
	
	selectWindow(targetLayer);
	run("Set Measurements...", "area mean standard min redirect=EDT-"+refLayer+" decimal=4");

	run("Analyze Particles...", "clear add");
	selectWindow("EDT-"+refLayer);
	run("Conversions...", "weighted");
	setOption("ScaleConversions", false);
	run("16-bit");
	run("From ROI Manager");

	// get number of pixels in reference layer
	selectWindow("EDT-"+refLayer);
	roiManager("Select", 0);
	roiManager("Measure");
	area = getResult("Area", 0);
	selectWindow("Results"); 
	run("Close");
	
	// get the other measures: average, std, min and max
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	selectWindow("EDT-"+refLayer);
	roiManager("Select", 0);
	roiManager("Measure");
	//area = getResult("Area", 0);

	if(useScale) {
		mean = getResult("Mean", 0)/pixPerMic;
		std = getResult("StdDev", 0)/pixPerMic;
		min = getResult("Min", 0)/pixPerMic;
		max = getResult("Max", 0)/pixPerMic;
	} else {
		mean = getResult("Mean", 0);
		std = getResult("StdDev", 0);
		min = getResult("Min", 0);
		max = getResult("Max", 0);
	}
	if(areaValue)
		measures = measures + area + ";" + mean + ";" + std + ";" + min + ";" + max;
	else
		measures = measures + mean + ";" + std + ";" + min + ";" + max;
	selectWindow(refLayer);
	run("Invert");

	selectWindow("Results"); 
	run("Close");
	selectWindow("ROI Manager"); 
	run("Close");

	return measures;
}

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

function saveOverlay(img, filePath, auto_adjust) {
	selectWindow(img);
	run("Duplicate...", "temp");
	rename("temp");	
	if(auto_adjust) autoAdjust();
	run("Restore Selection");
	run("Capture Image");
	saveAs("Tif", filePath);
	close("temp");
}

function findBiggestROI() {
	biggestInd = 0; // index of the first biggest ROI
	biggestROI = 0;
	
	for (i=0; i<nResults; i++) {
		area_i = getResult("Area", i);
		if(area_i > biggestROI) {
			biggestROI = area_i;
			biggestInd = i;
		}
	}

	return biggestInd;
}

function euclideanDist(x0,y0,x1,y1) {
	return Math.sqrt(Math.pow((x1-x0), 2) + Math.pow(y1-y0, 2));
}

function getEpithelium(roiName, wid, hei) {
	newImage("Temp", "8-bit black", wid, hei, 1);
	//roiName = replace(file, "/", "");
	roiManager("Open", maxProjPath + roiName + "_AB.zip");
	roiManager("Show All");
	
	roiManager("Select", 0);
	Roi.getCoordinates(x1, y1);
	run("Create Mask");								
	selectWindow("Temp");
	roiManager("Select", 1);				
	Roi.getCoordinates(x2, y2);
	run("Create Mask");

	// check the distance betweeen the ending points of each line and connect the points with minimum distance
	dA = euclideanDist(x1[0],y1[0],x2[0],y2[0]);
	dB = euclideanDist(x1[0],y1[0],x2[(x2.length-1)],y2[(y2.length-1)]);
	
	selectWindow("Mask");
	if(dA < dB) {
		makeLine(x1[0], y1[0], x2[0], y2[0]);
		run("Create Mask");
		makeLine(x1[(x1.length-1)], y1[(y1.length-1)], x2[(x2.length-1)], y2[(y2.length-1)]);
		run("Create Mask");
	} else {
		makeLine(x1[0], y1[0], x2[(x2.length-1)], y2[(y2.length-1)]);
		run("Create Mask");
		makeLine(x1[(x1.length-1)], y1[(y1.length-1)], x2[0], y2[0]);
		run("Create Mask");
	}
	run("Fill Holes");
	run("Select None");

	selectWindow("ROI Manager"); 
	run("Close");
	close("Temp");
}

/*************** main ************************/
dir = getFileList(path);
bufferMeasures = "";
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
overlayPath = path + "overlay_" + threshold_method + "_correc-factor-" + correc_factor + "_dist-" + distance_threshold + "_area-" + area_threshold + "_radDilation-" + radiusDilation + "_type-" + type + "_channel-" + channel + "_" + year + "-" + month + "-" + dayOfMonth + "/";
if(!File.exists(overlayPath)) File.makeDirectory(overlayPath);

for(cont=0; cont<dir.length; cont++) { 
	file = dir[cont];
	images = getFileList(path+file);
	dapiImg = false;
	channelImg = false;
	
	if(startsWith(file, "SNI") || startsWith(file, "Neg SNI")) { // exclude the overlay folder from the analysis if it is already created
		outpath = path + file + "output_" + threshold_method + "_correc-factor-" + correc_factor + "_dist-" + distance_threshold + "_area-" + area_threshold + "_radDilation-" + radiusDilation + "_type-" + type + "_channel-" + channel + "_" + year + "-" + month + "-" + dayOfMonth + "/";
		
		if(!File.exists(outpath)) File.makeDirectory(outpath);
	
		// open all channels and search for the one of interest
		for(cont2=0; cont2<images.length; cont2++) {
			image = images[cont2];
			
			if(matches(image, ".*" + channel + ".*") && startsWith(image, "Ne_")) { // check if file is an image
				print(image);
				folderImg = replace(file, "/", "");
				//aux = split(image,".");
				aux2 = substring(image, 3, lengthOf(image)); //split(image,"Ne_");

				// check if corresponding file exists in 'newTiffImages' folder
				if(File.exists(path + file + "newTiffImages/" + aux2)) {
					//run("Bio-Formats Importer", "open=[" + maxProjPath + folderImg + ".tif] color_mode=Default view=Hyperstack stack_order=XYCZT");
					run("Bio-Formats Importer", "open=[" + path + file + "newTiffImages/" + aux2 + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
					rename("inputImg");
					
					getPixelSize(unit, pw, ph, pd);
					w = getWidth();
					h = getHeight();

					// check if corresponding file exists in max projection folder
					if(File.exists(maxProjPath + folderImg + "_AB.zip")) {
						roiManager("Open", maxProjPath + folderImg + "_AB.zip");
						roiManager("Show All");
				
						countRois = roiManager("count");
						
						if(countRois == 2) { // then create contours of apical and basal layers and measure average height
							createLayer(0);
							//selectWindow(folderImg + ".tif");
							selectWindow("inputImg");
							createLayer(1);
							selectWindow("ROI Manager"); 
							run("Close");
							
							bufferMeasures = bufferMeasures + folderImg + ".tif;";
							bufferMeasures = calculateDistance("apical", "basal", bufferMeasures, true);
							bufferMeasures = bufferMeasures + ";";
							bufferMeasures = calculateDistance("basal", "apical", bufferMeasures, true);
							//bufferMeasures = bufferMeasures + "\n";
							run("Close All");

							// read original FITC channel again and perform segmentation with the distance transform
							run("Bio-Formats Importer", "open=[" + path + file + "newTiffImages/" + aux2 + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
							rename("original");
							wid = getWidth();
							hei = getHeight();
							run("Duplicate...", "mask");
							rename("mask");
							setThreshold(100, 65535); // rough threshold just to separate the background
							setOption("BlackBackground", true);
							run("Convert to Mask");
							run("Create Selection");
							run("Enlarge...", "enlarge=-8"); // cut some pixels of the border
			
							// read corresponding neuriteness image
							run("Bio-Formats Importer", "open=[" + path + file + image + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
							rename("neuriteness_original");
							run("Duplicate...", "neuriteness");
							rename("neuriteness");
							run("Median...", "radius=4");
							run("Restore Selection");
							run("Clear Outside");

							// generate epithelium contour + apical and basal layers to crop the epithelium region
							roiName = replace(file, "/", "");
							getEpithelium(roiName, wid, hei);
							
							selectWindow("mask");
							run("Close");
			
							selectWindow("Mask");
							rename("Mask-epithelium");
							saveAs("Tif", outpath + "mask_epithelium.tif");
							rename("Mask-epithelium");
							run("Create Selection");
			
							selectWindow("original");
							saveOverlay("original", outpath + "epithelium_overlay", true);
							close("epithelium_overlay.tif");
			
							// crop polygon mask on the neuriteness image, apply threshold and get threshold value
							selectWindow("neuriteness");
							run("Select None");
							run("Duplicate...", "mask_neurit");
							rename("mask_neurit");
			
							// save epithelium area
							selectWindow("Mask-epithelium");
							if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
							run("Set Measurements...", "area mean standard min redirect=None decimal=4");
							run("Analyze Particles...", "clear summarize");
							Table.rename("Summary", "Results");
							tot_area_epithelium = getResult("Total Area", 0);
							//print("total area epithelium: " + tot_area_epithelium);
							selectWindow("Results"); 
							run("Close");
			
							selectWindow("Mask-epithelium");
							run("Create Selection");
							selectWindow("mask_neurit");
							run("Restore Selection");
							run("Clear Outside");
							
							setAutoThreshold(threshold_method+" dark");
							getThreshold(lower, upper);
							run("Convert to Mask");
							close("mask_neurit");
			
							// apply correction factor to "lower" threshold
							selectWindow("neuriteness");
							run("Restore Selection");
							run("Clear Outside");
							setThreshold((lower*correc_factor), 255);
							setOption("BlackBackground", true);
							run("Convert to Mask");
							saveAs("Tif", overlayPath + "thresholded_neuriteness_" + image);
							saveAs("Tif", outpath + "thresholded_neuriteness.tif");
			
							// calculate euclidean distance transform on the thresholded image
							run("Invert");
							run("Exact Euclidean Distance Transform (3D)");
							rename("EDT");
							run("Conversions...", "weighted");
							setOption("ScaleConversions", false);
							run("16-bit");
			
							setThreshold(0, distance_threshold);
							run("Convert to Mask");
							
							// generate a mask for the apical layer to verify intersection with segmented regions
							newImage("Temp", "8-bit black", wid, hei, 1);
							roiManager("Open", maxProjPath + roiName + "_AB.zip");
							roiManager("Show All");
							
							roiManager("select", 0);
							rName = Roi.getName();
							if(!matches(rName, "A")) roiManager("select", 1);
							run("Create Mask");
							rename("mask_apical");
							selectWindow("ROI Manager"); 
							run("Close");
							close("Temp");
			
							// select only regions with area greater than area_threshold
							selectWindow("EDT");
							run("Analyze Particles...", "size="+area_threshold+"-Infinity show=Masks clear add");
							rename("Mask");
							setOption("BlackBackground", true);
							run("Convert to Mask");
							countSelecArea = roiManager("count");
							
							selectWindow("ROI Manager"); 
							run("Close");
			
							if(countSelecArea > 0) {
								// check intersections with the apical layer
								selectWindow("thresholded_neuriteness.tif");
								run("Invert");
								selectWindow("Mask");
								run("Select None");
								run("Morphological Filters", "operation=Dilation element=Square radius="+radiusDilation);
								//run("Fill Holes");
								rename("intact_net");
								close("Mask");
				
								run("Analyze Particles...", "size=0-Infinity clear add");
								nRois = roiManager("count");
								i = 0;
								indexToDelete = newArray(nRois);
								for(i=0; i<nRois; i++) {
									roiManager("Select", i);
									run("Create Mask");
									imageCalculator("AND create", "Mask","mask_apical");
									selectWindow("Result of Mask");
									run("Analyze Particles...", "clear summarize");
									Table.rename("Summary", "Results");
									tot_area = getResult("Total Area", 0);
									if(tot_area != 0) {
										indexToDelete[i] = 1;
										//run("Fill Holes");
									}
									else indexToDelete[i] = 0;
									close("Mask");
									close("Result of Mask");
								}
								
								// delete indexToDelete, the remaining selections are larger the threshold
								countSelected = 0;
								for(i=0; i<nRois; i++) {
									if(indexToDelete[i] == 0) {
										countSelected = countSelected + 1;
										roiManager("Select", i);
										run("Create Mask");
									}
								}
				
								selectWindow("ROI Manager"); 
								run("Close");
								selectWindow("Results"); 
								run("Close");
								
								if(countSelected != 0) {
									channelImg = true;
									
									close("intact_net");
									selectWindow("Mask");
									rename("intact_net");
									run("Select None");
		
									if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
									run("Analyze Particles...", "clear summarize");
									Table.rename("Summary", "Results");
									tot_area_segmented_net = getResult("Total Area", 0);
									//print("total_area_segmented_net: " + tot_area_segmented_net);
									selectWindow("Results"); 
									run("Close");
									
									run("Select None");
									run("Create Selection");
					
									selectWindow("thresholded_neuriteness.tif");
									run("Duplicate...", "segmented_net");
									rename("segmented_net");
									saveOverlay("segmented_net", outpath + "segmented_thresholded_neuriteness", false);
									
									close("segmented_network_" + image);
									selectWindow("segmented_net");
									run("Restore Selection");
									run("Clear Outside");
									run("Select None");
									if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
									else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
									run("Analyze Particles...", "clear summarize");
									Table.rename("Summary", "Results");
									tot_area_intact_net = getResult("Total Area", 0);
									//print("total_area_intact_net: " + tot_area_intact_net);
									selectWindow("Results"); 
									run("Close");
					
									saveOverlay("neuriteness_original", outpath + "neuriteness_distance_transform_overlay", true);
									saveAs("Tif", overlayPath + "neuriteness_distance_transform_" + image);
									saveOverlay("original", outpath + "original_distance_transform_overlay_original", true);
									saveAs("Tif", overlayPath + "original_distance_transform_" + image);

									// save intact_net
									selectWindow("intact_net");
									saveAs("Tif", outpath + "marker_of_interest_" + channel + "_channel.tif");
									run("Select None");

									//save in the buffer the mean intensity of the segmented neuriteness
									selectWindow("original");
									run("Restore Selection");
									run("Measure");
									orig_mean = getResult("Mean", 0);
									orig_std = getResult("StdDev", 0);
									selectWindow("Results"); 
									run("Close");

									//save in the buffer the mean intensity of the segmented neuriteness
									selectWindow("neuriteness_original");
									run("Restore Selection");
									run("Measure");
									neur_mean = getResult("Mean", 0);
									neur_std = getResult("StdDev", 0);
									selectWindow("Results"); 
									run("Close");
									
									bufferMeasures = bufferMeasures + ";" + tot_area_epithelium + ";" + tot_area_segmented_net + ";" + orig_mean + ";" + orig_std + ";" + tot_area_intact_net + ";" + neur_mean + ";" + neur_std + ";" + (tot_area_segmented_net/tot_area_epithelium) + ";" + (tot_area_intact_net/tot_area_epithelium) + ";" + (tot_area_intact_net/tot_area_segmented_net) + ";";

								} else {
									print("No regions were selected for this SNI");
									bufferMeasures = bufferMeasures + ";" + tot_area_epithelium + "\n";
								}
							} else {
								print("No regions with area bigger than the area threshold were found");
								bufferMeasures = bufferMeasures + ";" + tot_area_epithelium + "\n";
							}
							run("Close All");
							
						} else {
							run("Close All");
							print("The number of input layers is different than 2");
							bufferMeasures = bufferMeasures + folderImg + "\n";
							selectWindow("ROI Manager"); 
							run("Close");
						}
					} else {
						run("Close All");
						print("File " + maxProjPath + folderImg + "_AB.zip does not exist");
						bufferMeasures = bufferMeasures + folderImg + "\n";
					}
					run("Close All");
				} else {
					run("Close All");
					print(path + file + "newTiffImages/" + aux2 + " does not exist");
					bufferMeasures = bufferMeasures + folderImg + "\n";
				}
			} else if(matches(image, ".*DAPI.*")) { // check if file is an image
				dapiImg = true;

				run("Bio-Formats Importer", "open=[" + path + file + image + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
				rename("nuclei");

				w = getWidth();
				h = getHeight();

				getEpithelium(replace(file, "/", ""), w, h);
				selectWindow("Mask");
				rename("Mask-epithelium");
				run("Create Selection");

				selectWindow("nuclei");
				setAutoThreshold("Otsu dark");
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Morphological Filters", "operation=Dilation element=Square radius=6");
				run("Set Measurements...", "area mean standard min redirect=None decimal=4");
				run("Analyze Particles...", "size="+nucleiSize+"-Infinity show=Masks display clear add");
				selectWindow("Mask of nuclei-Dilation");
				run("Convert to Mask");

				run("Restore Selection");
				setBackgroundColor(0, 0, 0);
				run("Clear Outside");
				run("Select None");
				
				saveAs("Tif", outpath + "mask_nuclei.tif");
				selectWindow("ROI Manager"); 
				run("Close");
				selectWindow("Results"); 
				run("Close");
				run("Close All");
			}
		}

		if(dapiImg && channelImg) {
			run("Bio-Formats Importer", "open=[" + outpath + "mask_epithelium.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			rename("mask_epithelium");
			run("Create Selection");
				
			run("Bio-Formats Importer", "open=[" + outpath + "marker_of_interest_" + channel + "_channel.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			rename("mask_channel");
			run("Bio-Formats Importer", "open=[" + outpath + "mask_nuclei.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			rename("mask_nuclei");
			run("Restore Selection");
			setBackgroundColor(0, 0, 0);
			run("Clear Outside");
			run("Select None");
			saveAs("Tif", outpath + "mask_nuclei.tif");
			rename("mask_nuclei");
				
			run("Duplicate...", "temp_ROIs");
			rename("temp_ROIs");

			roiManager("Open", maxProjPath + folderImg + "_AB.zip");
			roiManager("Show All");
			
			createLayer(0);
			selectWindow("temp_ROIs");
			createLayer(1);
			selectWindow("ROI Manager"); 
			run("Close");
			close("temp_ROIs");

			// get area of parabasal layer
			imageCalculator("Add create", "mask_channel","mask_nuclei");
			selectWindow("Result of mask_channel");
			rename("parabasal_marker");
			run("Fill Holes");
			
			selectWindow("parabasal_marker");
			if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
			else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
			run("Set Measurements...", "area mean standard min redirect=None decimal=4");
			run("Analyze Particles...", "clear summarize");
			Table.rename("Summary", "Results");
			tot_area_parabasal = getResult("Total Area", 0);
			selectWindow("Results"); 
			run("Close");
			
			// get area of upper layer (between parabasal and superficial layer)
			imageCalculator("Subtract create", "mask_epithelium","parabasal_marker");
			selectWindow("Result of mask_epithelium");
			run("Analyze Particles...", "pixel show=Nothing display clear add");
			big_ind = findBiggestROI();
			roiManager("select", big_ind);
			run("Create Mask");
			rename("upper");
			selectWindow("Results"); 
			run("Close");
			selectWindow("ROI Manager"); 
			run("Close");
			close("Result of mask_epithelium");
			selectWindow("upper");
			
			if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
			else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
			run("Analyze Particles...", "clear summarize");
			Table.rename("Summary", "Results");
			tot_area_upper = getResult("Total Area", 0);
			selectWindow("Results"); 
			run("Close");
			
			// get area of lower layer (between upper + segmented marker and basal layer)
			selectWindow("mask_channel");
			imageCalculator("Add create", "upper","mask_channel");
			run("Fill Holes");

			imageCalculator("Subtract create", "mask_epithelium","Result of upper");
			rename("lower");
			if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
			else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
			run("Analyze Particles...", "clear summarize");
			Table.rename("Summary", "Results");
			tot_area_lower = getResult("Total Area", 0);
			selectWindow("Results"); 
			run("Close");

			// save layers
			selectWindow("parabasal_marker");
			saveAs("Tif", outpath + "parabasal_layer_marker_of_interest.tif");
			selectWindow("upper");
			saveAs("Tif", outpath + "upper_layer.tif");
			selectWindow("lower");
			saveAs("Tif", outpath + "lower_layer.tif");

			// add to buffer
			bufferMeasures = bufferMeasures + tot_area_parabasal + ";" + tot_area_upper + ";" + tot_area_lower + ";";

			// distance transform of segmented channel marker to apical layer
			bufferMeasures = calculateDistance("mask_channel", "apical", bufferMeasures, false);
			bufferMeasures = bufferMeasures + ";";
			close("EDT-mask_channel");
			selectWindow("apical");
			run("Remove Overlay");
			
			// distance transform of segmented channel marker to basal layer
			bufferMeasures = calculateDistance("mask_channel", "basal", bufferMeasures, false);
			bufferMeasures = bufferMeasures + ";";
			close("EDT-mask_channel");
			selectWindow("basal");
			run("Remove Overlay");

			bufferMeasures = bufferMeasures + "\n";
			//print(bufferMeasures);
			run("Close All");
			
		} else {
			run("Close All");
			if(!dapiImg && channelImg) bufferMeasures = bufferMeasures + "\n";
		}
	}
}
/*********************************************************/

// save summary files with measures
summaryPath = path + "summary_" + threshold_method + "_correc-factor-" + correc_factor + "_dist-" + distance_threshold + "_area-" + area_threshold + "_radDilation-" + radiusDilation + ".csv";
if(File.exists(summaryPath))
	File.delete(summaryPath);

// writing header and buffer in the summary file 
summaryFile = File.open(summaryPath);
print(summaryFile, "image;number of pixels in the basal layer;average height basal layer;std basal layer;min basal layer;max basal layer;number of pixels in the apical layer;height apical layer;std apical layer;min apical layer;max apical layer; area - epithelium; area - marker of interest (" + channel + " channel);MFI - marker of interest (" + channel + " channel);stdMFI - marker of interest (" + channel + " channel);area - neuriteness (" + channel + " channel);MFI - neuriteness (" + channel + " channel);stdMFI - neuriteness (" + channel + " channel);relative area - marker of interest (" + channel + " channel) per total epithelium;realtive area - neuriteness (" + channel + " channel) per total epithelium;realtive area - neuriteness (" + channel + " channel) per marker of interest (" + channel + " channel);area - marker of interest + parabasal layer (" + channel + " channel);area - upper layer (" + channel + " channel);area - lower layer (" + channel + " channel);average height between marker of interest (" + channel + " channel) and apical layer;std height between marker of interest (" + channel + " channel) and apical layer;min height between marker of interest (" + channel + " channel) and apical layer;max height between marker of interest (" + channel + " channel) and apical layer;average height between marker of interest (" + channel + " channel) and basal layer;std height between marker of interest (" + channel + " channel) and basal layer;min height between marker of interest (" + channel + " channel) and basal layer;max height between marker of interest (" + channel + " channel) and basal layer;\n");
print(summaryFile, bufferMeasures);
File.close(summaryFile);
