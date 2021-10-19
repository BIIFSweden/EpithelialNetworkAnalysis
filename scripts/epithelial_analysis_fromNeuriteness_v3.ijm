
// ===============================
// AUTHOR : Gisele Miranda 
// IVESTIGATORS : Annelie Tjernlund, Mathias Franzén Boger, Gabriela Edfeldt, Kristina Broliden
// CREATE DATE : 2020 - 06 - 18
// PURPOSE : Structural analysis of the cervical epithelial tissue
// NOTES : Required plugins - MorphoLibJ (https://imagej.net/MorphoLibJ)
// ===============================


/*************** parameters ******************************/

path = "/Users/gisele.miranda/Desktop/Annelie/test_13_08_2021/images/"; // images to be analyzed
maxProjPath = "/Users/gisele.miranda/Desktop/Annelie/test_13_08_2021/projection/"; // path to the folder with the max projection files
neur_channels = newArray("FITC", "cy3"); // select channels to be combined for neuriteness
threshold_method = "Otsu";
distance_threshold = 15;
correc_factor = 1.6;
area_threshold = 50000;
radiusDilation = 5;
nucleiSize = 1000;
pixPerMic = 3.07693;
useScale = true; // or 'false'
combine_neur = false;

/*********************************************************/

dir = getFileList(path);
bufferMeasuresFITC = "";
bufferMeasuresCy3 = "";
bufferMeasuresCy5 = "";
bufferMeasuresMcherry = "";
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
exp_id = threshold_method + "_correc-factor-" + correc_factor + "_dist-" + distance_threshold + "_area-" + area_threshold + "_radDilation-" + radiusDilation + "_" + year + "-" + (month+1) + "-" + dayOfMonth + "_" + hour + "-" + minute;
if(combine_neur) exp_id = exp_id + "_combined";
overlayPath = path + "overlay_" + exp_id + "/";
if(!File.exists(overlayPath)) File.makeDirectory(overlayPath);

for(cont=0; cont<dir.length; cont++) { // for each SNI
	sni = dir[cont];
	images = getFileList(path+sni);

	if(startsWith(sni, "SNI") || startsWith(sni, "Neg SNI")) { // check if input file is a SNI
		outpath = path + sni + "output_" + exp_id + "/";
		if(!File.exists(outpath)) File.makeDirectory(outpath);
		
		// crop epithelium region
		print(sni);
		sni = replace(sni, "/", "");

		// check if files exist
		fitcChannel = File.exists(path + sni + "/newTiffImages/" + sni + "_FITC_Extended.tif");
		projOutlines = File.exists(maxProjPath + sni + "_AB.zip");

		if(fitcChannel && projOutlines) { // check if FITC channel and the corresponding max projection exists for the segmentation of the epithelium

			// get epithelium
			getEpithelium(sni, outpath);
			
			// get average epithelium height
			epithHeight = getEpitheliumHeight(sni, outpath);
	
			// get combined Neuriteness, otherwise get individual neuriteness below
			if(combine_neur) getCombinedNeuriteness(sni, outpath);
	
			// segment nuclei
			segmentNuclei(sni, outpath);
		
			// iterate through all channels and search for the one(s) of interest
			for(cont2=0; cont2<images.length; cont2++) {
				image = images[cont2];
				neurit_file = false;
				
				if(containsSNI(image) && startsWith(image, "Ne_")) {
					aux = split(image, "_");
					print("image: " + image);
					channel = aux[(aux.length-2)];
	
					// check if neuriteness file exists
					if(combine_neur) {
						if(File.exists(outpath + "")) {
							open(outpath + "thresholded_neuriteness_combined.tif");
							rename("thresholded_neuriteness");
							neurit_file = true;
						} else print("file " + outpath + "thresholded_neuriteness.tif does not exist");
					} else {
						getNeuriteness(sni, image, channel, outpath);
						open(outpath + "thresholded_neuriteness_" + channel + ".tif");
						rename("thresholded_neuriteness");
						neurit_file = true;
					}
	
					// if yes, then start processing
					if(neurit_file) {
						run("Invert");
						run("Exact Euclidean Distance Transform (3D)");
						rename("EDT");
						run("Conversions...", "weighted");
						setOption("ScaleConversions", false);
						run("16-bit");
		
						setThreshold(0, distance_threshold);
						run("Convert to Mask");
		
						// generate a mask for the apical layer to verify intersection with the segmented regions
						newImage("Temp", "8-bit black", getWidth(), getHeight(), 1);
						roiManager("Open", maxProjPath + sni + "_AB.zip");
						roiManager("Show All");
						
						roiManager("select", 0);
						if(!matches(Roi.getName(), "A")) roiManager("select", 1);
						run("Create Mask");
						rename("mask_apical");
						selectWindow("ROI Manager"); 
						run("Close");
		
						// select only regions with area greater than area_threshold
						selectWindow("EDT");
						run("Set Measurements...", "area mean standard min redirect=None decimal=4");
						run("Analyze Particles...", "size="+area_threshold+"-Infinity show=Masks clear add");
						rename("Mask");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						countSelecArea = roiManager("count");
						
						selectWindow("ROI Manager"); 
						run("Close");
						close("Temp");
		
						if(countSelecArea > 0) {
							// check intersections with the apical layer
							selectWindow("thresholded_neuriteness");
							run("Invert");
							selectWindow("Mask");
							run("Select None");
							run("Morphological Filters", "operation=Dilation element=Square radius="+radiusDilation);
							rename("intact_net");
							close("Mask");

							open(outpath + "mask_epithelium.tif");
							run("Create Selection");
							selectWindow("intact_net");
							run("Restore Selection");
							run("Clear Outside");
							run("Select None");
							close("mask_epithelium.tif");
			
							run("Analyze Particles...", "size=0-Infinity clear add");
							nRois = roiManager("count");
							a = nRois-1;
							indexToDelete = newArray(nRois);
							for(i=a; i>=0; i--) {
								roiManager("Select", i);
								run("Create Mask");
								roiManager("Measure");
								area_i = getResult("Area", 0);
								
								imageCalculator("AND create", "Mask","mask_apical");
								selectWindow("Result of Mask");
								run("Analyze Particles...", "clear summarize");
								Table.rename("Summary", "Results");
								tot_area = getResult("Total Area", 0);
								//if(tot_area != 0 && area_i < area_threshold) {
								if(area_i < area_threshold) {
									indexToDelete[i] = 1;
								}
								else indexToDelete[i] = 0;
								
								close("Mask");
								close("Result of Mask");
							}
							
							// delete indexToDelete, the remaining selections are larger the threshold
							countSelected = 0;
							for(i=0; i<nRois; i++) {
								if(indexToDelete[i] == 0) {
									countSelected++;
									roiManager("Select", i);
									run("Create Mask");
								}
							}
			
							selectWindow("ROI Manager"); 
							run("Close");
							selectWindow("Results"); 
							run("Close");
							run("Select None");
		
							if(countSelected != 0) {
		
								selectWindow("Mask");
			
								if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
								run("Analyze Particles...", "clear summarize");
								Table.rename("Summary", "Results");
								tot_area_segmented_net = getResult("Total Area", 0);
								selectWindow("Results"); 
								run("Close");
								
								if(combine_neur) saveAs("Tif", outpath + "marker_of_interest.tif");
								else saveAs("Tif", outpath + "marker_of_interest_" + channel + ".tif");
								rename("intact_net");
								run("Create Selection");
				
								selectWindow("thresholded_neuriteness");
								run("Restore Selection");
								run("Flatten");
								if(combine_neur) saveAs("Tif", outpath + "segmented_thresholded_neuriteness.tif");
								else saveAs("Tif", outpath + "segmented_thresholded_neuriteness_" + channel + ".tif");
								close();
		
								selectWindow("thresholded_neuriteness");
								close("\\Others");
		
								run("Clear Outside");
								if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
								else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
								run("Analyze Particles...", "clear summarize");
								Table.rename("Summary", "Results");
								area_segmented_net = getResult("Total Area", 0);
								selectWindow("Results"); 
								run("Close");
		
								// get MFI from original channel image
								run("Bio-Formats Importer", "open=[" + path + sni + "/newTiffImages/" + replace(image, "Ne_", "") + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
								rename("original");
								run("Restore Selection");
								
								run("Measure");
								orig_mean = getResult("Mean", 0);
								orig_std = getResult("StdDev", 0);
								selectWindow("Results"); 
								run("Close");
		
								selectWindow("original");
								run("Flatten");
								saveAs("Tif", outpath + "segmented_original_" + channel + ".tif");
								saveAs("Tif", overlayPath + "segmented_original_" + channel + "_" + image);
								close();
		
								// get MFI from original neuriteness image
								run("Bio-Formats Importer", "open=[" + path + sni + "/" + image + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
								rename("neuriteness");
								run("Restore Selection");
		
								run("Measure");
								neur_mean = getResult("Mean", 0);
								neur_std = getResult("StdDev", 0);
								selectWindow("Results"); 
								run("Close");
		
								selectWindow("neuriteness");
								run("Flatten");
								saveAs("Tif", outpath + "segmented_neuriteness_original_" + channel + ".tif");
								saveAs("Tif", overlayPath + "segmented_neuriteness_original_" + channel + "_" + image);
		
								run("Close All");

								// get parabasal, upper and lower layers

								run("Bio-Formats Importer", "open=[" + outpath + "mask_epithelium.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
								rename("mask_epithelium");
								if(combine_neur) run("Bio-Formats Importer", "open=[" + outpath + "marker_of_interest.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
								else run("Bio-Formats Importer", "open=[" + outpath + "marker_of_interest_" + channel + ".tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
								rename("mask_channel");
								run("Bio-Formats Importer", "open=[" + outpath + "mask_nuclei.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
								rename("mask_nuclei");
				
								run("Duplicate...", "temp_ROIs");
								rename("temp_ROIs");

								roiManager("Open", maxProjPath + sni + "_AB.zip");
								roiManager("Show All");
								createLayer(0, outpath);
								selectWindow("temp_ROIs");
								createLayer(1, outpath);
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
								if(combine_neur) {
									selectWindow("parabasal_marker");
									saveAs("Tif", outpath + "parabasal_layer.tif");
									selectWindow("upper");
									saveAs("Tif", outpath + "upper_layer.tif");
									selectWindow("lower");
									saveAs("Tif", outpath + "lower_layer.tif");
								} else {
									selectWindow("parabasal_marker");
									saveAs("Tif", outpath + "parabasal_layer_" + channel + ".tif");
									selectWindow("upper");
									saveAs("Tif", outpath + "upper_layer_" + channel + ".tif");
									selectWindow("lower");
									saveAs("Tif", outpath + "lower_layer_" + channel + ".tif");
								}

								// add to buffer
								bufferChannel = "" + tot_area_parabasal + ";" + tot_area_upper + ";" + tot_area_lower + ";";
					
								// distance transform of segmented channel marker to apical layer
								dist = calculateDistance("mask_channel", "apical", false);
								bufferChannel = bufferChannel + dist;
								close("EDT-mask_channel");
								selectWindow("apical");
								run("Remove Overlay");
								
								// distance transform of segmented channel marker to basal layer
								dist = calculateDistance("mask_channel", "basal", false);
								bufferChannel = bufferChannel + ";" + dist;
								close("EDT-mask_channel");
								selectWindow("basal");
								run("Remove Overlay");
								run("Close All");

								// run watershed and get intact net
								area_inteact_net = 0;
								if(combine_neur) area_inteact_net = getIntactNet("", outpath, overlayPath, image);
								else area_inteact_net = getIntactNet(channel, outpath, overlayPath, image);

								// update buffers
								tot_area_epithelium = getEpitheliumArea(sni, outpath);
								bufferMeasures = image + ";" + epithHeight + ";" + tot_area_epithelium + ";" + tot_area_segmented_net + ";" + orig_mean + ";" + orig_std + ";" + area_segmented_net + ";" + neur_mean + ";" + neur_std + ";" + (tot_area_segmented_net/tot_area_epithelium) + ";" + (area_segmented_net/tot_area_epithelium) + ";" + (area_segmented_net/tot_area_segmented_net) + ";" + bufferChannel + ";" + area_inteact_net + ";" + (area_inteact_net/area_segmented_net) + "\n";
								if(matches(image, ".*FITC.*"))
									bufferMeasuresFITC = bufferMeasuresFITC + bufferMeasures;
								else if(matches(image, ".*cy3.*"))
									bufferMeasuresCy3 = bufferMeasuresCy3 + bufferMeasures;
								else if(matches(image, ".*Cy5.*"))
									bufferMeasuresCy5 = bufferMeasuresCy5 + bufferMeasures;
								else if(matches(image, ".*m cherry.*"))
									bufferMeasuresMcherry = bufferMeasuresMcherry + bufferMeasures;
							} else {
								print("No regions were selected for this SNI");
								run("Close All");
							}
						} else {
							print("No regions with area bigger than the area threshold were found");
							run("Close All");
						}
					} else run("Close All");
				}
			} 
		} else {
			if(!fitcChannel) print("File: " + path + sni + "/newTiffImages/" + sni + "_FITC_Extended.tif does not exist");
			if(!projOutlines) print("File: " + maxProjPath + sni + "_AB.zip does not exist");
		}
	}
}

// save summary files with measures 
for(countChannels=0; countChannels<neur_channels.length; countChannels++) {
	ch = neur_channels[countChannels];

	summaryPath = path + "summary" + ch + "_" + exp_id + ".csv";

	if(File.exists(summaryPath))
		File.delete(summaryPath);

	summaryFile = File.open(summaryPath);
	print(summaryFile, "image;number of pixels in the basal layer;average height basal layer;std basal layer;min basal layer;max basal layer;number of pixels in the apical layer;height apical layer;std apical layer;min apical layer;max apical layer;area - epithelium; area - marker of interest;MFI - marker of interest;stdMFI - marker of interest;area - neuriteness;MFI - neuriteness;stdMFI - neuriteness;relative area - marker of interest per total epithelium;realtive area - neuriteness per total epithelium;realtive area - neuriteness per marker of interest;area - marker of interest + parabasal layer;area - upper layer;area - lower layer;average height between marker of interest and apical layer;std height between marker of interest and apical layer;min height between marker of interest and apical layer;max height between marker of interest and apical layer;average height between marker of interest and basal layer;std height between marker of interest and basal layer;min height between marker of interest and basal layer;max height between marker of interest and basal layer;area - intact net;relative area - intact net per neuriteness total area\n");
	print(summaryFile, bufferMeasures);
	File.close(summaryFile);
}

/****************** functions ****************************/

function getIntactNet(channel, out_path, overlay_path, img) {

	// open neuriteness image
	if(matches(channel, ""))
		open(out_path + "thresholded_neuriteness_combined.tif");
	else
		open(out_path + "thresholded_neuriteness_" + channel + ".tif");
		
	run("Median...", "radius=4");
	rename("inputImg");
	run("Duplicate...", "markerImg");
	rename("markerImg");
	run("Invert");
	
	run("Marker-controlled Watershed", "input=inputImg marker=markerImg binary calculate use");
	run("16-bit");
	setThreshold(2, 65535);
	run("Convert to Mask");
	saveAs("Tif", out_path + "watershed_" + channel + ".tif");
	rename("watershed");

	run("Morphological Filters", "operation=Dilation element=Square radius=6");
	imageCalculator("AND create", "inputImg","watershed-Dilation");
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	run("Set Measurements...", "area mean standard min redirect=None decimal=4");
	selectWindow("Result of inputImg");
	run("Analyze Particles...", "clear summarize"); // get area of the intact net
	Table.rename("Summary", "Results");
	tot_area = getResult("Total Area", 0);

	selectWindow("Results"); 
	run("Close");

	run("Morphological Reconstruction", "marker=[Result of inputImg] mask=inputImg type=[By Dilation] connectivity=8");
	run("Create Selection");
	
	selectWindow("inputImg");
	run("Restore Selection");
	run("Properties... ", "  width=3");
	run("Flatten");
	saveAs("Tif", out_path + "intact_net_" + channel + ".tif");
	saveAs("Tif", overlay_path + "intact_net" + channel + "_" + img);

	run("Close All");

	return(tot_area);
}

function getEpithelium(sni, out_path) {

	// read original FITC channel and retrieve the binary mask of the tissue, excluding the black background
	run("Bio-Formats Importer", "open=[" + path + sni + "/newTiffImages/" + sni + "_FITC_Extended.tif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename("mask_tissue");

	// rough threshold just to separate the background
	if(bitDepth() == 8)
		setThreshold(5, 256);
	else 
		setThreshold(5, 65535);
	
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Create Selection");
	run("Enlarge...", "enlarge=-3"); // cut some pixels of the border
	run("Clear Outside");
	run("Select None");
	
	// create image that will be the binary mask of the epithelium
	newImage("mask_epithelium", "8-bit black", getWidth(), getHeight(), 1);
	roiManager("Open", maxProjPath + sni + "_AB.zip");
	roiManager("Show All");
	
	roiManager("Select", 0);
	Roi.getCoordinates(x1, y1);
	run("Create Mask");
	selectWindow("mask_epithelium");
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

	// final crop using both masks: tissue + epithelium
	selectWindow("mask_tissue");
	run("Create Selection");
	selectWindow("Mask");
	run("Restore Selection");
	run("Clear Outside");
	run("Select None");
	
	saveAs("Tif", out_path + "mask_epithelium.tif");

	selectWindow("ROI Manager"); 
	run("Close");
	run("Close All");
}

function getEpitheliumArea(sni, out_path) {

	open(out_path + "mask_epithelium.tif");
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	run("Set Measurements...", "area mean standard min redirect=None decimal=4");
	run("Analyze Particles...", "clear summarize");
	Table.rename("Summary", "Results");
	tot_area = getResult("Total Area", 0);
	close("mask_epithelium.tif");

	selectWindow("Results"); 
	run("Close");

	return(tot_area);
}

function euclideanDist(x0,y0,x1,y1) {
	return Math.sqrt(Math.pow((x1-x0), 2) + Math.pow(y1-y0, 2));
}

function containsSNI(querySNI) {
	for(i=0; i<neur_channels.length; i++) 
		if(matches(querySNI, ".*" + neur_channels[i] + ".*"))
			return true;
			
	return false;
}

function getCombinedNeuriteness(sni, out_path) {
	// iterate through all channels and search for the one(s) of interest (neur_channels)
	sniList = getFileList(path+sni);
	n = 0;

	// open mask_epithelium to crop neuriteness image
	open(out_path + "mask_epithelium.tif");
	rename("mask_epithelium");
	run("Create Selection");
	
	for(i=0; i<sniList.length; i++) {
		image = sniList[i];
		
		if(containsSNI(image) && startsWith(image, "Ne_")) { // if channel of interest
			
			// crop polygon mask on the neuriteness image, apply threshold and get threshold value
			run("Bio-Formats Importer", "open=[" + path + sni + "/" + image + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
			run("Median...", "radius=4");
			rename("neuriteness_"+n);
			run("Restore Selection");
			run("Clear Outside");
			run("Select None");

			run("Duplicate...", "mask_neuriteness_" + n);
			rename("mask_neuriteness_" + n);
			run("Restore Selection");
						
			setAutoThreshold(threshold_method+" dark");
			getThreshold(lower, upper);
			run("Convert to Mask");
			close("mask_neuriteness_" + n);
			
			// apply correction factor to "lower" threshold
			selectWindow("neuriteness_" + n);
			run("Restore Selection");
			setThreshold((lower*correc_factor), 255);
			setOption("BlackBackground", true);
			run("Convert to Mask");
			
			n++;
		}
	}

	selectWindow("neuriteness_0");
	rename("neuriteness");
	
	for(i=1; i<neur_channels.length; i++) {
		imageCalculator("Add create", "neuriteness","neuriteness_"+i);
		close("neuriteness");
		close("neuriteness_"+i);
		selectWindow("Result of neuriteness");
		rename("neuriteness");
	}
	saveAs("Tif", out_path + "thresholded_neuriteness_combined.tif");
	run("Close All");
}

function getNeuriteness(sni, image, channel, out_path) {
	
	// open mask_epithelium to crop neuriteness image
	open(out_path + "mask_epithelium.tif");
	rename("mask_epithelium");
	run("Create Selection");
			
	// crop polygon mask on the neuriteness image, apply threshold and get threshold value
	open(path + sni + "/" + image);
	run("Median...", "radius=4");
	rename("neuriteness");
	run("Restore Selection");
	run("Clear Outside");
	run("Select None");

	run("Duplicate...", "mask_neuriteness");
	rename("mask_neuriteness");
	run("Restore Selection");
				
	setAutoThreshold(threshold_method+" dark");
	getThreshold(lower, upper);
	run("Convert to Mask");
	close("mask_neuriteness");
	
	// apply correction factor to "lower" threshold
	selectWindow("neuriteness");
	run("Restore Selection");
	setThreshold((lower*correc_factor), 255);
	setOption("BlackBackground", true);
	run("Convert to Mask");

	saveAs("Tif", out_path + "thresholded_neuriteness_" + channel + ".tif");
	run("Close All");
}

function createLayer(ind, out_path) {
	roiManager("Select", ind);
	rName = Roi.getName();
	run("Create Mask");
	
	open(out_path + "mask_epithelium.tif");
	rename("epithelium");
	run("Create Selection");
	selectWindow("Mask");
	run("Restore Selection");
	run("Clear Outside");
	run("Select None");

	run("Set Measurements...", "area mean standard min redirect=None");
	run("Set Scale...", "distance=1 known=1 unit=unit");
	run("Analyze Particles...", "size=500-Infinity show=Masks");
	setOption("BlackBackground", true);
	run("Convert to Mask");

	if(rName == "A") rename("apical");
	else rename("basal");

	close("Mask");
	close("epithelium");
}

function getEpitheliumHeight(sni, out_path) {
	measures = ";;;;"; // initialize buffer of measures. If it cannot be calculated, then return string csv format
	
	// check if corresponding file exists in max projection folder
	if(File.exists(maxProjPath + sni + "_AB.zip")) {

		open(maxProjPath + sni + ".tif");
		// create image that will be the binary mask of the epithelium
		newImage("mask_epithelium", "8-bit black", getWidth(), getHeight(), 1);
		roiManager("Open", maxProjPath + sni + "_AB.zip");
		roiManager("Show All");
		close(maxProjPath + sni + ".tif");

		countRois = roiManager("count");
		
		if(countRois == 2) { // then create contours of apical and basal layers and measure average height
			createLayer(0, out_path);
			selectWindow("mask_epithelium");
			createLayer(1, out_path);
			selectWindow("ROI Manager"); 
			run("Close");
			
			AB_dist = calculateDistance("apical", "basal", true);
			BA_dist = calculateDistance("basal", "apical", true);
			run("Close All");

			measures = AB_dist + ";" + BA_dist;
		} else {
			print("The number of input layers is different than 2");
		}
	} else {
		print("File " + maxProjPath + sni + "_AB.zip does not exist");
	}

	return measures;
}

function calculateDistance(refLayer, targetLayer, areaValue) {
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

	selectWindow(refLayer);
	run("Invert");

	selectWindow("Results"); 
	run("Close");
	selectWindow("ROI Manager"); 
	run("Close");

	buffer = "";
	if(areaValue)
		buffer = "" + area + ";" + mean + ";" + std + ";" + min + ";" + max;
	else
		buffer = "" + mean + ";" + std + ";" + min + ";" + max;

	return buffer;
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

function segmentNuclei(sni, out_path) {

	// open mask_epithelium to crop nuclei image
	open(out_path + "mask_epithelium.tif");
	rename("mask_epithelium");
	run("Create Selection");

	// find nuclei image
	sniList = getFileList(path+sni);
	for(i=0; i<sniList.length; i++) {
		if(matches(sniList[i], ".*DAPI.*")) {
			run("Bio-Formats Importer", "open=[" + path + sni + "/" + sniList[i] + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			rename("nuclei");

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
}