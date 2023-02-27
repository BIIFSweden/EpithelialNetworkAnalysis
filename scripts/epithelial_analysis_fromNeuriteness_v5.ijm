
// ===============================
// AUTHOR : Gisele Miranda 
// IVESTIGATORS : Annelie Tjernlund, Mathias Franzén Boger, Gabriela Edfeldt, Kristina Broliden
// CREATE DATE : 2020 - 06 - 18
// PURPOSE : Structural analysis of the cervical epithelial tissue
// NOTES : Required plugins - MorphoLibJ (https://imagej.net/MorphoLibJ)
// 							- Multi Otsu Threshold (https://imagej.net/plugins/multi-otsu-threshold)
// ===============================


/********************** parameters ***********************/

// path to the input directory containing the SNIs
path = "/Users/giselemiranda/Downloads/Gisele Cy5 (ZO1)/HIV/";
// path to the input directory containing the maximum intensity projection that were used to draw AB lines
maxProjPath = "/Users/giselemiranda/Downloads/Gisele Cy5 (ZO1)/HIV_Max/";
// channel to be analyzed - names should be used according to the nomenclature of the files: cy3, FITC, Cy5 and m cherry
channel_of_interest = "Cy5";
// method chosen to threshold the NEURITENESS image
threshold_method = "Otsu";
// correction factor applied to the segmented NEURITENESS
correc_factor = 0.8; //1.5;
// distance threshold used to segment the REGION OF INTEREST, defined based on the EDT applied over NEURITENESS IMAGE
distance_threshold = 7; //10;
// area threshold to filter the segmented regions
area_threshold = 70000;
// area threshold to filter the connected components of the intact net
area_connected_threshold = 5000;
// size of the structuring element for Dilation
radiusDilation = 5;
// minimum nuclei size used to post process the nuclei segmentation
nucleiSize = 1000;
// scale used to calibrate the images (pixels per micron)
pixPerMic = 3.07693;
// if scale should be used, otherwise results will be given in pixels
useScale = true;
// used to set up minimum and maximum threshold values to segment the neuriteness image
minThr = 0; //0;
maxThr = 1;
// name of the experiment that will be used to name the output folder
exp_name = "Cy5_pipeline5";

/*********************************************************/

dir = getFileList(path);
bufferMeasures = "";
exp_id = threshold_method + "_correc-factor-" + correc_factor + "_dist-" + distance_threshold + "_area-" + area_threshold + "_radDilation-" + radiusDilation + "_" + exp_name;
overlayPath = path + "overlay_" + exp_id + "/";
if(!File.exists(overlayPath)) File.makeDirectory(overlayPath);

for(countSNI=0; countSNI<dir.length; countSNI++) { // for each SNI
	print("countSNI: " + countSNI);
	sni = dir[countSNI];
	images = getFileList(path+sni);
	
	// check if input folder is a SNI. It should start with "SNI" or "Neg SNI" (for negative control)
	if(startsWith(sni, "SNI") || startsWith(sni, "Neg SNI")) {
		outpath = path + sni + "output_" + exp_id + "/";
		if(!File.exists(outpath)) File.makeDirectory(outpath);
		
		// crop epithelium region
		sni = replace(sni, "/", "");
		print(sni);

		// check if original image, MIP image, AB lines and neuriteness exist
		img_path = path + sni + "/newTiffImages/" + sni + "_" + channel_of_interest + "_Extended.tif"; // original image (fluorescent channel)
		neur_path = path + sni + "/" + "Ne_" + sni + "_" + channel_of_interest + "_Extended.tif"; // processed neuriteness
		img_path_max = maxProjPath + sni + ".tif"; // MIP of original channels
		img_path_max_AB = maxProjPath + sni + "_AB.zip"; // AB lines
		dapi_path = path + sni + "/newTiffImages/" + sni + "_DAPI_Extended.tif"; // DAPI path

		bufferMeasures = bufferMeasures + sni + ";" + channel_of_interest + ";";
		
		if(File.exists(img_path) && File.exists(neur_path) && File.exists(img_path_max) && File.exists(img_path_max_AB) && File.exists(dapi_path)) { //check if FITC channel and the corresponding max projection exists for the segmentation of the epithelium

			// get epithelium
			getEpithelium(sni, img_path, outpath);
			
			// get average epithelium height
			epithHeight = getEpitheliumHeight(sni, img_path_max, outpath);
			
			// get total area of the epithelium
			tot_area_epithelium = getEpitheliumArea(sni, outpath);
	
			// segment nuclei, if dapi channel exists
			segmentNuclei(dapi_path, outpath);
			
			// get MFI of whole epithelium
			epithMFI = getEpitheliumMFI(img_path, outpath);
			
			// open corresponding neuriteness image and start processing it
			getNeuriteness(channel_of_interest, neur_path, outpath); // segmentation
			open(outpath + "thresholded_neuriteness_" + channel_of_interest + ".tif");
			rename("thresholded_neuriteness");
			run("Invert");
			run("Exact Euclidean Distance Transform (3D)");
			rename("EDT");
			run("Conversions...", "weighted");
			setOption("ScaleConversions", false);
			run("16-bit");
			setThreshold(0, distance_threshold);
			run("Convert to Mask");
			
			// select only regions with area greater than area_threshold
			selectWindow("EDT");
			run("Set Measurements...", "area mean standard min redirect=None decimal=4");
			run("Analyze Particles...", "size="+area_threshold+"-Infinity show=Masks clear add"); // filtering
			rename("Mask");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			countSelecArea = roiManager("count");
			selectWindow("ROI Manager"); 
			run("Close");

			// fill buffer of measures for the epithelium
			bufferMeasures = bufferMeasures + epithHeight + ";" + tot_area_epithelium + ";" + epithMFI + ";";

			// if there is a segmented ROI, then create fill buffer to store measures of the segmented ROI (bufferChannel)
			bufferChannel = "";
			
			if(countSelecArea > 0) { // then the region of intesest is not empty
				// check intersections with the apical layer
				selectWindow("thresholded_neuriteness");
				run("Invert");
				selectWindow("Mask");
				run("Select None");
				run("Morphological Filters", "operation=Dilation element=Square radius="+radiusDilation);
				rename("segmented_region");
				run("Duplicate...", "segmented_region_final");
				rename("segmented_region_final");
				selectWindow("segmented_region");
				close("Mask");
				
				//generate a mask for the apical layer to verify intersection with the segmented regions
				newImage("Temp", "8-bit black", getWidth(), getHeight(), 1);
				roiManager("Open", maxProjPath + sni + "_AB.zip");
				roiManager("Show All");
				roiManager("select", 0);
				if(!matches(Roi.getName(), "A")) roiManager("select", 1);
				run("Create Mask");
				rename("mask_apical");
				selectWindow("ROI Manager"); 
				run("Close");
				close("Temp");

				// open mask_epithelium and crop neuriteness image
				open(outpath + "mask_epithelium.tif");
				run("Create Selection");
				selectWindow("segmented_region");
				run("Restore Selection");
				run("Clear Outside");
				run("Select None");
				
				// iterate over particles of the segmented region to check intersection with apical layer
				run("Analyze Particles...", "size=0-Infinity clear add");
				nRois = roiManager("count");
				indexToDelete = newArray(nRois);
				for(i=(nRois-1); i>=0; i--) {
					roiManager("Select", i);
					run("Create Mask");
					roiManager("Measure");
					area_i = getResult("Area", 0);
					
					imageCalculator("AND create", "Mask","mask_apical");
					selectWindow("Result of Mask");
					run("Analyze Particles...", "clear summarize");
					Table.rename("Summary", "Results");
					area_intersection = getResult("Total Area", 0);
					
					if(matches(channel_of_interest, "cy3") || matches(channel_of_interest, "FITC")) {
						if(area_intersection != 0 && area_i < area_threshold) indexToDelete[i] = 1;
						else indexToDelete[i] = 0;
					} else {
						if(area_i < area_threshold) indexToDelete[i] = 1;
						else indexToDelete[i] = 0;	
					}
					
					close("Mask");
					close("Result of Mask");
				}
				
				// delete indexToDelete, the remaining selections that are larger the threshold
				countSelectedROIs = 0;
				for(i=0; i<nRois; i++) {
					if(indexToDelete[i] == 0) {
						countSelectedROIs++;
						roiManager("Select", i);
						run("Create Mask");
					}
				}
				
				print("countSelectedROIs: " + countSelectedROIs);
				run("Select None");
				selectWindow("ROI Manager"); 
				run("Close");
				selectWindow("Results"); 
				run("Close");
				close("mask_epithelium.tif");
				
				// if there are segmented regions, then calculate derived measures
				if(countSelectedROIs != 0) { 
					
					// obs: include this code if holes in the thresholded neuriteness mask should not be included in the final segmented region
					// create selection over the final mask and crop the corresponding boundary on the original neuriteness
					/*selectWindow("Mask");
					//run("Create Selection");
					//selectWindow("segmented_region_final");
					//run("Restore Selection");
					//run("Clear Outside");
					//run("Select None");
					//selectWindow("segmented_region_final");*/
					
					selectWindow("Mask");

					// total segmented area
					if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
					run("Analyze Particles...", "clear summarize");
					//run("Analyze Particles...", "clear include summarize");
					Table.rename("Summary", "Results");
					tot_area_segmented = getResult("Total Area", 0);
					selectWindow("Results"); 
					run("Close");
					
					saveAs("Tif", outpath + "segmented_region_" + channel_of_interest + ".tif");
					rename("intact_net");
					run("Create Selection");
	
					selectWindow("thresholded_neuriteness");
					run("Restore Selection");
					run("Flatten");
					saveAs("Tif", outpath + "segmented_region_over_neuriteness_" + channel_of_interest + ".tif");
					close();

					selectWindow("thresholded_neuriteness");
					close("\\Others");
					run("Clear Outside");
					run("Select None");
					saveAs("Tif", outpath + "thresholded_neuriteness_segmented_region_" + channel_of_interest + ".tif");
					
					if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
					else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
					run("Analyze Particles...", "clear summarize");
					Table.rename("Summary", "Results");
					tot_area_segmented_neuriteness = getResult("Total Area", 0);
					selectWindow("Results");
					run("Close");

					// get MFI from original channel image
					open(img_path);
					rename("original");
					run("Restore Selection");
					
					run("Measure");
					orig_mean = getResult("Mean", 0);
					orig_std = getResult("StdDev", 0);
					selectWindow("Results"); 
					run("Close");

					selectWindow("original");
					run("Flatten");
					saveAs("Tif", outpath + "segmented_region_overlay_original_" + channel_of_interest + ".tif");
					saveAs("Tif", overlayPath + "segmented_region_overlay_original_" + channel_of_interest + "_" + sni);
					
					// save segmented region over neuriteness image
					open(neur_path);
					rename("neuriteness");
					run("Restore Selection");

					selectWindow("neuriteness");
					run("Flatten");
					saveAs("Tif", outpath + "segmented_region_overlay_neuriteness_" + channel_of_interest + ".tif");
					saveAs("Tif", overlayPath + "segmented_region_overlay_neuriteness_" + channel_of_interest + "_" + sni);

					run("Close All");

					// get parabasal, upper and lower layers
					open(outpath + "mask_epithelium.tif");
					rename("mask_epithelium");
					open(outpath + "segmented_region_" + channel_of_interest + ".tif");
					rename("mask_channel");
					open(outpath + "mask_nuclei.tif");
					rename("mask_nuclei");
					
					run("Duplicate...", "temp_ROIs");
					rename("temp_ROIs");

					roiManager("Open", img_path_max_AB);
					roiManager("Show All");
					createLayer(0, outpath);
					selectWindow("temp_ROIs");
					createLayer(1, outpath);
					selectWindow("ROI Manager"); 
					run("Close");
					close("temp_ROIs");

					// get area of parabasal layer by adding the mask_channel to the mask_nuclei
					imageCalculator("Add create", "mask_channel","mask_nuclei");
					selectWindow("Result of mask_channel");
					rename("parabasal_marker");
					run("Fill Holes");

					selectWindow("parabasal_marker");
					if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
					else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
					run("Set Measurements...", "area mean standard min redirect=None decimal=4");
					run("Analyze Particles...", "clear include summarize");
					Table.rename("Summary", "Results");
					tot_area_parabasal = getResult("Total Area", 0);
					selectWindow("Results"); 
					run("Close");

					// get area of upper layer (between parabasal and superficial layer)
					imageCalculator("Subtract create", "mask_epithelium","parabasal_marker");
					selectWindow("Result of mask_epithelium");
					run("Analyze Particles...", "pixel show=Nothing display clear include add");
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
					run("Analyze Particles...", "clear include summarize");
					Table.rename("Summary", "Results");
					tot_area_lower = getResult("Total Area", 0);
					selectWindow("Results"); 
					run("Close");

					selectWindow("parabasal_marker");
					saveAs("Tif", outpath + "parabasal_layer_" + channel_of_interest + ".tif");
					selectWindow("upper");
					saveAs("Tif", outpath + "upper_layer_" + channel_of_interest + ".tif");
					selectWindow("lower");
					saveAs("Tif", outpath + "lower_layer_" + channel_of_interest + ".tif");

					// add to buffer
					bufferChannel = bufferChannel + tot_area_segmented + ";" + orig_mean + ";" + orig_std + ";" + tot_area_segmented_neuriteness + ";" + (tot_area_segmented/tot_area_epithelium) + ";" + (tot_area_segmented_neuriteness/tot_area_epithelium) + ";" + (tot_area_segmented_neuriteness/tot_area_segmented) + ";" + tot_area_parabasal + ";" + tot_area_upper + ";" + tot_area_lower + ";";
					
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

					tot_area_neuriteness = 0;
					intact_net_measures = 0;
					flooded_area = 0;
					
					tot_area_neuriteness = getNeuritenessArea(channel_of_interest, outpath);
					intact_net_measures = getIntactNet_ConnectedComponent(img_path, neur_path, channel_of_interest, sni, outpath, overlayPath);
					//this functions also returns the area, however we're not using here now
					getIntactNet(channel_of_interest, sni, outpath, overlayPath);
					// getFloodedArea depends on the results generated on function getIntactNet_ConnectedComponent
					flooded_area = getFloodedArea(channel_of_interest, outpath, overlayPath, img_path_max_AB, sni);
					if (flooded_area == -1) flooded_area = tot_area_epithelium; // then flooded_area will be equal the area of the epithelium
					
					intact_net_measures = split(intact_net_measures," ");
					area_intact_net_cc_holes = parseFloat(intact_net_measures[0]);
					area_intact_net_cc = parseFloat(intact_net_measures[1]);
					MFI_intact_net_cc = parseFloat(intact_net_measures[2]);
					std_MFI_intact_net_cc = parseFloat(intact_net_measures[3]);
					
					// update buffers
					bufferMeasures = bufferMeasures + bufferChannel + ";" + area_intact_net_cc + ";" + area_intact_net_cc_holes + ";" + tot_area_neuriteness + ";" + (area_intact_net_cc/tot_area_neuriteness) + ";" + MFI_intact_net_cc + ";" + std_MFI_intact_net_cc + ";" + flooded_area + ";" + (flooded_area/tot_area_epithelium) + "\n";
				} else {
					print("No regions were selected for this SNI");
					run("Close All");
					bufferMeasures = bufferMeasures + "\n"; // complete csv buffer with empty symbols
				}
			} else {
				print("No regions with area bigger than the area threshold were found");
				run("Close All");
				bufferMeasures = bufferMeasures + "\n"; // complete csv buffer with empty symbols
			}
		} else {
			run("Close All");
			bufferMeasures = bufferMeasures + "\n";
			if(!File.exists(img_path)) print("File: " + path + sni + "/newTiffImages/" + sni + "_" + channel_of_interest + "_Extended.tif does not exist");
			if(!File.exists(img_path_max)) print("File: " + maxProjPath + sni + ".tif does not exist");
			if(!File.exists(img_path_max_AB)) print("File: " + maxProjPath + sni + "_AB.zip does not exist");
			if(!File.exists(neur_path)) print("File: " + path + sni + "/" + "Ne_" + sni + "_" + channel_of_interest + "_Extended.tif does not exist");
			if(!File.exists(dapi_path)) print("File: " + path + sni + "/newTiffImages/" + sni + "_DAPI_Extended.tif does not exist");
		}
	}
}

summaryPath = path + "summary_" + channel_of_interest + "_" + exp_id + ".csv";

if(File.exists(summaryPath))
	File.delete(summaryPath);

summaryFile = File.open(summaryPath);
print(summaryFile, "sni;channel;number of pixels in the basal layer;average height basal layer;std basal layer;min basal layer;max basal layer;number of pixels in the apical layer;height apical layer;std apical layer;min apical layer;max apical layer;area - epithelium; MFI - epithelium; area - marker of interest;MFI - marker of interest;stdMFI - marker of interest;area - neuriteness;relative area - marker of interest per total epithelium;realtive area - neuriteness per total epithelium;realtive area - neuriteness per marker of interest;area - marker of interest + parabasal layer;area - upper layer;area - lower layer;average height between marker of interest and apical layer;std height between marker of interest and apical layer;min height between marker of interest and apical layer;max height between marker of interest and apical layer;average height between marker of interest and basal layer;std height between marker of interest and basal layer;min height between marker of interest and basal layer;max height between marker of interest and basal layer;area - intact net;area - intact net with holes;neuriteness total area;relative area - intact net per neuriteness total area;MFI intact net (original);stdMFI intact net (original);flooded area;relative area - flooded area per total epithelium\n");
print(summaryFile, bufferMeasures);
File.close(summaryFile);

/****************** functions ****************************/

function getEpithelium(sni, input_path, out_path) {

	// read original FITC channel and retrieve the binary mask of the tissue, excluding the black background
	run("Bio-Formats Importer", "open=[" + input_path + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename("mask_tissue");

	// rough threshold to separate the background
	if(bitDepth() == 8) setThreshold(1, 255);
	else setThreshold(5, 65535);
	
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Fill Holes");
	run("Create Selection");
	run("Enlarge...", "enlarge=-3"); // cut some pixels of the border to remove the bright square around the SNI 
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

function getEpitheliumHeight(sni, input_path, out_path) {
	measures = "";

	open(input_path);
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

	return measures;
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

function getEpitheliumMFI(input_path, out_path) {
	mfi = 0;
	
	// check if corresponding file exists in max projection folder
	if(File.exists(out_path + "mask_epithelium.tif")) {
		
		open(out_path + "mask_epithelium.tif");
		run("Create Selection");
		open(input_path);
		rename("original");
		
		run("Restore Selection");
		run("Set Measurements...", "area mean redirect=None decimal=4");
		run("Measure");
		mfi = getResult("Mean", 0);
	} 
	
	selectWindow("Results"); 
	run("Close");
	
	run("Close All");

	return mfi;
}

function segmentNuclei(dapi_path, out_path) {
	// open mask_epithelium to crop nuclei image
	open(out_path + "mask_epithelium.tif");
	rename("mask_epithelium");
	run("Create Selection");
	
	open(dapi_path);
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

function getNeuriteness(channel, input_path, out_path) {
	
	// open mask_epithelium to crop neuriteness image
	open(out_path + "mask_epithelium.tif");
	rename("mask_epithelium");
	run("Create Selection");
			
	// crop polygon mask on the neuriteness image, apply threshold and get threshold value
	open(input_path);
	run("Median...", "radius=4");
	rename("neuriteness");
	run("Restore Selection");
	run("Clear Outside");
	run("Select None");

	run("Duplicate...", "mask_neuriteness");
	rename("mask_neuriteness");
	run("Multi OtsuThreshold", "numlevels=3");
	
	/*selectWindow("Region 1");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");

	selectWindow("Region 2");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	
	imageCalculator("Add create", "Region 1","Region 2");*/

	logString = getInfo("log");
	buffer = split(logString, "\n");
	buffer = split(buffer[buffer.length-1], " ");
	
	lower = "";
	if(matches(channel, "cy3") || matches(channel, "FITC"))
		lower = split(buffer[2], "=");
	else
		lower = split(buffer[3], "=");
	lower = lower[1];
	lower = replace(lower, ",", "");
	lower = parseInt(lower);
	
	lower = maxOf((lower*correc_factor), minThr*255);
	upper = minOf(255, maxThr*255);
	
	// apply correction factor to "lower" threshold
	selectWindow("neuriteness");
	setThreshold(lower, upper);
	setOption("BlackBackground", true);
	run("Convert to Mask");

	saveAs("Tif", out_path + "thresholded_neuriteness_" + channel + ".tif");
	run("Close All");
}

// returns the area of the entire neuriteness network
function getNeuritenessArea(channel, out_path) {
	// open neuriteness image
	open(out_path + "thresholded_neuriteness_" + channel + ".tif");
	
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	run("Set Measurements...", "area mean standard min redirect=None decimal=4");
	run("Analyze Particles...", "clear summarize");
	Table.rename("Summary", "Results");
	tot_area = getResult("Total Area", 0);

	selectWindow("Results"); 
	run("Close");
	
	run("Close All");
	
	return(tot_area);
}

function getIntactNet_ConnectedComponent(img_path, neur_path, channel, sni, out_path, overlay_path) {
	// open original channel image
	//run("Bio-Formats Importer", "open=[" + path + sni + "/newTiffImages/" + img_file + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
	open(img_path);
	rename("original");
	if(!useScale) run("Set Scale...", "distance=1 known=1 unit=micron");
	
	open(out_path + "thresholded_neuriteness_segmented_region_" + channel + ".tif");
		
	rename("inputImg");
	run("Duplicate...", "markerImg");
	rename("markerImg");
	
	// filter connected components by size and create filtered mask
	run("Analyze Particles...", "size="+area_connected_threshold+"-Infinity show=Masks clear include add"); // include holes
	count = roiManager("count");
	selectWindow("ROI Manager"); 
	run("Close");
	
	if(count != 0) {
		
		run("Convert to Mask");
		run("Create Selection");
		selectWindow("markerImg");
		run("Restore Selection");
		run("Clear Outside");
		run("Select None");
	
		selectWindow("Mask of markerImg");
		if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		run("Measure");
		area_with_holes = getResult("Area", 0);
		selectWindow("Results"); 
		run("Close");

		// measure MFI of filtered connected components and save LUT image (glasbey) - not taking into account holes in the binary mask
		selectWindow("markerImg");
		run("Duplicate...", "markerImg_intactNet");
		rename("markerImg_intactNet");
		
		run("Set Measurements...", "area mean standard min redirect=original decimal=4");
		run("Analyze Particles...", "size="+area_connected_threshold+"-Infinity show=[Count Masks] clear");
		run("glasbey");
		saveAs("PNG", out_path + "intact_net_connected_" + channel + ".png");
		saveAs("PNG", overlay_path + "intact_net_connected_" + channel + "_" + sni);
		
		selectWindow("markerImg");
		run("Create Selection");
		selectWindow("original");
		run("Restore Selection");
		run("Measure");
		
		if(nResults > 0) {
			tot_area = getResult("Area", 0);
			mfi = getResult("Mean", 0);
			stdMfi = getResult("StdDev", 0);
			selectWindow("Results"); 
			run("Close");
	
			// concatenate all measures 
			m = d2s(area_with_holes,4) + " " + d2s(tot_area,4) + " " + d2s(mfi,4) + " " + d2s(stdMfi,4);
		} 
	} else m = "0 0 0 0";
	
	run("Close All");
	
	print(m);
	return m;
}

function getIntactNet(channel, sni, out_path, overlay_path) {

	open(out_path + "thresholded_neuriteness_segmented_region_" + channel + ".tif");
		
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
	saveAs("Tif", overlay_path + "watershed_" + channel + "_" + sni);
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
	
	saveAs("Tif", out_path + "intact_net_watershed_" + channel + ".tif");
	saveAs("Tif", overlay_path + "watershed_" + channel + "_" + sni);

	run("Close All");
	return(tot_area);
}

function getFloodedArea(channel, out_path, overlay_path, maxProj_path, sni) {
	open(out_path + "mask_epithelium.tif");
	
	wid = getWidth();
	hei = getHeight();
	newImage("temp", "8-bit black", wid, hei, 1);
	
	roiManager("Open", maxProj_path);
	roiManager("Show All");
	
	roiManager("Select", 0);
	run("Create Mask");
	rename("apical");
	
	selectWindow("temp");
	roiManager("Select", 1);
	run("Create Mask");
	rename("basal");
	
	selectWindow("temp");
	close();
	selectWindow("ROI Manager"); 
	run("Close");

	if(File.exists(out_path + "intact_net_connected_" + channel + ".png")) {
		open(out_path + "intact_net_connected_" + channel + ".png");
		setThreshold(1, 65535, "raw");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		rename("neuriteness");
			
		selectWindow("mask_epithelium.tif");
		run("Duplicate...", "epithelium_outline");
		rename("epithelium_outline");
		run("Outline");
		
		// create marker for whatershed
		selectWindow("apical");
		run("Morphological Filters", "operation=Dilation element=Square radius=20");
		rename("marker");
		imageCalculator("AND create", "mask_epithelium.tif","marker");
		rename("marker_processed");
		run("Erode");
		
		imageCalculator("Add create", "epithelium_outline","neuriteness");
		rename("mask_watershed");
		run("Invert");
		run("Morphological Reconstruction", "marker=marker_processed mask=mask_watershed type=[By Dilation] connectivity=4");
		
		saveAs("Tif", overlay_path + "flooded_region_" + channel + "_" + sni);
		saveAs("Tif", out_path + "flooded_region_" + channel + ".tif");
		rename("marker_processed-rec");
		
		if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		run("Set Measurements...", "area mean standard min redirect=None decimal=4");
		run("Analyze Particles...", "clear summarize");
		Table.rename("Summary", "Results");
		area_flooded = getResult("Total Area", 0);
		selectWindow("Results"); 
		run("Close");
	} else area_flooded = -1;
	
	run("Close All");
	return(area_flooded);
}

function euclideanDist(x0,y0,x1,y1) {
	return Math.sqrt(Math.pow((x1-x0), 2) + Math.pow(y1-y0, 2));
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

