import Iba1_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import loci.common.DebugTools;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


/**
* Perform maximum intensity z-projection
* Detect microglial cells soma with Cellpose and give their number
* Segment microglial cells with thresholding ang give their area and integrated intensity
* @author Héloïse Monnet
*/
public class Iba1 implements PlugIn {

    private Iba1_Tools.Tools tools = new Tools();
       
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            }
            
            String imageDir = IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                return;
            }
            
            // Find images with fileExt extension
            String fileExt = tools.findImageType(imageDir);
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create OME-XML metadata store of the latest schema version
            DebugTools.setRootLevel("warn");
            ServiceFactory factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String channel = tools.dialog(imageDir, channels);
            if (channel == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            // Create output folder
            String outDirResults = imageDir + File.separator + "Results_" + tools.cellThMethod + "_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter fwResults = new FileWriter(outDirResults + "results.csv", false);
            BufferedWriter results = new BufferedWriter(fwResults);
                results.write("Image name\tImage vol (µm3)\tImage-ROI vol (µm3)\tIba1 bg\tSomas number\t"
                        + "Cells volume (µm3)\tCells bg-corr mean intensity\tCells bg-corr integrated intensity\n");
            results.flush();
            
            for (String f: imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // Check if ROIs file exists, keep rois to clear regions containing "artefacts"
                String roiName = imageDir + File.separator + rootName; 
                roiName = new File(roiName + ".zip").exists() ? roiName + ".zip" : roiName + ".roi";
                ArrayList<Roi> rois = new ArrayList<>();
                if (new File(roiName).exists()) {
                    RoiManager rm = new RoiManager(false);
                    rm.reset();
                    rm.runCommand("Open", roiName);
                    Collections.addAll(rois, rm.getRoisAsArray());
                }
                
                // Open Iba1 channel
                tools.print("- Opening Iba1 channel -");
                int index = ArrayUtils.indexOf(channels, channel);
                ImagePlus imgIba1 = BF.openImagePlus(options)[index];
                
                // Segment Iba1 cells
                tools.print("- Segmenting Iba1 cells -");
                Object3DInt cellObj = tools.segmentation(imgIba1, rois);
                
                // Detect Iba1 soma with Cellpose
                tools.print("- Detecting Iba1 somas -");
                Objects3DIntPopulation somaPop = tools.cellposeDetection(imgIba1, rois, cellObj);
                
                // Computing Iba1 background noise
                tools.print("- Computing Iba1 background noise -");
                double bg = tools.computeBackgroundNoise(imgIba1);
                
                // Write results
                tools.print("- Writing results -");
                double imgVol = imgIba1.getWidth() * imgIba1.getHeight() * imgIba1.getNSlices() * tools.pixVol;
                double roisVol = tools.getRoisVolume(rois, imgIba1);
                MeasureVolume mv = new MeasureVolume(cellObj);
                MeasureIntensity mi = new MeasureIntensity(cellObj, ImageHandler.wrap(imgIba1));
                results.write(rootName+"\t"+imgVol+"\t"+(imgVol-roisVol)+"\t"+bg+"\t"+somaPop.getNbObjects()+"\t"+
                        mv.getVolumeUnit()+"\t"+(mi.getValueMeasurement(MeasureIntensity.INTENSITY_AVG) - bg)+"\t"+
                        (mi.getValueMeasurement(MeasureIntensity.INTENSITY_SUM) - bg*mv.getVolumePix())+"\n");
                results.flush();
                        
                // Draw results
                tools.print("- Drawing results -");
                tools.drawResults(somaPop, cellObj, imgIba1, outDirResults+rootName+".tif");
                
                tools.closeImage(imgIba1);
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Iba1.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("All done!");
    }
}
