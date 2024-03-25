package Iba1_Tools;

import Iba1_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Iba1_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.Analyzer;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.atomic.DoubleAccumulator;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;


/**
 * @author Héloïse Monnet
 */
public class Tools {
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/Iba1";
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    public Calibration cal = new Calibration();
    public double pixVol;
    
    // Soma detection
    public String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    public final String cellposeModelPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\" : "";
    public String cellposeModel = "cyto2_Iba1_microglia";
    public int cellposeDiam = 60;
    public double cellposeStitchTh = 0.5;
    public double minSomaVol = 50;
    public double maxSomaVol = 800;
    
    // Cells segmentation
    public String cellThMethod = "Otsu";
    public double minCellVol = 0.5;
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ2 not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(String imagesFolder) {
        String ext = "";
        String[] files = new File(imagesFolder).list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "nd2" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }

        
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);     
    }
    
    
    /**
     * Generate dialog box
     */
    public String dialog(String imagesDir, String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channel", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Iba1: ", channels, channels[2]);
        
        gd.addMessage("Somas detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min volume (µm3): ", minSomaVol, 2);
        gd.addNumericField("Max volume (µm3): ", maxSomaVol, 2);
        
        gd.addMessage("Cells segmentation", Font.getFont("Monospace"), Color.blue);
        String[] thMethods = AutoThresholder.getMethods();
        gd.addChoice("Threshold method: ", thMethods, cellThMethod);
        gd.addNumericField("Min volume (µm3): ", minCellVol, 2);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm): ", cal.pixelHeight, 3);
        gd.addNumericField("Z calibration (µm): ", cal.pixelDepth, 3);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String channel = gd.getNextChoice();
        
        minSomaVol = gd.getNextNumber();
        maxSomaVol = gd.getNextNumber();
        
        cellThMethod = gd.getNextChoice();
        minCellVol = gd.getNextNumber();
        
        cal.pixelHeight = cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixVol = cal.pixelHeight*cal.pixelWidth*cal.pixelDepth;
        
        if (gd.wasCanceled())
            channel = null;
        return(channel);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
       
    /**
     * Detect objects in 3D using 2D-stitched version of Cellpose
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus imgIn, ArrayList<Roi> rois, Object3DInt obj) {
        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModelPath+cellposeModel, 1, cellposeDiam, cellposeEnvDir);
        settings.setStitchThreshold(cellposeStitchTh);
        settings.useGpu(true);

        // Run Cellpose
        ImagePlus img = imgIn.duplicate();
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
        ImagePlus imgOut = cellpose.run();        
        imgOut.setCalibration(cal);
        
        // Fill ROIs in black
        if (!rois.isEmpty())
            fillImg(imgOut, rois);
        
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        System.out.println("Nb objects detected:"+pop.getNbObjects());
        popFilterZ(pop);
        popFilterSize(pop, minSomaVol, maxSomaVol);
        popFilterColoc(pop, obj);
        System.out.println("Nb objects remaining after filtering: "+ pop.getNbObjects());
        
        closeImage(img);
        closeImage(imgOut);
        return(pop);
    }
    
    
    /**
     * Fill ROIs in black in image
     */
    public ImagePlus fillImg(ImagePlus img, ArrayList<Roi> rois) {
        img.getProcessor().setColor(Color.BLACK);
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setSlice(s);
            for (Roi r : rois) {
                img.setRoi(r);
                img.getProcessor().fill(img.getRoi());
            }
        }
        img.deleteRoi();
        return(img);
    } 
    
    
    /**
     * Remove objects that appear in only one z-slice
     */
    public void popFilterZ(Objects3DIntPopulation pop) {
        pop.getObjects3DInt().removeIf(p -> (p.getObject3DPlanes().size() == 1));
        pop.resetLabels();
    }
    
    
    /**
     * Remove objects in population with size < min and size > max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.setVoxelSizeXY(cal.pixelWidth);
        pop.setVoxelSizeZ(cal.pixelDepth);
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    public void popFilterColoc(Objects3DIntPopulation somaPop, Object3DInt cellMask) {
        Objects3DIntPopulation cellPop = new Objects3DIntPopulation();
        cellPop.addObject(cellMask);
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(somaPop, cellPop);
        somaPop.getObjects3DInt().removeIf(soma -> coloc.getValueObjectsPair(soma, cellMask) < 0.25*soma.size());
        somaPop.resetLabels();
    }
    
    
    /**
     * Segment objects in 2D with median filtering + thresholding
     */
    public Object3DInt segmentation(ImagePlus imgIn, ArrayList<Roi> rois) {
        ImagePlus imgMed = median3DSliceBySlice(imgIn, 2);
        ImagePlus imgOut = threshold(imgMed, cellThMethod);
        imgOut.setCalibration(cal);
        
        // Fill ROIs in black
        if (!rois.isEmpty())
            fillImg(imgOut, rois);
        
        Objects3DIntPopulation pop = getPopFromImage(imgOut);
        System.out.println("Nb objects detected:"+pop.getNbObjects());
        popFilterSize(pop, minCellVol, Double.MAX_VALUE);
        System.out.println("Nb objects remaining after size filtering: "+ pop.getNbObjects());
        
        Object3DInt obj = convertPopToObj(pop, imgIn);
        
        closeImage(imgMed);
        closeImage(imgOut);
        return(obj);
    }
      
    
    /**
     * 2D median filtering slice by slice using CLIJ2
     */ 
    public ImagePlus median3DSliceBySlice(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DSliceBySliceSphere(imgCL, imgCLMed, sizeXY, sizeXY);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCL);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
    
    /**
     * Automatic thresholding using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    
    /**
     * Return population of 3D objects population from binary image
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageInt labels = new ImageLabeller().getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        labels.closeImagePlus();
        return(pop);
    }
    
    
    /**
     * Convert Objects3DIntPopulation to Object3DInt
     */
    public Object3DInt convertPopToObj(Objects3DIntPopulation pop, ImagePlus img) {
        ImageHandler imh = ImageHandler.wrap(img).createSameDimensions();
        for (Object3DInt obj: pop.getObjects3DInt())
            obj.drawObject(imh, 255);
        return(new Object3DInt(imh));
    }
    
        
    /**
     * Compute image background noise:
     * z-project over min intensity + read median intensity
     */
    public double computeBackgroundNoise(ImagePlus img) {
      ImagePlus imgProj = zProject(img, ZProjector.MIN_METHOD);
      double bg = imgProj.getProcessor().getStatistics().median;
      System.out.println("Background noise (median of the min projection) = " + bg);
      closeImage(imgProj);
      return(bg);
    }
    
    
    /**
     * Z-projection a stack
     */
    public ImagePlus zProject(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    
    /**
     * Compute ROIs total volume
     */
    public double getRoisVolume(ArrayList<Roi> rois, ImagePlus img) {
        double roisVol = 0;
        for(Roi roi: rois) {
            PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
            poly.setLocation(0, 0);
            
            img.resetRoi();
            img.setRoi(poly);

            ResultsTable rt = new ResultsTable();
            Analyzer analyzer = new Analyzer(img, Analyzer.AREA, rt);
            analyzer.measure();
            roisVol += rt.getValue("Area", 0);
        }

        return(roisVol * img.getNSlices() * cal.pixelDepth);
    }
    
    
    /**
     * Find total volume of objects in population
     */
    public double findPopVolume(Objects3DIntPopulation pop) {
        DoubleAccumulator sumVol = new DoubleAccumulator(Double::sum,0.d);
        pop.getObjects3DInt().parallelStream().forEach(obj -> { 
            sumVol.accumulate(new MeasureVolume(obj).getVolumeUnit());
        });
        return(sumVol.doubleValue());
    }
    
    
    /**
     * Draw results
     */
    public void drawResults(Objects3DIntPopulation somaPop, Object3DInt cellObj, ImagePlus img,  String name) {
        ImageHandler imhSoma = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imhCell = imhSoma.createSameDimensions();
        
        // Draw soma in red and cells in blue
        for(Object3DInt soma: somaPop.getObjects3DInt())
            soma.drawObject(imhSoma); //, 255 TODO
        cellObj.drawObject(imhCell, 255);

        ImagePlus[] imgColors = {imhSoma.getImagePlus(), null, imhCell.getImagePlus(), img};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
        
        imhSoma.closeImagePlus();
        imhCell.closeImagePlus();
    }
    
}
