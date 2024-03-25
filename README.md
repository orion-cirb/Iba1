# Iba1 

* **Developed for:** Julie
* **Team:** Rouach
* **Date:** March 2024
* **Software:** Fiji

### Images description

3D images taken with a x40 objective on  a spinning_disk microscope

* 1 channel: Iba1 microglia

With each image can be provided a *.roi* or *.zip* file containing one or multiple ROI(s).

### Plugin description

* Detect microglia soma with Cellpose
* Segment entire microglia with median filtering + thresholding
* Compute background noise of Iba1 channel
* Give soma number + microglia total volume + microglia background-corrected mean and integrated intensity
* If ROI(s) provided, remove from the analysis microglia that are inside

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *cyto2_Iba1_microglia* (homemade) model

### Version history

Version 1 released on March 25, 2024.
