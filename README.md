# Iba1 

* **Developed for:** Julie
* **Team:** Rouach
* **Date:** March 2024
* **Software:** Fiji

### Images description

3D images taken with a x40 objective on  a spinning-disk microscope

1 channel: Iba1 microglia

With each image can be provided a *.roi* or *.zip* file containing one or multiple ROI(s).

### Plugin description

* Detect microglial somas with Cellpose
* Segment microglial cells with median filtering + thresholding
* Compute background noise of Iba1 channel
* Give microglial somas number + microglial cells volume + microglial cells background-corrected mean and integrated intensity
* If ROI(s) provided, remove from the analysis microglia that are inside

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *cyto2_Iba1_microglia* (homemade) model

### Version history

Version 1 released on March 26, 2024.
