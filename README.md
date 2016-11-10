## HTC_examples

### Part 0: Setup

Assuming that you have already installed [miniconda3](http://conda.pydata.org/miniconda.html) and that it resides in a subdirectory named "miniconda3" in your HTC home directory ...

`$ cd ~/.`

`$ export PATH=/home/<username>/miniconda3/bin:$PATH`

`$ python --version`

* make sure that you are using Python v3.x

`$ conda list`

* look for numpy, pandas, h5py, hdf5, and matplotlib
* install any of these that you don't already have

`$ conda install pandas h5py [etc.]`

Once you have all the right Python packages installed:

```
$ tar -czf python3.tar.gz miniconda3
$ git clone https://github.com/megarcia/HTCondor_examples.git
$ cp python3.tar.gz HTCondor_examples/
$ cd HTCondor_examples
$ chmod 755 */*.sh
```

### Part 1: CSV_processing example

This example illustrates the basic HTCondor job unit.

`$ cd csv_processing`

The data file *GHCND_WIS_2015.csv* is named to indicate source, geographic coverage, and temporal coverage. These are [Global Historic Climate Network-Daily (GHCND)](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn) surface weather observations for all 1st-order and cooperative stations in the rectangle containing Wisconsin for 2015. The dataset has a header line with column names, and includes a number of fields such as location, measurements, and data quality flags. The Python script *clean_GHCND.py* was written to "clean" this raw dataset, removing flagged measurements and stations without location information, etc. This is the same routine as the Python script *process_NCEI_00.py* that is used in my [WxCD project](https://megarcia.github.io/WxCD).

`$ condor_submit_dag clean_GHCND_dag.sub`

You can check the job status with

`$ condor_q $USER`

When the job is completed, you'll find a lot of processing (cleaning) information in the *clean_GHCND.out* file. The original input datafile is left intact, and you'll find new *.csv* files with various info on the errors encountered and the cleaned dataset. The results can be further used for temporal examination, spatial analysis, etc.

### Part 2: Image_processing example

This example illustrates running HTCondor workflows.

Before you start, you need to get the two Landsat image files from my [Google Drive](https://drive.google.com/open?id=0B4-FFhSfVlLyQnNrbVlDeUhyZG8). They are too large for me to keep on GitHub, and I can't figure out a place to put them for you to download directly to the HTCondor system. Go to that link and then right click and choose "Download" for each file to save them to your own computer (they're safe, I promise). Then in your HTCondor account:
```
$ cd ~/HTCondor_examples
$ cd image_processing
$ mkdir images
```
Now FTP the *.h5* image files from your computer to your new `~/HTCondor_examples/image_processing/images` subdirectory.

The datafiles in your new *images* subdirectory are Landsat 5 images of northeastern Minnesota for dates in 2010 and 2011 almost exactly 1 year apart. Within that time, a large forest fire occurred near the middle of the image. (The effects of another large, earlier forest fire are also visible in the northern part of both images.) We often examine the extent and severity of forest fires using calculated vegetation indices and before-and-after image differencing methods. Each of these datafiles contains metadata information and images that I have already processed from raw data to show surface reflectance measurements in six spectral bands (blue, green, red, near infrared, and two in the shortwave infrared range).

#### Part 2a: Short workflow

From these spectral band images we will calculate two vegetation indices, both of which are designed to indicate vegetation health and vigor. We'll do this first in a shorter example for the 2010 image, which will show before-fire forest conditions. Since the image is around the beginning of October, we'll see some effects of autumn leaf changes. We'll also see the condition of the forest areas that burned in 2006 and 2007 in the northern portion of the image.

The *workflow_short_dag.sub* contains information for four jobs ("nodes"). You'll see that each node is defined individually with a *JOB* line and then some variable values in a corresponding *VARS* line, just as in Part 1 above. In addition, the first nodes has a *SCRIPT POST* line that runs a small executable (paired *.sh* and *.py* scripts) to prepare data for the second node. The second node is a call to another *_dag.sub* file (actually generated by the first node's *SCRIPT POST* line) that has ten nodes of its own, thus the special *SUBDAG EXTERNAL* language. The third and fourth nodes each have a *SCRIPT PRE* line that runs a small executable (paired *.sh* and *.py* scripts) with different input values in each case.

`$ condor_submit_dag workflow_short_dag.sub`

In this case, Node A breaks up the image into ten chunks, then creates the *_dag.sub* file for Node B. The Node B *SUBDAG* calculates vegetation indices for these ten chunks in a distributed (parallel) schedule, since each does not depend on another. When all of Node B is completed, a script packages the resulting vegetation index chunks for use in Node C, which stitches the chunks back together into a full-size map and generates a figure for viewing. The two instances of Node C, one for each vegetation index, are scheduled to run separately with no dependence on each other.

You can check the job status soon after submission to see that Node A is queued:
```
$ condor_q $USER

OWNER    BATCH_NAME                         SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
xxxxxxxx workflow_short_dag.sub+XXXXXXXX  XX/XX XX:XX      _      _      1      4 XXXXXXXX.0
```

Then check a short time later, Node A is running:
```
$ condor_q $USER

OWNER    BATCH_NAME                         SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
xxxxxxxx workflow_short_dag.sub+XXXXXXXX  XX/XX XX:XX      _      1      _      4 XXXXXXXX.0
```

And another short time later, Node B with multiple SUBDAG nodes is queued:
```
$ condor_q $USER

OWNER    BATCH_NAME                         SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
xxxxxxxx workflow_short_dag.sub+XXXXXXXX  XX/XX XX:XX      1      _     10     14 XXXXXXXX.0 ... XXXXXXXX.0
```

And another short time later, those Node B elements are running (and some may be done already):
```
$ condor_q $USER

OWNER    BATCH_NAME                         SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
xxxxxxxx workflow_short_dag.sub+XXXXXXXX  XX/XX XX:XX      3      6      2     14 XXXXXXXX.0 ... XXXXXXXX.0
```

And another short time later, Node B is done and the two separate instances of Node C are queued or running:
```
$ condor_q $USER

OWNER    BATCH_NAME                         SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
xxxxxxxx workflow_short_dag.sub+XXXXXXXX  XX/XX XX:XX      2      1      1      4 XXXXXXXX.0 ... XXXXXXXX.0
```

And finally, if you catch it just before the job is completed:
```
$ condor_q $USER

OWNER    BATCH_NAME                         SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
xxxxxxxx workflow_short_dag.sub+XXXXXXXX  XX/XX XX:XX      3      1      _      4 XXXXXXXX.0 ... XXXXXXXX.0
```

Since I did not invoke any script(s) to manage (delete, move, etc.) all the resulting files, you'll find a bunch of intermediate *20101003_X_rows_XXXX-XXXX.h5* and *20101003_X.tar.gz* files that were created and used along the way. There will be many script *.err*, *.log*, and *.out* files. You'll also find the final output *20101003_ndii.h5* and *20101003_nbr.h5* datafiles and their corresponding *.png* figure plots. Your Normalized Burn Ratio (NBR) plot should look like this:

![20101003_nbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20101003_nbr.png)

The result for the Normalized Difference Infrared Index (NDII) is calculated from a different Landsat band combination and shows different values, but generally similar features, over the image area. Both NDII and NBR are oriented on indicating both greenness and moisture content of the vegetation. Healthy vegetation is both green and high in moisture content, with high NDII or NBR values. Dry vegetation (as during a drought period) will likely have lower moisture content but may still be green. Senescing vegetation (as during the autumn) may be changing color but could remain filled with moisture. Dead vegetation is neither green nor (typically) full of moisture, with low NDII or NBR values.

Was this faster, splitting the image into chunks for processing and then stitching the results back together, than processing the whole image at once? In this case, probably not: the vegetation index calculations were relatively quick operations, and really only about half of the image is processed since the area over Lake Superior is not considered. If you have a large image but are concerned only with specific areas or if your calculation steps take a longer time, perhaps with more analysis and decisions involved, splitting up your computational domain can be quite useful. This was just an example to show a workflow, not necessarily the fastest path to the results.

#### Part 2b: Long workflow

In the case of a forest fire, we may want an idea on the state of the forest before and after the event but we are also interested in the net *change* over the event. Was it a small fire that consumed ground litter but did not affect too many trees, or was it a large and intense fire that reduced entire stands to ash? In order to find this out, we use *change detection* methods. Specifically, we will use before-and-after differences to look at the fire's effects on the forest. This is exactly the method used by the US Forest Service to measure and evaluate forest fire impacts.

This longer workflow uses both repetition and then extension of the shorter workflow discussed above. We will do that entire shorter workflow twice, once on each of the images available in this example. We will then combine the results to calculate a difference map. We will calculate temporal differences for both NDII and NBR, though the NBR is more often used to evaluate forest fire impacts and has a specialized calculation that we will also perform.

If you want to remove the results of the short workflow example, do that now. Then,

`$ condor_submit_dag workflow_long_dag.sub`

The two images are processed in parallel to their respective vegetation indices, since that part of the workflow for each date does not depend on the results of the other. You should get the above NBR map and the corresponding map for the 2011 image, which looks like this:

![20111006_nbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20111006_nbr.png)

You can see clearly the size and shape of the burned area near the center of this after-fire image. After these vegetation index calculations are completed for both dates, we then use both maps for each vegetation index to calculate and plot the differences. The difference-NBR (dNBR) map looks like this:

![20101003_20111006_dnbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20101003_20111006_dnbr.png)

Our difference calculation takes into account a quirk of change detection using NBR, which often reverses the dates to give an idea of the severity of the change agent (the fire, in this case) which is the opposite of a value for the change in the vegetation health. (I myself prefer looking at NDII, but this is a better example of how fires are examined.) The US Forest Service also typically goes one step farther, adjusting the dNBR relative to its pre-fire value in case areas in the fire zone were already damaged from some earlier event. This is called the relative-dNBR (RdNBR) and for this fire event looks like this:

![20101003_20111006_rdnbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20101003_20111006_rdnbr.png)
