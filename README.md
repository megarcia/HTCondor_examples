# HTC_examples

## Part 0: Setup

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

### Part 1a: CSV processing example

This example illustrates the basic HTCondor job unit.

`$ cd csv_processing`

The data file *GHCND\_WIS\_2015.csv* is named to indicate source, geographic coverage, and temporal coverage. These are [Global Historic Climate Network-Daily (GHCND)](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn) surface weather observations for all 1st-order and cooperative stations in the rectangle containing Wisconsin for 2015. The dataset has a header line with column names, and includes a number of fields such as location, measurements, and data quality flags. I wrote the Python script *clean\_GHCND.py* to "clean" this raw dataset, removing flagged measurements and stations without location information, etc. This is the same procedure as *process\_NCEI\_00.py* in my [WxCD project](https://megarcia.github.io/WxCD).

`$ condor_submit clean_GHCND.sub`

You can check the job status with

`$ condor_q $USER`

When the job is completed, you'll find a lot of processing (cleaning) information in the *clean\_GHCND.out* file. The original input datafile is left intact, and you'll find new *.csv* files with various info on the errors encountered and the cleaned dataset. The results can be further used for temporal examination, spatial analysis, etc.

### Part 1b: Thinking ahead

Let's say you have another year of data in *GHCND\_WIS\_2016.csv* that you want to process in the same way. You have a few ways forward here.

The easiest option might be to edit the *clean\_GHCND.sub* file with the new file name, then run the job as we did for the 2015 file. At the end of the job you'll receive new *.csv* output files for the 2016 dataset, as desired, but you'll also get back a new *clean\_GHCND.out* file that will overwrite your 2015 *clean\_GHCND.out* file that you may have wanted to keep as a record of the processing decisions. You can manually rename the old one to *clean\_GHCND\_2015.out* before running the new one, and then manually rename the new output as well. For a couple of files, this gets the job done, but what if you have 30 years of data in separate files?

The next-easiest option might be to name those output files as desired in *clean\_GHCND.sub* itself. You can edit the file names that are delivered as job *log*, *error*, and *output* on their respective lines so that the names are specific to the dataset, e.g. *output = clean\_GHCND\_2015.out*. You might want to rename the submit file to *clean\_GHCND\_2015.sub* given how specific it is now on the inside. Then you can copy it to *clean\_GHCND\_2016.sub*, edit that new submit file everywhere the date appears, and run your 2016 dataset as its own separate unit. You could then repeat this copy/edit/submit process 28 more times for the remainder of your 30-year dataset.

Instead of all that extra work, the HTCondor system with DAGMan capabilities lets us use variables in submit files. If the only thing that really changes from one job to the next is the year, and each year of input data is contained in its own *.csv* file, we might make that our submit file variable with a value that is assigned elsewhere. In this case, edit *clean\_GHCND.sub* so that *output = clean\_GHCND\_\$(year).out* and the log and error lines are changed likewise. Edit the *arguments* and *transfer\_input\_files* lines to use that variable as well. Your *clean\_GHCND.sub* file should end up looking like this:

```
# UW-Madison HTCondor submit file
# clean_GHCND.sub
universe = vanilla
log = clean_GHCND_$(year).log
error = clean_GHCND_$(year).err
executable = clean_GHCND.sh
arguments = GHCND_WIS_$(year).csv
output = clean_GHCND_$(year).out
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = ../python3.tar.gz,clean_GHCND.py,GHCND_WIS_$(year).csv
request_cpus = 1
request_memory = 8GB
request_disk = 8GB
requirements = (OpSys == "LINUX") && (OpSysMajorVer == 7)
queue 1
```

Note that the names of *clean\_GHCND.sh* and *clean\_GHCND.py* remain unchanged, since we want to use those scripts as they're already written. Now we need to assign the values for the *year* variable, for which we create a DAG file. At its most basic, to run the single job that we processed above, we could create a file *clean\_GHCND\_dag.sub* with only two lines:

```
JOB A_2015 clean_GHCND.sub
VARS A_2015 year="2015"
```

Note first that the *JOB* declaration must come first, that the *VARS* declaration must have a job name (*A\_2015*) that matches the *JOB* line, and that the variable name itself must match what you're using in the *.sub* file that you declared on the *JOB* line. Note that the variable value must be in quotes, whether it's a number or an alphanumeric string (such as a file name or a directory path).

At this point you might want to clean up your directory with

`$ rm *.err *.log *.out`

Then type a slightly different command:

`$ condor_submit_dag clean_GHCND_dag.sub`

When the job is complete, we should get the same output that we received the first time through, but with the output files named a little more specifically, as we'll want when processing more than one year of data. We could add another job to *clean\_GHCND\_dag.sub*:

```
JOB A_2015 clean_GHCND.sub
VARS A_2015 year="2015"

JOB A_2016 clean_GHCND.sub
VARS A_2016 year="2016"
```

If you had that *GHCND\_WIS\_2016.csv* datafile in place, you could submit your new two-job DAG file with the same command as above. Since the two data files and processes wouldn't depend on each other, they could run simultaneously, and you'd get individualized output from each job as it's completed. From this, to process all 30 years of the hypothetical dataset, you can edit *clean\_GHCND\_dag.sub* accordingly and get all of the processing done in the same time it would take to process just one year of data.

## Part 2: Image processing example

This example illustrates running HTCondor workflows using these DAGMan capabilities.

Before you start, you need to get the two Landsat image files to be processed. They are too large (~100 MB each) for me to keep in this project on GitHub, and I don't yet have a place on the web to put them for you to download directly to the HTCondor system. First, in your HTCondor account:

```
$ cd ~/HTCondor_examples
$ cd image_processing
$ mkdir images
```

The next step depends on where you are working. **If you are working on one of the UW-Madison CHTC submit nodes (including *townsend-submit*)**, you are welcome to copy the images from my own *squid* directory:

`$ cp /squid/megarcia/HTCondor_examples/images/* images/`

**If you are working on an HTCondor system somewhere else**, open a tab in your browser to my [squid directory](http://proxy.chtc.wisc.edu/SQUID/megarcia/HTCondor_examples/images/) and right click and choose "Download" for each file to save them to your own computer (they're safe, I promise). Then SFTP or otherwise transfer the two *.h5* image files from your computer to your new `~/HTCondor_examples/image_processing/images` subdirectory.

The datafiles now in your new *images* subdirectory are Landsat 5 images of northeastern Minnesota for dates in 2010 and 2011 almost exactly 1 year apart. Within that time, a large forest fire occurred near the middle of the image. We often examine the extent and severity of forest fires using calculated vegetation indices and before-and-after image differencing methods. Each of these datafiles contains metadata information and images that I have already processed from raw data to surface reflectance values in six spectral bands: blue, green, red, near infrared (NIR), and two in the shortwave infrared (SWIR) range. We see RGB digital images all the time, and can generate them from the first three bands in the Landsat dataset. Part of the power of Landsat data is really in the NIR and SWIR bands, beyond human vision but very useful for looking at forests and other land cover. Landsat sensors also provide a couple additional spectral bands, including one in the thermal infrared (TIR) range that is useful for calculating surface temperatures and for detecting even small forest fires in remote areas.

### Part 2a: Short workflow

From these spectral band images we will calculate two vegetation indices, both of which were formulated using NIR and SWIR bands to indicate vegetation health and vigor. We'll do this first in a shorter example for just the 2010 image, which will show before-fire forest conditions. Since the image is around the beginning of October, we'll see some effects of autumn leaf changes. We'll also see the condition of the forest areas that burned in 2006 and 2007 in the northern portion of the image.

The *workflow\_short\_dag.sub* contains information for four jobs ("nodes"). You'll see that each node is defined individually with a *JOB* line and then some variable values in a corresponding *VARS* line, as in Example 1b above. In addition, the first node has a *SCRIPT POST* line that runs a small executable (paired *.sh* and *.py* scripts) to prepare data for the second node. The second node is a call to another *\_dag.sub* file (which was actually generated by the first node's *SCRIPT POST* line) that has ten nodes of its own, so it's named with the special *SUBDAG EXTERNAL* syntax. The third and fourth nodes each have a *SCRIPT PRE* line that runs a small executable (another pair of *.sh* and *.py* scripts) with different input values in each case. The *PARENT ... CHILD ...* lines tell the system the order in which the nodes are to be run. If you look at the contents of the *SUBDAG* submit file after it's generated by the system, you'll see that it has no *PARENT ... CHILD ...* lines, since each of its nodes can run independently and therefore simultaneously.

`$ condor_submit_dag workflow_short_dag.sub`

In this case, Node A breaks up the image into ten chunks, then creates the *\_dag.sub* file for Node B. The Node B *SUBDAG* calculates vegetation indices for these ten chunks in a distributed (parallel) schedule, since each does not depend on another. When all of Node B is completed, another script packages the resulting vegetation index chunks for use in Node C, which stitches the chunks back together into a full-sized map and generates a figure for viewing. The two instances of Node C, one for each vegetation index, are also scheduled to run separately with no dependence on each other.

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

Since I did not invoke any script(s) to manage (delete, move, etc.) all the resulting files, you'll find a bunch of intermediate *20101003\_X\_rows\_XXXX-XXXX.h5* and *20101003\_X.tar.gz* files that were created and used along the way. There will be many script *.err*, *.log*, and *.out* files. You'll also find the final output *20101003\_ndii.h5* and *20101003\_nbr.h5* datafiles and their corresponding *.png* figure plots. Your Normalized Burn Ratio (NBR) plot should look like this:

![20101003_nbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20101003_nbr.png)

You can see the effects of a couple of earlier forest fires in the northern part of this image. The result for the Normalized Difference Infrared Index (NDII) is calculated from a different Landsat band combination and shows different values, but generally similar features, over the image area. Both NDII and NBR are oriented on indicating both greenness (using the NIR spectral band) and moisture content (using one of the SWIR spectral bands) of the vegetation. Healthy vegetation is both green and high in moisture content, with high NDII and NBR values. Dry vegetation (as during a drought period) will likely have lower moisture content but may still be green. Senescing vegetation (as during the autumn) may be changing color but could remain filled with moisture. Dead vegetation is neither green nor (typically) full of moisture, with low NDII and NBR values.

Was this workflow, splitting the image into chunks for processing and then stitching the results back together, faster than processing the whole image at once? In this case, probably not: most vegetation index calculations are relatively quick operations, for which Python is optimized to handle for whole arrays at once. If you have a large image but are concerned only with specific areas, or if your calculation steps take a longer time with more analysis and decisions involved, splitting up your computational domain can be quite useful. This was just an example to show a processing workflow, not necessarily the fastest path to this particular result. If you have 200 images for this location (as I do) with lots of vegetation indices and annual combinations to examine, building a workflow to automate the processing will save you a good amount of time in the end.

### Part 2b: Long workflow

In the case of a forest fire, we may want an idea on the condition of the forest before and after the event, but we are also interested in the forest *change* over the event. Was it a small fire that consumed ground litter but did not affect too many trees, or was it a large and intense fire that reduced entire forest stands to ash? In order to find this out, we use *change detection* methods. Specifically, we will use before-and-after differences to look at the fire's effects on the forest. This is exactly the method used by the US Forest Service to measure and evaluate forest fire severity and impacts.

This longer workflow uses both repetition and then extension of the shorter workflow in Example 2a. We will do that entire shorter workflow twice, once for each of the images you're given. We will then combine the results to calculate a difference map for each vegetation index. We will calculate temporal differences for both NDII and NBR, though the NBR is more often used to evaluate forest fire impacts and has a specialized calculation that we will also perform.

If you want to remove the results of the short workflow example, do that now (if you don't, this step will simply overwrite them). Then,

`$ condor_submit_dag workflow_long_dag.sub`

The two images are processed in parallel to their respective vegetation indices, since that part of the workflow for each date does not depend on the results of the other. You should get the above 2010 NBR map and a corresponding 2011 NBR map, which looks like this:

![20111006_nbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20111006_nbr.png)

You can see the size and shape of the burned area near the center of this after-fire image. After these vegetation index calculations are completed, we can use the dates together to calculate and examine the change over time. The difference-NBR (dNBR) map looks like this:

![20101003_20111006_dnbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20101003_20111006_dnbr.png)

Our difference calculation takes into account a quirk of change detection using NBR, which often reverses the dates to give an idea of the severity of the change agent (the fire, in this case) which is the opposite of a value for the change in the vegetation health. (I prefer looking at NDII, but this is a better example of how fires are examined in practice.) The US Forest Service also typically goes one step farther, adjusting the dNBR relative to its pre-fire value, in case areas in the fire zone were already damaged from some earlier event. This is called the relative-dNBR (RdNBR) and for this fire event looks like this:

![20101003_20111006_rdnbr.png](https://github.com/megarcia/HTCondor_examples/blob/master/image_processing/results/20101003_20111006_rdnbr.png)

Most of the image has values near 0, indicating very little change in the intervening year. In this event, the [Pagami Creek Fire in August 2011](http://www.fs.usda.gov/detail/superior/home/?cid=stelprdb5341928) turned intense and the damage to the affected forest area was generally severe. High RdNBR values cover much of the fire area except for the northwestern lobe, where the fire was started by lightning and left more moderate forest damage. The early fire was fought somewhat successfully until high winds [carried it southeast and then rapidly east over just a few days](http://inciweb.nwcg.gov/incident/map/2534/7/).
