## Notes on problems I encounted with XIOS3

Using XIOS3 is not plug-and-play.

First: xmllint is a useful tool for eliminating errors with parsing .xml files.

Second:
You need to ensure that your tasks match the resources requested (match the nodes of a job with the number of xios procs). The XIOS nprocs assigned to the job also needs to scale with the amount of output that is requested in the .xml files. 
A common issue with XIOS3 is the job hanging but no error in the ocean.output, restart files are written (not using xios) but no data files. No error messages.
If no xios\_client\_\*.out files are made xios3 isn't running.
It may be that XIOS is getting killed too soon before running down and saving output. The hanging can be cause by memory issues from not assigning enough processors to either xios nprocs or not enough nodes to the job your submitting.
For exmaple with NAARC, setting nprocs in iodef.xml to 42 wasn't enought for the month 3d output requested and there was an issue with no output being written.

The servers in the submitted job (S num\_servers) needs to match nprocs assigned for xios3. Then the mapping (C) can be worked out by running the job and checking the log files if number of nodes is high enough and nproc is high enough. 
===>>> : W A R N I N G
 
          ===============
 
    The number of mpi processes:  6131
    exceeds the maximum number of ocean subdomains =  6122 
Or, the brute force solution is to count the 0's in the map and match that to the tasks allocated in iodef.xml

Third:
There may be a issue with slow start up time of XIOS3.
If you are getting output files but it takes a long time to write output. Check timing\_step.nc to see if there is a large lag at the write time between files. 

If the first time step take a long time it could be that the chunks in the input netcdf files (intial conditions, ancillaries, forcing) is not set up to match the number of cores assigned to nemo. 
If you have ~6200 cores trying to read the data they will all read as much data as it takes to fill the arrays. Tiles might be around 24 x 23 in size, but "initcd\_m01.nc" file has chunks that are: so:\_ChunkSizes = 10, 515, 618 ;
which means that every tile is going to be reading in a huge amount of data to memory. Check with ncdump -h -s. Rechunk the data.

Forth:
Slow wiritng output. You can also set the  sync_freq="1d"  to be longer, e.g.  sync_freq="15d", sync_freq="1mo"  or  sync_freq="45d". Whilst the job is running, sync_freq forces it to pause until the tasks have all caught up. You should see performance gain for a longer sync_freq, at the slight expense of increased memory usage. The only downside is that if your model crashes not all your data will have been written to disk - for example if you are writing daily output and set sync_freq = "1mo", if your model crashes on day 25 you may not have written output for days 1-24.
