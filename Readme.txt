Repression_Lag
This is a repository that contains Matlab scripts and a GUI interface that I have written to do image segmentation on confocal image stacks of early drosophila embryos stained in an in situ hybridization reaction.  The scripts do things like reliably segment nuclei from an image of fluorescently labeled DNA and filter and segment images of fluorescently labeled mRNA probes to  find sites of mRNA transcription. This can be done for multiple types of different mRNA and then one can find where they co-localize and determine which nuclei are expressing them. The code does this for thousands of cells at a time, and it is designed to be run in a semi-automated fashion allowing for the rapid analysis of hundreds of embryos.
The segmentation is implemented using two GUIs that split the different parts of the analysis up into a series of steps and allow for user input at each step, if necessary, or they can be run in an automated fashion. The first is imviewer4chan which takes a raw confocal stack and projects it onto a 2D image with multiple channels that can be segmented. EmbProcV3 then does the segmentation.
I hope that this will be of use for others working on the early drosophila embryo. I also hope that by putting it on GitHub this that will facilitate the collaborative development of code to solve common problems. Feel free to download the code, or better yet make your own account (which is free) and fork this repository. This way you and go crazy and change and add things for the better that can be seen and used by others. I am also more than happy to help troubleshoot problems with the code and will be actively involved in further development. 
I developed a lot of this code for the analysis reported in my recent publication:
Bothma et al., The Snail Repressor Inhibits Release, Not Elongation, of Paused Pol II in the Drosophila Embryo, Current Biology (2011), doi:10.1016/j.cub.2011.08.019

If you find any of the code useful in your own work please cite it.

Jacques Bothma
UC Berkeley

jpbothma@gmail.com


