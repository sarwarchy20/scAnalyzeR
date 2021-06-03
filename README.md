#                                                              scAnalyzeR
**scAnalyzeR:** a comprehensive software package with graphical user interface for single-cell RNA sequencing analysis

# How to run the scAnalyzeR?
There are two different ways to setting up the pipeline in your own machine: 
#### Way 1: using from Docker image (strongly recommended)
1.	Download and install [Docker](https://www.docker.com/products/docker-desktop)
2.	Pull the docker image by running the following command:
  ```
   docker pull gscdocker/scanalyzer:1.0 
  ```
3.	Run the docker image locally on your computer and access the link:
  To run the docker image, execute the following command:
  ```
  docker run -d --rm -p 3838:3838 gscdocker/scanalyzer:1.0
  ```
  After running the docker image successfully, open the following link on a web browser (e.g., Firefox, Google Chrome) to access the pipeline: 
  ```
  http://localhost:3838/
  ```
#### Way 2: using from source
Firstly, you need to download and install following softwares (install R then RStudio):
Download and install R and RStudio on your machine,
i.	Download and install R (v-3.6.2 or above): https://cran.r-project.org/ 
ii.	Download and install RStudio: RStudio (v-1.1.456 or above): https://rstudio.com/products/rstudio/download/
After installing the R and RStudio on your machine successfully, then, you need to clone this (https://github.com/sarwar-chy/scAnalyzeR/archive/master.zip) repository as well as unzip it.
Now, please run the following script on RStudio to install the renv R package:
```
install.packages("renv")
```
Next, run the following scripts on RStudio to install all the dependent R packages. Please replace ~ with the location of your scAnalyzeR-master unzipped folder: 
```
renv::consent(provided=TRUE)
setwd("~/scAnalyzeR-master")
renv::restore() 
```
Finally, run the app using RStudio by running the script below:
shiny::runApp('~/scAnalyzeR-master/')
After successfully running the scAnalyzeR, the GUI will be displayed automatically.
<br/>
# Interface with user manual <br/>
https://github.com/sarwar-chy/scAnalyzeR/blob/master/user_manual/User_manual_scAnalyzeR.pdf





