#                                                              scAnalyzeR
**scAnalyzeR:** a comprehensive software package with graphical user interface for single-cell RNA sequencing analysis

# How to run the scAnalyzeR?
Firstly, you need to download and install the following softwares(install R then RStudio):
###### Download and install R and RStudio on your machine,
`1. Download and install R (v-3.6.1 or above):` https://cran.r-project.org/ <br/>
`2. Download and install RStudio: RStudio (v-1.1.456 or above):` https://rstudio.com/products/rstudio/ 

After installed the R and RStudio on your machine , then you also need to install the shiny package. 
Please run the following code on RStudio to install the shiny package: <br/>

`install.packages("shiny")` <br/>
Now, your machine is ready for running the scAnalyzeR app.<br/>

## There are many ways to run the scAnalyzeR:<br>
**Easiest way is to use runGitHub**<br/>
Run the following code of line on RStudio.<br/>
`shiny::runGitHub("scAnalyzeR", "sarwar-chy")` <br/>
or <br/>
**Run a tar or zip file directly** <br/>
`shiny::runUrl("https://github.com/sarwar-chy/scAnalyzeR/archive/master.tar.gz")` <br/>
`shiny::runUrl("https://github.com/sarwar-chy/scAnalyzeR/archive/master.zip")` <br/>

**If you want to use it locally, you can simply clone this(https://github.com/sarwar-chy/scAnalyzeR/archive/master.zip) repository as well as unzip it and run the app using RStudio.** <br/>
To run the app on your local computer, use RStudio to run the following code and replace ~ with the location of your scAnalyzeR-master folder.<br/>
`shiny::runApp('~/scAnalyzeR-master/')`
<br/>
# Interface with user manual <br/>
https://github.com/sarwar-chy/scAnalyzeR/blob/master/user_manual/User_manual_scAnaluzeR.pdf





