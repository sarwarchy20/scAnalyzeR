#                                                              scAnalyzeR
**scAnalyzeR:** a comprehensive software package with graphical user interface for single-cell RNA sequencing analysis

# How to run the scAnalyzeR?
Firstly, you need to download and install the following softwares:
###### A. if you are a windows user,
`1. Download and install R (v-3.6.1 or above): https://cran.r-project.org/bin/windows/base/`
`2. Download and install RStudio: RStudio (v-1.1.456 or above): https://rstudio.com/products/rstudio/download/#download `
###### B. if you are a Linux user,
`1. Download and install R (v-3.6.1 or above): https://cran.r-project.org/bin/linux/
2. Download and install RStudio: RStudio (v-1.1.456 or above): https://rstudio.com/products/rstudio/download/#download`
###### C. if you are a MacOS user,
`1. Download and install R (v-3.6.1 or above): https://cran.r-project.org/bin/macosx/
2. Download and install RStudio: RStudio (v-1.1.456 or above):https://rstudio.com/products/rstudio/download/#download`

After installed the R and RStudio on your machine , then you also need to install the shiny package. 
Please run the following code on RStudio to install the shiny package:

`install.packages("shiny")`
Now, your machine is ready for running the scAnalyzeR app.

## There are many ways to run the scAnalyzeR:
**Easiest way is to use runGitHub**
Run the following code of line on RStudio.
`shiny::runGitHub("scAnalyzeR", "sarwar-chy")`

or 
# Run a tar or zip file directly
`shiny::runUrl("https://github.com/sarwar-chy/scAnalyzeR/archive/master.tar.gz")
shiny::runUrl("https://github.com/sarwar-chy/scAnalyzeR/archive/master.zip")
`
**If you want to use it locally, you can simply clone this repository and run the app using RStudio.**
To run the app on your local computer, use RStudio to run the following code and replace ~ with the location of your scAnalyzeR folder.

`shiny::runApp("~/scAnalyzeR")`

# Interface with user manual






