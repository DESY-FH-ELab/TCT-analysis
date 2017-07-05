# TCT-Analysis

### Content:

1 - Requirements

2 - Installation instructions

3 - Test files, Test sensor data


### 1. Requirements

Minimal requirements are: a C++11 compliant compiler, ROOT v5.x or later (including v6.x), and CMake 2.6.x or later. Also Qt4 or Qt5, in case of GUI compilation.


### 2. Installation instructions

Following text is the short version of the compilation guide, in case of problems, please check the [Developers Guide](https://github.com/DESY-FH-ELab/TCT-analysis/blob/master/developers_guide/tctanalysis_developers_guide.pdf).

  * Install CMake.
  * IF GUI needed - install Qt4 or Qt5
  * Install ROOT. 
  * Source ROOT 5.x before compilation. 
  * IF LeCroyRAW data converter needed - put the external LeCroyConverter lib to the external/LeCroyConverter/lib/libLeCroy.so
  * Go to /TCT-analysis/build
  * `> cd <.../TCT-analysis/>`
  * Let CMake create the makefile for you (available optionds -DWITH_GUI=ON, -DWITH_LECROY_RAW=ON)
  * `> cmake ..`
  * Now compile:
  * `> make install`

The oscilloscope data analysis part of program assumes the following folder structure for the data files to be read:
`<path-to-data>/<sample-name>/<temperature>/<bias-volt>/<data-file.txt>`


The program assumes the following folder structure for the results to be written:
`<path-to-results>`
here, a subfolder structure is created: `<sample-name>/<temperature>/<root-file.root>`

Look at the header of the sample files included in the repository.
If your headers differ from that, change the Read() method in acquisition.cc accordingly

  * (console version) Execute the binary from the TCT-analysis/build folder: `> ./tct-analysis -af <path-to-analysis-file>`
e.g.
 ./tct-analysis -af ../testanalysis/test_ana.txt




### 3. Test files, Test Sensor data

The test data files are located at .../TCT-analysis/testdata/S57/295K/500V/.
The test data file in TCT data format is located at .../TCT-analysis/testdata/lpnhe/.
The test sensor files are located at .../TCT-analysis/testsensor/.

