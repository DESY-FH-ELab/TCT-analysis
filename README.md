Content:
=======

1 - Requirements

2 - Installation instructions

3 - Test files, Test sensor data



1. Requirements
---------------

Minimal requirements are: a C++11 compliant compiler, ROOT v5.x or later (including v6.x), and CMake 2.6.x or later.


2. Installation instructions
----------------------------

  * Install CMake.
  * Install ROOT. 
  * Source ROOT 5.x before compilation. 
  * Go to /TCT-analysis/build
  * `> cd <.../TCT-analysis/>`
  * Let CMake create the makefile for you
  * `> cmake ..`
  * Now compile:
  * `> make install`

The program assumes the following folder structure for the data files to be read:
`<path-to-data>/<sample-name>/<temperature>/<bias-volt>/<data-file.txt>`


The program assumes the following folder structure for the results to be written:
`<path-to-results>`
here, a subfolder structure is created: `<sample-name>/<temperature>/<root-file.root>`

Look at the header of the sample files included in the repository.
If your headers differ from that, change the Read() method in acquisition.cc accordingly

  * Execute the binary: `> ./tct-analysis -r <path-to-project-directory>`
e.g.
`> ./tct-analysis -r /home/<user>/<my-path>/TCT-analysis/`
or if executed from /TCT-analysis/build
`> ./tct-analysis -r ..`



3. Test files, Test Sensor data
-------------------------------

The test data files are located at .../TCT-analysis/testdata/S57/295K/500V/.
The test sensor files are located at .../TCT-analysis/testsensor/.

