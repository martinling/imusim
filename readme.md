Welcome to the Dockerized version of IMUSim!

All necessary dependencies have been downloaded using either the aptitude or pip package managers.
Due to display issues with Docker, data visualization and plotting do not currently work as intended.
However, the dependencies necessary for Mayavi are present in the container, meaning that once a working x-server is established, IMUSim should be fully functional.
At present, basic functionality has been established and the data from the GitHub repository is copied into the folder on build rather than being cloned to account for changes made to the code that fix compatibility with Vicon CSV files.

The Dockerfile included in this version of the IMUSim image will:
-Pull a base Ubuntu image
-Copy files from the current directory into the /tmp folder
-Set /tmp as the working directory in the container
-Copy the file "test.csv" to the directory /imusim-py3/imusim/maths
-Install all necessary Linux packages (namely, python, pip, and several Mayavi dependencies)
-Set an environment variable to solve an otherwise obscure Unicode error
-Specify that the frontend is uninteractive (not necessary, but reduces the amount of red text on build as Docker will stop trying to find a frontend that is not there)
-Use pip to install and build all IMUSim dependencies other than Mayavi into the /tmp directory
-Use pip to install and build Mayavi into the /tmp directory
-Set the working directory to /imusim-py3, which contains all of the code from the GitHub repository that has been edited as necessary
-Manually build the missing C files for the imusim/maths files using cython, specifying the language level as 3
-Run IMUSim's setup.py file
-Set the working directory to imusim/maths (python must be run from this directory to avoid ModuleNotFound and Import errors)
-Create a folder called mount

The purpose of this image is to process CSV files containing Vicon capture data and use them to simulate an IMU. Thus, the volume tag should be used when the container is run so that the container has access to all necessary files from the host machine. For information on mounting and volumes, go to https://docs.docker.com/v17.09/engine/admin/volumes/bind-mounts/

The file "test.csv" will be included as a template for appropriate CSV format. Any CSV file to be processed by IMUSim should match the format of this file, which can be opened using the nano text editor from the Linux terminal in the container (the command is nano "test.csv").
