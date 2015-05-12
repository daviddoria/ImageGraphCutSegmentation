Getting the code
----------------
After you have cloned this repository, you will need to initialize the submodules:
git submodule update --init --recursive

If you ever update with 'git pull origin master', you must then do 'git submodule update --recursive'.

This pulls in all of the dependencies including Mask (which includes ITKHelpers and then Helpers).

Overview
--------
This software allows the user to perform a foreground/background segmentation of an image.
This implementation is based on "Graph Cuts and Efficient N-D Image Segmentation" by Yuri Boykov (IJCV 2006).

Build notes
------------------
This code depends on c++0x/11 additions to the c++ language. For Linux, this means it must be built with the flag
gnu++0x (or gnu++11 for gcc >= 4.7).

Dependencies
------------
- ITK >= 4
- Boost 1.51 
You can tell this project's CMake to use a local boost build with: cmake . -DBOOST_ROOT=/home/doriad/build/boost_1_51
