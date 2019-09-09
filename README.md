# MCSFD

Install

Install OpenMP fisrt

    sudo apt-get install libomp-dev

Then clone Eigen to ./external
    
    cd external 
    git clone https://github.com/eigenteam/eigen-git-mirror.git

If you are using Windows, download glm headers folder to ./include.

If you are using linux, just

    sudo apt-get install libglm-dev

Then you can just do

    ./run.sh c x

    or

    mkdir build
    cd build
    cmake ..
    make

    cd ..
    ./bin/main

If your system is Windows, then
use 
    
    cmake -G"Visual Studio 14 2015 Win64" .. 

Switch to other examples by using set(MAINFUNC examplexxx.cpp) in ./CMakeLists.txt

If you use visual studio, you need copy 2tet.ele, 2tet.node and all xml files to visual studio project root.