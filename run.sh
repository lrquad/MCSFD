cmake_if=${1:-"false"}
stoprun=${2:-"false"}

if [ $cmake_if = "r" ]; then
    rm -r build/
    mkdir build
fi

cd build

if [ "$cmake_if" = "r" ] || [ "$cmake_if" = "c" ]; then
    cmake -DCMAKE_BUILD_TYPE=Release ..
fi

make
cd ..

if [ "$stoprun" = "x" ]; then
catchsegv ./bin/main -o hai
fi