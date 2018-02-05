#!/bin/bash
CC=$1
echo "Compiling ginac_mg5_interface using the compiler $CC"
echo "$> $CC -fPIC -shared -o ginac_mg5_interface.so ginac_mg5_interface.cpp -lcln -lginac --std=c++11"
$CC -fPIC -shared -o ginac_mg5_interface.so ginac_mg5_interface.cpp -lcln -lginac --std=c++11