#!/bin/bash

c++ -fPIC -shared -o ginacg.so ginacg.cpp -lcln -lginac --std=c++11