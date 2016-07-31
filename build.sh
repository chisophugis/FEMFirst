#clang++ -g -std=c++11 -fsanitize=address main.cpp -fcolor-diagnostics
#clang++ -g -std=c++11 -fsanitize=undefined main.cpp -fcolor-diagnostics
clang++ -O3 -g -std=c++11 main.cpp -fcolor-diagnostics
