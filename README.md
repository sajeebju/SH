# SH
SH wave propagation modeling using staggered-grid Finite Difference Method and seismogram calculation

Compile just the sh2d.cpp\
clang++ -o sh2d sh2d.cpp -std=c++11 

Compile for seismogram calculation\
g++ -std=c++11 -Wall -Wextra -c sh2dwave.cpp -o sh2dwave.o \
g++ -std=c++11 -Wall -Wextra -c seismo.cpp -o seismo.o \
g++ sh2dwave.o seismo.o -o sh2dwave_program \
./sh2dwave_program


![Image](https://github.com/user-attachments/assets/2067a9c0-05d4-439a-8ce8-39fef89e145e)
