# SH
SH wave propagation modeling using staggered-grid Finite Difference Method and seismogram calculation

Compile for seismogram calculation\
g++ -std=c++11 -Wall -Wextra -c sh2dwave.cpp -o sh2dwave.o \
g++ -std=c++11 -Wall -Wextra -c seismo.cpp -o seismo.o \
g++ sh2dwave.o seismo.o -o sh2dseis_program \
./sh2dseis_program


SH Wave Propagation using absorbing boundary condition

![Image](https://github.com/user-attachments/assets/d9790660-916e-4224-ace6-45fc14db45f2)
