MPIWRAP_DEBUG=verbose LD_PRELOAD=/home/jtu/valgrind/lib/valgrind/libmpiwrap-amd64-linux.so \
/home/jtu/mpich3.2.1/bin/mpiexec -n 4 /home/jtu/valgrind/bin/valgrind \
 --tool=memcheck --trace-children=yes --leak-check=full --show-reachable=no --track-origins=yes \
./dist/Debug/iditm2d -snes_type newtontr -ksp_type fgmres -pc_type asm \
-fp_trap -snes_rtol 1.0e-10 -snes_stol 1.0e-10 -snes_atol 1.0e-10 \
-kps_atol 1.0e-11 -kps_rtol 1.0e-11 -ksp_stol 1.0e-11
