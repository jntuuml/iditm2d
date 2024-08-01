/home/jtu/mpich3.2.1/bin/mpiexec -n 4 ./dist/Debug/iditm2d \
-snes_type newtontr -ksp_type fgmres -pc_type asm -snes_monitor \
-fp_trap -snes_rtol 1.0e-10 -snes_stol 1.0e-10 -snes_atol 1.0e-10  \
-ksp_rtol 1.0e-11 -ksp_atol 1.0e-11 -snes_trtol 1.0e-10 -start_in_debugger
