/home/jtu/mpich3.2.1/bin/mpiexec -n 4 ./dist/Release/iditm2d \
-snes_type newtontr -ksp_typ fgmres -pc_type asm -snes_mgonitor \
-snes_rtol 1.0e-10 -snes_stol 1.0e-10 -snes_atol 1.0e-10  \
-ksp_rtol 1.0e-11 -ksp_atol 1.0e-11 -snes_trtol 1.0e-12
