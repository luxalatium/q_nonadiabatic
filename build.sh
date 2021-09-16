ml reset
ml intel/18.0.2
ml impi/18.0.2

make clean-all
make QTRMPI=1 NADBPOT_TN1=1
mv qtr_mpi qtr_mpi_tn1

#make clean-all
#make QTRMPI=1 NADBPOT_TN2=1
#mv qtr_mpi qtr_mpi_tn2

#make clean-all
#make QTRMPI=1 NADBPOT_GM3=1
#mv qtr_mpi qtr_mpi_gm3

#make clean-all
#make QTRMPI=1 NADBPOT_MM4=1
#mv qtr_mpi qtr_mpi_mm4

#make clean-all
#make QTRMPI=1 NADBPOT_GM5=1
#mv qtr_mpi qtr_mpi_gm5

#make clean-all
#make QTRMPI=1 NADBPOT_PBC=1
#mv qtr_mpi qtr_mpi_pbc_2

