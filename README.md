# Q_nonadiabatic

The package was tested on TACC's Stampede2 using 

* intel/18.0.2
* impi/18.0.2

To compile the code:

```
cd src
make QTRMPI=1
```

To run the program (Stampede2 SKX node):

```
module load intel/18.0.2
module load impi/18.0.2

export OMP_NUM_THREADS=1

time ibrun ./qtr_mpi ini.test
```

