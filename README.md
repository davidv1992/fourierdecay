This repository contains the code used to create the plots in TODO. It depends on the GMP, MPFR and MPFI libraries, which should be installed before running the code.

To recreate figure 1, run
```
cargo run --example plota -- --na 800 --n <TODO> --mina 0.00004 --maxa 0.3 --delta 1.00145 -m 1 --precision 256 | tee fig1.dat
```
then plot the data with
```
./gen_plot_single.py --amin 0.00004 --amax 0.3 --rmax 3.0 fig1.dat 1.698e33 -o fig1.png
```

To recreate figure 2, run
```
cargo run --example plotweak -- --na 800 --n <TODO> --ni <TODO> --mina 0.003 --maxa 0.3 --delta 1.036 --precision 256 | tee fig2.dat
```
then plot the data with
```
./gen_plot_single.py --amin 0.003 --maxa 0.3 fig2.dat 2.327e29 -o fig2.png
```

To recreate figure3, run
```
cargo run --example plotphoton -- --na 800 --delta 1 | tee fig3.dat
```
then plot the data with
```
./gen_plot_single.py --rmax 10 fig3.dat 3.805e24 -o fig3.png
```

The code for calculating the bessel functions in the archt/archt folder is originaly by Andrew R. Booker, Andreas Strombergsson and Holger Then. When used in preparing a paper, it should be cited as:

A. R. Booker, A. Strombergsson, and H. Then,
Bounds and algorithms for the K-Bessel function of imaginary order,
LMS J. Comp. Math. 16 (2013) 78--108.
