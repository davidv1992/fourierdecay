This repository contains the code used to create the plots in ["Unruh detectors, Feynman diagrams, acceleration and decay" by Wim Beenakker and David Venhoek](). It depends on the GMP, MPFR and MPFI libraries, which should be installed before running the code.

To recreate figure 1, run
```
cargo run --example plota -- --na 800 --n 1000001 --maxa 30 | tee plotdata.txt
```
then plot the data with
```
python3 gen_plot_single.py
```

To recreate figure 2, run
```
cargo run --example plota -- --na 800 --n 1000001 -d 2 | tee plotdata1.txt
cargo run --example check2d -- --na 800 | tee plotdata2.txt
```
then plot the data with
```
python3 gen_plot_compare.py
```

The code for calculating the bessel functions in the archt/archt folder is originaly by Andrew R. Booker, Andreas Strombergsson and Holger Then. When used in preparing a paper, it should be cited as:

A. R. Booker, A. Strombergsson, and H. Then,
Bounds and algorithms for the K-Bessel function of imaginary order,
LMS J. Comp. Math. 16 (2013) 78--108.
