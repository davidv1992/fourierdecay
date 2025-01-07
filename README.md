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
```