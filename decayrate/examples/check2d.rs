use archt::Ctx;
use clap::Parser;
use decayrate::{gamma0, gammma_a_d2};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Calculate relative change in decay rate of a 2-particle decay with given parameters for multiple accelerations in 2d, using the closed form
#[derive(Parser, Debug, Copy, Clone)]
struct Args {
    /// Mass difference between decaying and heaviest particle
    #[arg(long, default_value_t = 1.5)]
    delta: f64,

    /// Mass of lightest particle
    #[arg(short, default_value_t = 1.0)]
    m: f64,

    /// Minimum acceleration
    #[arg(long, default_value_t = 0.03)]
    mina: f64,

    /// Maximum acceleration
    #[arg(long, default_value_t = 3.0)]
    maxa: f64,

    /// Number of accelerations to consider between 0 and maxa
    #[arg(long, default_value_t = 128)]
    na: usize,

    /// Precision of bessel function calculations
    #[arg(long, default_value_t = 53)]
    precision: u32,
}

fn main() {
    let args = Args::parse();

    let base = gamma0(args.delta, args.m, 2.0);
    let result: Vec<_> = (0..args.na)
        .into_par_iter()
        .map(|i| args.mina * ((args.maxa / args.mina).powf((i as f64) / ((args.na - 1) as f64))))
        .map(|a| {
            (
                a,
                gammma_a_d2(&mut Ctx::new(args.precision), a, args.delta, args.m) / base,
            )
        })
        .collect();
    println!("{:?}", result);
}
