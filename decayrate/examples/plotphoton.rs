use clap::Parser;
use decayrate::gamma_a_photon;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Calculate relative change in decay rate of a n-particle decay with given parameters for multiple accelerations
#[derive(Parser, Debug, Clone)]
struct Args {
    /// Mass difference between decaying and heaviest particle
    #[arg(long, default_value_t = 1.5)]
    delta: f64,

    /// Minimum acceleration
    #[arg(long, default_value_t = 0.03)]
    mina: f64,

    /// Maximum acceleration
    #[arg(long, default_value_t = 3.0)]
    maxa: f64,

    /// Number of accelerations to consider between mina and maxa
    #[arg(long, default_value_t = 128)]
    na: usize,
}

fn main() {
    let args = Args::parse();

    let result: Vec<_> = (0..args.na)
        .into_par_iter()
        .map(|i| args.mina * ((args.maxa / args.mina).powf((i as f64) / ((args.na - 1) as f64))))
        .map(|a| {
            (
                a,
                gamma_a_photon(a, args.delta) / gamma_a_photon(0.0, args.delta),
            )
        })
        .collect();

    println!("{:?}", result);
}
