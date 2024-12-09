use archt::Ctx;
use clap::Parser;
use decayrate::{gamma0, gamma_a};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Calculate relative change in decay rate of a 2-particle decay with given parameters for multiple accelerations
#[derive(Parser, Debug, Copy, Clone)]
struct Args {
    /// Mass difference between decaying and heaviest particle
    #[arg(long, default_value_t=1.5)]
    delta: f64,

    /// Mass of lightest particle
    #[arg(short, default_value_t=1.0)]
    m: f64,

    /// Number of space-time dimensions
    #[arg(short, default_value_t=4.0)]
    d: f64,

    /// Minimum acceleration
    #[arg(long, default_value_t=0.03)]
    mina: f64,

    /// Maximum acceleration
    #[arg(long, default_value_t=3.0)]
    maxa: f64,

    /// Number of accelerations to consider between mina and maxa
    #[arg(long, default_value_t=128)]
    na: usize,

    /// Maximum ammount of points to use in numeric integration
    #[arg(long, default_value_t=1001)]
    n: usize,
}

fn main() {
    let mut args = Args::parse();

    // Ensure n is odd, integrater only works for odd n because it uses simpson's rule
    if args.n % 2 == 0 {
        args.n -= 1;
    }

    let base = gamma0(args.delta, args.m, args.d);
    let result: Vec<_> = (0..args.na)
        .into_par_iter()
        .map(|i| args.mina * ((args.maxa/args.mina).powf((i as f64)/((args.na - 1) as f64))))
        .map(|a| (a, gamma_a(&mut Ctx::new(), a, args.delta, args.m, args.d, args.n) / base))
        .collect();
    println!("{:?}", result);
}
