use archt::Ctx;
use clap::Parser;
use decayrate::{gamma0_weakdecay, gamma_a};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Calculate relative change in decay rate of a n-particle decay with given parameters for multiple accelerations
#[derive(Parser, Debug, Clone)]
struct Args {
    /// Mass difference between decaying and heaviest particle
    #[arg(long, default_value_t = 1.5)]
    delta: f64,

    /// Mass of the electron
    #[arg(short, default_value_t = 1.0)]
    m: f64,

    /// Minimum acceleration
    #[arg(long, default_value_t = 0.03)]
    mina: f64,

    /// Maximum acceleration
    #[arg(long, default_value_t = 3.0)]
    maxa: f64,

    /// Number of accelerations to consider between mina and maxa
    #[arg(long, default_value_t = 128)]
    na: usize,

    /// Maximum ammount of points to use in numeric integration
    #[arg(long, default_value_t = 1001)]
    n: usize,

    /// Maximum ammount of points to use in numeric integration of inner decay rate
    #[arg(long, default_value_t = 101)]
    ni: usize,

    /// Precision of the bessel function calculations
    #[arg(long, default_value_t = 53)]
    precision: u32,
}

fn main() {
    let mut args = Args::parse();

    // Ensure n is odd, integrater only works for odd n because it uses simpson's rule
    if args.n % 2 == 0 {
        args.n -= 1;
    }

    let base = gamma0_weakdecay(args.delta, args.m, args.ni);
    let result: Vec<_> = (0..args.na)
        .into_par_iter()
        .map(|i| args.mina * ((args.maxa / args.mina).powf((i as f64) / ((args.na - 1) as f64))))
        .map(|a| {
            (
                a,
                gamma_a(
                    &mut Ctx::new(args.precision),
                    a,
                    args.delta,
                    args.m,
                    args.n,
                    |dp| gamma0_weakdecay(dp, args.m, args.ni),
                ) / base,
            )
        })
        .collect();
    println!("{:?}", result);
}
