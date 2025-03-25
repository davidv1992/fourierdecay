use archt::Ctx;
use clap::Parser;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Calculate relative change in decay rate of a n-particle decay with given parameters for multiple accelerations
#[derive(Parser, Debug, Clone)]
struct Args {
    /// Mass difference between decaying and heaviest particle
    #[arg(long, default_value_t = 1.5)]
    delta: f64,

    /// Acceleration
    #[arg(short, default_value_t = 1.0)]
    a: f64,

    /// Number of points to consider for u
    #[arg(long, default_value_t = 128)]
    nu: usize,

    /// Precision of the bessel function calculations
    precision: u32,
}

fn main() {
    let args = Args::parse();

    let r = 2.0 * args.delta / args.a;

    let result: Vec<_> = (0..args.nu)
        .into_par_iter()
        .map(|i| (i as f64) / ((args.nu - 1) as f64))
        .map_with(Ctx::new(args.precision), |ctx, u| {
            let (l, r) = ctx.kbessel_bound(r, 2.0 / (args.a * u));
            (u, l, r)
        })
        .collect();
    println!("{:?}", result);
}
