use archt::Ctx;
use clap::{Parser, ValueEnum};
use decayrate::{gamma0, gamma_a};

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, ValueEnum)]
enum ExtremumType {
    Maximum,
    Minimum,
}

/// Calculate relative change in decay rate of a 2-particle decay with given parameters for multiple accelerations
#[derive(Parser, Debug, Copy, Clone)]
struct Args {
    /// Mass difference between decaying and heaviest particle
    #[arg(long, default_value_t = 1.5)]
    delta: f64,

    /// Mass of lightest particle
    #[arg(short, default_value_t = 1.0)]
    m: f64,

    /// Number of space-time dimensions
    #[arg(short, default_value_t = 4.0)]
    d: f64,

    /// Lower bound of the search range
    #[arg(long, default_value_t = 0.15)]
    mina: f64,

    /// Upper bound of the search range
    #[arg(long, default_value_t = 0.3)]
    maxa: f64,

    /// Relative precision of the result
    #[arg(long, default_value_t = 1e-3)]
    precision: f64,

    /// Maximum ammount of points to use in numeric integration
    #[arg(long, default_value_t = 1001)]
    n: usize,

    /// Type of extremum to look for
    #[arg(long, default_value = "maximum")]
    r#type: ExtremumType,

    /// Precision of the bessel function calculations
    #[arg(long, default_value_t = 53)]
    bessel_precision: u32,
}

fn main() {
    let mut args = Args::parse();

    // Ensure n is odd, integrater only works for odd n because it uses simpson's rule
    if args.n % 2 == 0 {
        args.n -= 1;
    }

    let mut low = args.mina;
    let mut high = args.maxa;

    let mut ctx = Ctx::new(args.bessel_precision);

    while high - low > low * args.precision {
        let step = (high - low) / 3.;
        let a = low + step;
        let b = high - step;

        let fa = gamma_a(&mut ctx, a, args.delta, args.m, args.n, |dp| {
            gamma0(dp, args.m, args.d)
        });
        let fb = gamma_a(&mut ctx, b, args.delta, args.m, args.n, |dp| {
            gamma0(dp, args.m, args.d)
        });

        match (args.r#type, fa < fb) {
            (ExtremumType::Maximum, true) => low = a,
            (ExtremumType::Maximum, false) => high = b,
            (ExtremumType::Minimum, true) => high = b,
            (ExtremumType::Minimum, false) => low = a,
        }
    }

    let result = (high + low) / 2.0;
    let result_delta = gamma_a(&mut ctx, result, args.delta, args.m, args.n, |dp| {
        gamma0(dp, args.m, args.d)
    }) / gamma0(args.delta, args.m, args.d);

    println!(
        "Extrema at {}, relative decay rate {}",
        result, result_delta
    );
}
