use clap::Parser;
use decayrate::gamma0_multiparticle;

/// Calculate relative change in decay rate of a n-particle decay with given parameters for multiple accelerations
#[derive(Parser, Debug, Clone)]
struct Args {
    /// Lower bound delta
    #[arg(long, default_value_t = 0.0)]
    delta_low: f64,

    /// Upper bound delta
    #[arg(long, default_value_t = 1.0)]
    delta_high: f64,

    /// Number of space-time dimensions
    #[arg(short, default_value_t = 4.0)]
    d: f64,

    /// Number of particles
    #[arg(long, default_value_t = 3)]
    np: usize,

    /// Number of sample points
    #[arg(long, default_value_t = 128)]
    nd: usize,

    /// Maximum ammount of points to use in numeric integration
    #[arg(long, default_value_t = 1001)]
    n: usize,
}

fn main() {
    let mut args = Args::parse();

    if args.np == 0 {
        println!("Missing particles");
        return;
    }

    // Ensure n is odd, integrater only works for odd n because it uses simpson's rule
    if args.n % 2 == 0 {
        args.n -= 1;
    }

    let masses = vec![0.0; args.np];

    let result: Vec<_> = (0..args.nd)
        .into_iter()
        .map(|i| {
            args.delta_low
                + (i as f64) * (args.delta_high - args.delta_low) / ((args.nd - 1) as f64)
        })
        .map(|delta| (delta, gamma0_multiparticle(delta, &masses, args.d, args.n)))
        .collect();
    println!("{:?}", result);
}
