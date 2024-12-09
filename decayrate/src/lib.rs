use std::f64::consts::PI;

use archt::Ctx;

pub fn gamma0(delta: f64, m: f64, d: f64) -> f64 {
    (delta * delta - m * m).powf((d - 3.0) / 2.0)
}

pub fn gamma_a(ctx: &mut Ctx, a: f64, delta: f64, m: f64, d: f64, n: usize) -> f64 {
    assert_eq!(n % 2, 1);

    let prep = 2.0 * m / (a * PI);
    let r = 2.0 * delta / a;
    let xprop = 2.0 / a;
    let mut summands = vec![0.0; n + 1];
    for i in 1..n {
        let uinv = (n as f64) / (i as f64);
        summands[i] = prep
            * uinv.powi(2)
            * gamma0(uinv, m, d)
            * ctx.kbessel(r, xprop * uinv)
            * (if i % 2 == 1 { 2.0 } else { 1.0 });
    }

    summands.iter().sum::<f64>() / (1.5 * ((n - 1) as f64))
}

pub fn gammma_a_d2(ctx: &mut Ctx, a: f64, delta: f64, m: f64) -> f64 {
    let r = delta / a;
    (1.0 / (a * PI)) * ((ctx.kbessel(r, m / a)).powi(2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let mut ctx = Ctx::new();
        println!("{}", gamma0(2., 1., 4.));
        println!("{}", gamma_a(&mut ctx, 1., 2., 1., 4., 1001));
        println!("{}", gamma_a(&mut ctx, 0.1, 2., 1., 4., 1001));
        println!("{}", gamma_a(&mut ctx, 0.01, 2., 1., 4., 1001));
        println!("{}", gamma_a(&mut ctx, 0.001, 2., 1., 4., 1001));
    }
}
