use std::f64::consts::PI;

use archt::Ctx;

pub fn gamma0(delta: f64, m: f64, d: f64) -> f64 {
    (delta * delta - m * m).powf((d - 3.0) / 2.0)
}

pub fn gamma0_multiparticle(delta: f64, m: &[f64], d: f64, n: usize) -> f64 {
    if m.iter().sum::<f64>() >= delta {
        return 0.0;
    }

    match m.split_first() {
        Some((m, [])) => (delta * delta - m * m).powf((d - 3.0) / 2.0),
        Some((m, rest)) => {
            let mut summands = vec![0.0; n];
            let sum_rest = rest.iter().sum::<f64>();
            let limit = ((delta - sum_rest - m) * (delta - sum_rest + m)).sqrt();
            let step = limit / ((n - 1) as f64);

            assert_eq!(n % 2, 1);

            for i in 1..n {
                let p = i as f64 * step;
                let e = (m.powi(2) + p.powi(2)).sqrt();
                summands[i] = (p.powf(d - 2.0) * gamma0_multiparticle(delta - e, rest, d, n) / e)
                    * (if i % 2 == 1 { 2.0 } else { 1.0 });
            }

            summands[0] /= 2.0;
            summands[n - 1] /= 2.0;

            summands.iter().sum::<f64>() / (1.5 * ((n - 1) as f64)) * limit
        }
        None => 0.0,
    }
}

pub fn gamma0_weakdecay(delta: f64, m: f64, n: usize) -> f64 {
    assert_eq!(n % 2, 1);
    let limit = (delta * delta - m * m) / (2.0 * delta);
    let step = limit / ((n - 1) as f64);
    let mut summands = vec![0.0; n];

    for i in 0..n {
        let p = i as f64 * step;
        summands[i] =
            p * p * (delta * delta + m * m + p * p - 2.0 * delta * ((m * m + p * p).sqrt()));
    }

    summands[n - 1] /= 2.0;
    summands[0] /= 2.0;

    summands.iter().sum::<f64>() / (1.5 * ((n - 1) as f64)) * limit
}

#[cfg(test)]
fn gamma_a_old(ctx: &mut Ctx, a: f64, delta: f64, m: f64, d: f64, n: usize) -> f64 {
    assert_eq!(n % 2, 1);

    let prep = 2.0 * m / (a * PI);
    let r = 2.0 * delta / a;
    let xprop = 2.0 * m / a;
    let mut summands = vec![0.0; n + 1];
    for i in 1..n {
        let uinv = (n as f64) / (i as f64);
        summands[i] = prep
            * uinv.powi(2)
            * gamma0(m * uinv, m, d)
            * ctx.kbessel(r, xprop * uinv)
            * (if i % 2 == 1 { 2.0 } else { 1.0 });
    }

    summands.iter().sum::<f64>() / (1.5 * ((n - 1) as f64))
}

pub fn gamma_a(
    ctx: &mut Ctx,
    a: f64,
    delta: f64,
    e_min: f64,
    n: usize,
    gamma_0: impl Fn(f64) -> f64,
) -> f64 {
    assert_eq!(n % 2, 1);

    let prep = 2.0 * e_min / (a * PI);
    let r = 2.0 * delta / a;
    let xprop = 2.0 * e_min / a;
    let mut summands = vec![0.0; n + 1];
    for i in 1..n {
        let uinv = (n as f64) / (i as f64);
        summands[i] = prep
            * uinv.powi(2)
            * gamma_0(e_min * uinv)
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
        let mut ctx = Ctx::new(53);

        assert_eq!(
            gamma_a(&mut ctx, 1., 2., 1., 1001, |dp| gamma0(dp, 1., 4.)),
            gamma_a_old(&mut ctx, 1., 2., 1., 4., 1001)
        );
        assert_eq!(
            gamma_a(&mut ctx, 0.1, 2., 1., 1001, |dp| gamma0(dp, 1., 4.)),
            gamma_a_old(&mut ctx, 0.1, 2., 1., 4., 1001)
        );
    }

    #[test]
    fn multiparticle() {
        let a = gamma0_multiparticle(10.0, &[1.0, 2.0], 4.0, 1001);
        let b = gamma0_multiparticle(10.0, &[2.0, 1.0], 4.0, 1001);
        assert!((a - b).abs() < 1e-3);

        assert_eq!(
            gamma0(10.0, 1.0, 4.0),
            gamma0_multiparticle(10.0, &[1.0], 4.0, 101)
        );

        let mut ctx = Ctx::new(53);
        let c = gamma_a(&mut ctx, 1.5, 1.5, 1.0, 1001, |dp| {
            gamma0_multiparticle(dp, &[0.4, 0.6], 4., 101)
        });
        let d = gamma_a(&mut ctx, 1.5, 1.5, 1.0, 1001, |dp| {
            gamma0_multiparticle(dp, &[0.6, 0.4], 4., 101)
        });
        assert!((c - d).abs() < 1e-5);

        println!("{}", gamma0_multiparticle(1.5, &[1.0, 0.0, 0.0], 4.0, 101));
    }
}
