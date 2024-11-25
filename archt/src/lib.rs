use std::{ffi::c_ulong, ptr::null_mut};
use mpfi_sys::mpfi_ptr;

extern "C" {
    //void clear_kbessel(mpfi_ptr f)
    fn clear_kbessel(f: mpfi_ptr);
    //mpfi_ptr init_kbessel(unsigned long d)
    fn init_kbessel(d: c_ulong) -> mpfi_ptr;
    //double kbessel_wrap(mpfi_ptr f, double r, double x)
    fn kbessel_wrap(f: mpfi_ptr, r: f64, x: f64) -> f64;
}

#[derive(Debug)]
pub struct Ctx {
    ctx: mpfi_ptr,
}

impl Ctx {
    pub fn new() -> Self {
        let ctx = unsafe { init_kbessel(53) };
        assert!(ctx != null_mut());
        Ctx {
            ctx,
        }
    }

    pub fn kbessel(&mut self, r: f64, x: f64) -> f64 {
        unsafe { kbessel_wrap(self.ctx, r, x) }
    }
}

impl Drop for Ctx {
    fn drop(&mut self) {
        unsafe {clear_kbessel(self.ctx)};
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanity() {
        let mut ctx = Ctx::new();
        let a = ctx.kbessel(1.0, 1.0);
        assert!(a - 1.39 < 0.005);
    }

    #[test]
    fn repeated() {
        let mut ctx = Ctx::new();
        for i in 1..100 {
            for j in 1..100 {
                let _ = ctx.kbessel(i as f64 * 1e-2, j as f64 * 1e-2);
            }
        }
    }
}
