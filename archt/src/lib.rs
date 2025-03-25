use mpfi_sys::mpfi_ptr;
use std::{ffi::c_ulong, ptr::null_mut};

extern "C" {
    //void clear_kbessel(mpfi_ptr f)
    fn clear_kbessel(f: mpfi_ptr);
    //mpfi_ptr init_kbessel(unsigned long d)
    fn init_kbessel(d: c_ulong) -> mpfi_ptr;
    //double kbessel_wrap(mpfi_ptr f, double r, double x)
    fn kbessel_wrap(f: mpfi_ptr, r: f64, x: f64) -> f64;
    fn kbessel_both(f: mpfi_ptr, r: f64, x: f64, left: *mut f64, right: *mut f64);
}

#[derive(Debug)]
pub struct Ctx {
    precision: u32,
    ctx: mpfi_ptr,
}

impl Ctx {
    pub fn new(precision: u32) -> Self {
        let ctx = unsafe { init_kbessel(precision as _) };
        assert!(ctx != null_mut());
        Ctx { precision, ctx }
    }

    pub fn kbessel(&mut self, r: f64, x: f64) -> f64 {
        unsafe { kbessel_wrap(self.ctx, r, x) }
    }

    pub fn kbessel_bound(&mut self, r: f64, x: f64) -> (f64, f64) {
        let mut left: f64 = 0.0;
        let mut right: f64 = 0.0;
        unsafe { kbessel_both(self.ctx, r, x, &mut left, &mut right) };
        (left, right)
    }
}

impl Clone for Ctx {
    fn clone(&self) -> Self {
        Self::new(self.precision)
    }
}

unsafe impl Send for Ctx {}

impl Drop for Ctx {
    fn drop(&mut self) {
        unsafe { clear_kbessel(self.ctx) };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanity() {
        let mut ctx = Ctx::new(53);
        let a = ctx.kbessel(1.0, 1.0);
        assert!(a - 1.39 < 0.005);
    }

    #[test]
    fn repeated() {
        let mut ctx = Ctx::new(53);
        for i in 1..100 {
            for j in 1..100 {
                let _ = ctx.kbessel(i as f64 * 1e-2, j as f64 * 1e-2);
            }
        }
    }
}
