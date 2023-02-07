mod slsqp;
use crate::slsqp::raw_slsqp;

use std::os::raw::c_void;
use std::slice;

// mod cobyla_solver;
// mod cobyla_state;
// pub use crate::cobyla_solver::*;
// pub use crate::cobyla_state::*;

// use std::os::raw::c_void;
// use std::slice;

/// A trait for an objective function to be minimized
///
/// An objective function takes the form of a closure `f(x: &[f64], user_data: &mut U) -> f64`
///
/// * `x` - n-dimensional array
/// * `user_data` - user defined data
pub trait ObjFn<U>: Fn(&[f64], &mut U) -> f64 {}
impl<F, U> ObjFn<U> for F where F: Fn(&[f64], &mut U) -> f64 {}

/// A trait for a constraint function which should be positive eventually
///
/// A constraint function takes the form of a closure `f(x: &[f64]) -> f64`
/// The algorithm makes the constraint positive eventually.
///
/// For instance if you want an upper bound MAX for x,
/// you have to define the constraint as `|x| MAX - x`.
/// Conversly for a lower bound you would define `|x| x - MIN`
///
/// * `x` - n-dimensional array
pub trait CstrFn: Fn(&[f64]) -> f64 {}
impl<F> CstrFn for F where F: Fn(&[f64]) -> f64 {}

/// Packs a function with a user defined parameter set of type `U`
/// and constraints to be made positive eventually by the optimizer
struct FunctionCfg<'a, F: ObjFn<U>, G: CstrFn, U> {
    pub func: F,
    pub cons: &'a [G],
    pub data: U,
}

/// Callback interface for Cobyla C code to evaluate objective and constraint functions
fn function_raw_callback<F: ObjFn<U>, G: CstrFn, U>(
    n: ::std::os::raw::c_long,
    m: ::std::os::raw::c_long,
    x: *const f64,
    con: *mut f64,
    data: *mut ::std::os::raw::c_void,
) -> f64 {
    // prepare args
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    // recover FunctionCfg object from supplied params and call
    let f = unsafe { &mut *(data as *mut FunctionCfg<F, G, U>) };
    let res = (f.func)(argument, &mut f.data);

    for i in 0..m as isize {
        unsafe {
            *con.offset(i) = (f.cons[i as usize])(argument);
        }
    }

    // Important: we don't want f to get dropped at this point
    #[allow(clippy::forget_ref)]
    std::mem::forget(f);
    res
}

/// Minimizes a function using the SLSQP method.
///
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn fmin_slsqp<'a, F: ObjFn<U>, G: CstrFn, U>(
    func: F,
    x0: &'a mut [f64],
    cons: &[G],
    args: U,
    rhobeg: f64,
    rhoend: f64,
    maxfun: i32,
    iprint: i32,
) -> (i32, &'a [f64]) {
    todo!()
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    /////////////////////////////////////////////////////////////////////////
    // First problem (see cobyla.c case 1)

    fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
        10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
    }

    #[test]
    fn test_fmin_slsqp() {
        let mut x = vec![1., 1.];

        let mut cons: Vec<&dyn CstrFn> = vec![];
        let cstr1 = |x: &[f64]| x[0];
        cons.push(&cstr1);

        // x_opt = [0, 0]
        let (status, x_opt) = fmin_slsqp(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 1);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        assert_abs_diff_eq!(x.as_slice(), [0., 0.].as_slice(), epsilon = 1e-4);
    }
}
