mod slsqp;

// use std::os::raw::c_void;
// use std::slice;

// use std::os::raw::c_void;
// use std::slice;

// /// A trait for an objective function to be minimized
// ///
// /// An objective function takes the form of a closure `f(x: &[f64], user_data: &mut U) -> f64`
// ///
// /// * `x` - n-dimensional array
// /// * `user_data` - user defined data
// pub trait ObjFn<U>: Fn(&[f64], &mut U) -> f64 {}
// impl<F, U> ObjFn<U> for F where F: Fn(&[f64], &mut U) -> f64 {}

// /// A trait for a constraint function which should be positive eventually
// ///
// /// A constraint function takes the form of a closure `f(x: &[f64]) -> f64`
// /// The algorithm makes the constraint positive eventually.
// ///
// /// For instance if you want an upper bound MAX for x,
// /// you have to define the constraint as `|x| MAX - x`.
// /// Conversly for a lower bound you would define `|x| x - MIN`
// ///
// /// * `x` - n-dimensional array
// pub trait CstrFn: Fn(&[f64]) -> f64 {}
// impl<F> CstrFn for F where F: Fn(&[f64]) -> f64 {}

// /// Packs a function with a user defined parameter set of type `U`
// /// and constraints to be made positive eventually by the optimizer
// struct FunctionCfg<'a, F: ObjFn<U>, G: CstrFn, U> {
//     pub func: F,
//     pub cons: &'a [G],
//     pub data: U,
// }

// /// Callback interface for Cobyla C code to evaluate objective and constraint functions
// fn function_raw_callback<F: ObjFn<U>, G: CstrFn, U>(
//     n: ::std::os::raw::c_long,
//     m: ::std::os::raw::c_long,
//     x: *const f64,
//     con: *mut f64,
//     data: *mut ::std::os::raw::c_void,
// ) -> f64 {
//     // prepare args
//     let argument = unsafe { slice::from_raw_parts(x, n as usize) };
//     // recover FunctionCfg object from supplied params and call
//     let f = unsafe { &mut *(data as *mut FunctionCfg<F, G, U>) };
//     let res = (f.func)(argument, &mut f.data);

//     for i in 0..m as isize {
//         unsafe {
//             *con.offset(i) = (f.cons[i as usize])(argument);
//         }
//     }

//     // Important: we don't want f to get dropped at this point
//     #[allow(clippy::forget_ref)]
//     std::mem::forget(f);
//     res
// }

// /// Minimizes a function using the SLSQP method.
// ///
// #[allow(clippy::useless_conversion)]
// #[allow(clippy::too_many_arguments)]
// pub fn fmin_slsqp<'a, F: ObjFn<U>, G: CstrFn, U>(
//     func: F,
//     x0: &'a mut [f64],
//     cons: &[G],
//     args: U,
//     rhobeg: f64,
//     rhoend: f64,
//     maxfun: i32,
//     iprint: i32,
// ) -> (i32, &'a [f64]) {
//     todo!()
// }

#[cfg(test)]
mod tests {
    use crate::slsqp::{nlopt_slsqp, nlopt_stopping};
    use approx::assert_abs_diff_eq;

    // use super::*;

    ////////////////////////////////////////////////////////////////////////////////
    /// Nlopt slsqp

    fn nlopt_raw_paraboloid(
        _n: libc::c_uint,
        x: *const libc::c_double,
        gradient: *mut libc::c_double,
        _func_data: *mut libc::c_void,
    ) -> libc::c_double {
        unsafe {
            let r1 = *x.offset(0) + 1.0;
            let r2 = *x.offset(1);

            if !gradient.is_null() {
                *gradient.offset(0) = 20.0 * r1;
                *gradient.offset(1) = 2. * r2;
            }

            10.0 * (r1 * r1) + (r2 * r2) as libc::c_double
        }
    }

    #[test]
    fn test_nlopt_slsqp_minimize() {
        let mut x = vec![1., 1.];
        let mut lb = vec![-10.0, -10.0];
        let mut ub = vec![10.0, 10.0];
        let mut minf = f64::INFINITY;
        let mut nevals_p = 0;
        let mut force_stop = 0;

        let mut stop = nlopt_stopping {
            n: 2,
            minf_max: -f64::INFINITY,
            ftol_rel: 0.0,
            ftol_abs: 0.0,
            xtol_rel: 0.0,
            xtol_abs: std::ptr::null(),
            x_weights: std::ptr::null(),
            nevals_p: &mut nevals_p,
            maxeval: 1000,
            maxtime: 0.0,
            start: 0.0,
            force_stop: &mut force_stop,
            stop_msg: "".to_string(),
        };

        let res = unsafe {
            nlopt_slsqp(
                2,
                Some(nlopt_raw_paraboloid),
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                lb.as_mut_ptr(),
                ub.as_mut_ptr(),
                x.as_mut_ptr(),
                &mut minf,
                &mut stop,
            )
        };

        println!("status = {:?}", res);
        println!("x = {:?}", x);

        assert_abs_diff_eq!(x.as_slice(), [-1., 0.].as_slice(), epsilon = 1e-4);
    }
}
