#![doc = include_str!("../README.md")]

mod slsqp;

pub use slsqp::Func;

use crate::slsqp::{
    nlopt_constraint, nlopt_constraint_raw_callback, nlopt_function_raw_callback, nlopt_stopping,
    NLoptConstraintCfg, NLoptFunctionCfg,
};

use std::os::raw::c_void;

/// Failed termination status of the optimization process
#[derive(Debug, Clone, Copy)]
pub enum FailStatus {
    Failure,
    InvalidArgs,
    OutOfMemory,
    RoundoffLimited,
    ForcedStop,
    UnexpectedError,
}

/// Successful termination status of the optimization process
#[derive(Debug, Clone, Copy)]
pub enum SuccessStatus {
    Success,
    StopValReached,
    FtolReached,
    XtolReached,
    MaxEvalReached,
    MaxTimeReached,
}

/// Outcome when optimization process fails
type FailOutcome = (FailStatus, Vec<f64>, f64);
/// Outcome when optimization process succeeds
type SuccessOutcome = (SuccessStatus, Vec<f64>, f64);

/// Minimizes a function using the SLSQP method.
///
/// # Arguments
///
/// * `func` - the function to minimize
/// * `x0` - the initial guess (will be matated to reflect the argmin result at the end)
/// * `cons` - slice of constraint function intended to be negative at the end
/// * `args` - user data pass to objective and constraint functions
/// * `bounds` - x domain specified as a n-vector of tuple `(lower bound, upper bound)`  
/// * `ftol_rel` - relative tolerance on function value, algorithm stops when `func(x)` changes by less than `ftol_rel * func(x)`
/// * `ftol_abs` - absolute tolerance on function value, algorithm stops when `func(x)` change is less than `ftol_rel`
/// * `xtol_rel` - relative tolerance on optimization parameters, algorithm stops when all `x[i]` changes by less than `xtol_rel * x[i]`
/// * `xtol_abs` - relative tolerance on optimization parameters, algorithm stops when `x[i]` changes by less than `xtol_abs[i]`
/// * `maxeval` - maximum number of objective function evaluation
///     
/// # Returns
///
/// The status of the optimization process, the argmin value and the objective function value
///
/// # Implementation note:
///
/// This implementation is a translation of [NLopt](https://github.com/stevengj/nlopt) 2.7.1
/// See also [NLopt SLSQP](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#slsqp) documentation.
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn minimize<'a, F: Func<U>, G: Func<U>, U: Clone>(
    func: F,
    xinit: &[f64],
    cons: &[G],
    args: U,
    bounds: &[(f64, f64)],
    ftol_rel: f64,
    ftol_abs: f64,
    xtol_rel: f64,
    xtol_abs: &[f64],
    maxiter: usize,
) -> Result<SuccessOutcome, FailOutcome> {
    let fn_cfg = Box::new(NLoptFunctionCfg {
        objective_fn: func,
        user_data: args.clone(), // move user_data into FunctionCfg
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;
    let mut cstr_tol = 2e-4;

    let mut cstr_cfg = cons
        .iter()
        .map(|c| {
            let c_cfg = Box::new(NLoptConstraintCfg {
                constraint_fn: c as &dyn Func<U>,
                user_data: args.clone(), // move user_data into FunctionCfg
            });
            let c_cfg_ptr = Box::into_raw(c_cfg) as *mut c_void;

            nlopt_constraint {
                m: 1,
                f: Some(nlopt_constraint_raw_callback::<F, U>),
                pre: None,
                mf: None,
                f_data: c_cfg_ptr,
                tol: &mut cstr_tol,
            }
        })
        .collect::<Vec<_>>();

    let mut x = vec![0.; xinit.len()];
    x.copy_from_slice(xinit);
    let n = x.len() as u32;
    let m = cons.len() as u32;

    let lbs: Vec<f64> = bounds.iter().map(|b| b.0).collect();
    let ubs: Vec<f64> = bounds.iter().map(|b| b.1).collect();
    let x_weights = vec![0.; n as usize];
    let mut minf = f64::INFINITY;
    let mut nevals_p = 0;
    let mut force_stop = 0;
    let mut stop = nlopt_stopping {
        n,
        minf_max: -f64::INFINITY,
        ftol_rel,
        ftol_abs,
        xtol_rel,
        xtol_abs: xtol_abs.as_ptr(),
        x_weights: x_weights.as_ptr(), // unused
        nevals_p: &mut nevals_p,       // unused
        maxeval: maxiter as i32,
        maxtime: 0.0,                // unused
        start: 0.0,                  // unused
        force_stop: &mut force_stop, // unused
        stop_msg: "".to_string(),    // unused
    };

    let status = unsafe {
        slsqp::nlopt_slsqp::<U>(
            n.into(),
            Some(nlopt_function_raw_callback::<F, U>),
            fn_cfg_ptr,
            m.into(),
            cstr_cfg.as_mut_ptr(),
            0,
            std::ptr::null_mut(),
            lbs.as_ptr(),
            ubs.as_ptr(),
            x.as_mut_ptr(),
            &mut minf,
            &mut stop,
        )
    };

    // Convert the raw pointer back into a Box with the B::from_raw function,
    // allowing the Box destructor to perform the cleanup.
    unsafe {
        let _ = Box::from_raw(fn_cfg_ptr as *mut NLoptFunctionCfg<F, U>);
    };

    match status {
        -1 => Err((FailStatus::Failure, x, minf)),
        -2 => Err((FailStatus::InvalidArgs, x, minf)),
        -3 => Err((FailStatus::OutOfMemory, x, minf)),
        -4 => Err((FailStatus::RoundoffLimited, x, minf)),
        -5 => Err((FailStatus::ForcedStop, x, minf)),
        1 => Ok((SuccessStatus::Success, x, minf)),
        2 => Ok((SuccessStatus::StopValReached, x, minf)),
        3 => Ok((SuccessStatus::FtolReached, x, minf)),
        4 => Ok((SuccessStatus::XtolReached, x, minf)),
        5 => Ok((SuccessStatus::MaxEvalReached, x, minf)),
        6 => Ok((SuccessStatus::MaxTimeReached, x, minf)),
        _ => Err((FailStatus::UnexpectedError, x, minf)),
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    ////////////////////////////////////////////////////////////////////////////////
    /// Nlopt slsqp

    fn raw_paraboloid(
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
    fn test_slsqp_minimize() {
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
            slsqp::nlopt_slsqp::<()>(
                2,
                Some(raw_paraboloid),
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

    fn paraboloid(x: &[f64], gradient: Option<&mut [f64]>, _data: &mut ()) -> f64 {
        println!("{:?}", x);
        let r1 = x[0] + 1.0;
        let r2 = x[1];
        if let Some(g) = gradient {
            g[0] = 20.0 * r1;
            g[1] = 2. * r2;
        }
        10. * r1 * r1 + r2 * r2
    }

    #[test]
    fn test_paraboloid() {
        let xinit = vec![1., 1.];

        let mut cons: Vec<&dyn Func<()>> = vec![];
        let cstr1 = |x: &[f64], gradient: Option<&mut [f64]>, _user_data: &mut ()| {
            if let Some(g) = gradient {
                g[0] = -1.;
                g[1] = 0.;
            }
            -x[0]
        };
        cons.push(&cstr1 as &dyn Func<()>);

        // x_opt = [0, 0]
        match minimize(
            paraboloid,
            &xinit,
            &cons,
            (),
            &[(-10., 10.)],
            1e-4,
            0.0,
            0.0,
            &[0.0, 0.0],
            200,
        ) {
            Ok((_, x, _)) => {
                let exp = [0., 0.];
                for (act, exp) in x.iter().zip(exp.iter()) {
                    assert_abs_diff_eq!(act, exp, epsilon = 1e-3);
                }
            }
            Err((status, _, _)) => {
                panic!("{}", format!("Error status : {:?}", status));
            }
        }
    }
}
