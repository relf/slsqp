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

/// Tolerances used as termination criteria.
/// For all, condition is disabled if value is not strictly positive.
/// ```rust
/// # use crate::slsqp::StopTols;
/// let stop_tol = StopTols {
///     ftol_rel: 1e-4,
///     xtol_abs: vec![1e-3; 3],   // size should be equal to x dim
///     ..StopTols::default()      // default stop conditions are disabled
/// };  
/// ```
#[derive(Debug, Clone, Default)]
pub struct StopTols {
    /// Relative tolerance on function value, algorithm stops when `func(x)` changes by less than `ftol_rel * func(x)`
    pub ftol_rel: f64,
    /// Absolute tolerance on function value, algorithm stops when `func(x)` change is less than `ftol_rel`
    pub ftol_abs: f64,
    /// Relative tolerance on optimization parameters, algorithm stops when all `x[i]` changes by less than `xtol_rel * x[i]`
    pub xtol_rel: f64,
    /// Relative tolerance on optimization parameters, algorithm stops when `x[i]` changes by less than `xtol_abs[i]`
    pub xtol_abs: Vec<f64>,
}

/// Minimizes a function using the SLSQP method.
///
/// # Arguments
///
/// * `func` - the function to minimize
/// * `xinit` - n-vector the initial guess
/// * `bounds` - x domain specified as a n-vector of tuple `(lower bound, upper bound)`  
/// * `cons` - slice of constraint function intended to be negative at the end
/// * `args` - user data pass to objective and constraint functions
/// * `maxeval` - maximum number of objective function evaluation
///     
/// ## Returns
///
/// The status of the optimization process, the argmin value and the objective function value
///
/// ## Panics
///
/// When some vector arguments like `bounds`, `xtol_abs` do not have the same size as `xinit`
///
/// ## Example
/// ```
/// # use approx::assert_abs_diff_eq;
/// use slsqp::{minimize, Func};
///
/// fn paraboloid(x: &[f64], gradient: Option<&mut [f64]>, _data: &mut ()) -> f64 {
///     let r1 = x[0] + 1.0;
///     let r2 = x[1];
///     if let Some(g) = gradient {
///         g[0] = 20.0 * r1;
///         g[1] = 2. * r2;
///     }
///     10. * r1 * r1 + r2 * r2
/// }
///
/// // Initial guess
/// let mut x = vec![1., 1.];
///
/// // Constraints definition to be negative eventually: here `x_0 > 0`
/// // hence the `-x[0]` value returned below
/// let cstr1 = |x: &[f64], gradient: Option<&mut [f64]>, _user_data: &mut ()| {
///     if let Some(g) = gradient {
///         g[0] = -1.;
///         g[1] = 0.;
///     }
///     -x[0]
/// };
/// let cons: Vec<&dyn Func<()>> = vec![&cstr1];
///
/// match minimize(
///     paraboloid,
///     &mut x,
///     &[(-10., 10.), (-10., 10.)],
///     &cons,
///     (),
///     100,
///     None
/// ) {
///     Ok((status, x_opt, y_opt)) => {
///         println!("status = {:?}", status);
///         println!("x_opt = {:?}", x_opt);
///         println!("y_opt = {}", y_opt);
/// #       assert_abs_diff_eq!(y_opt, 10.0, epsilon=1e-8);
///     }
///     Err((e, _, _)) => panic!("Optim error: {:?}", e),
/// }
/// ```
///
/// ## Implementation note:
///
/// This implementation is a translation of [NLopt](https://github.com/stevengj/nlopt) 2.7.1
/// See also [NLopt SLSQP](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#slsqp) documentation.
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn minimize<F: Func<U>, G: Func<U>, U: Clone>(
    func: F,
    xinit: &[f64],
    bounds: &[(f64, f64)],
    cons: &[G],
    args: U,
    maxeval: usize,
    stop_tol: Option<StopTols>,
) -> Result<SuccessOutcome, FailOutcome> {
    let fn_cfg = Box::new(NLoptFunctionCfg {
        objective_fn: func,
        user_data: args.clone(),
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;
    let mut cstr_tol = 0.0; // no cstr tolerance

    let mut cstr_cfg = cons
        .iter()
        .map(|c| {
            let c_cfg = Box::new(NLoptConstraintCfg {
                constraint_fn: c as &dyn Func<U>,
                user_data: args.clone(),
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

    if bounds.len() != x.len() {
        panic!(
            "{}",
            format!(
                "Minimize Error: Bounds and x should have same size! Got {} for bounds and {} for x.",
                bounds.len(),
                x.len()
            )
        )
    }
    let lbs: Vec<f64> = bounds.iter().map(|b| b.0).collect();
    let ubs: Vec<f64> = bounds.iter().map(|b| b.1).collect();

    let x_weights = vec![0.; n as usize];
    let mut minf = f64::INFINITY;
    let mut nevals_p = 0;
    let mut force_stop = 0;

    let stop_tol = stop_tol.unwrap_or_default();
    let xtol_abs = if stop_tol.xtol_abs.is_empty() {
        std::ptr::null()
    } else if stop_tol.xtol_abs.len() != n as usize {
        panic!(
            "{}",
            format!(
                "Minimize Error: xtol_abs should have x dim size ({}), got {}",
                n,
                stop_tol.xtol_abs.len()
            )
        );
    } else {
        stop_tol.xtol_abs.as_ptr()
    };
    let mut stop = nlopt_stopping {
        n,
        minf_max: -f64::INFINITY,
        ftol_rel: stop_tol.ftol_rel,
        ftol_abs: stop_tol.ftol_abs,
        xtol_rel: stop_tol.xtol_rel,
        xtol_abs,
        x_weights: x_weights.as_ptr(), // unused
        nevals_p: &mut nevals_p,       // unused
        maxeval: maxeval as i32,
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

        let stop_tol = StopTols {
            ftol_rel: 1e-4,
            ..StopTols::default()
        };

        // x_opt = [0, 0]
        match minimize(
            paraboloid,
            &xinit,
            &[(-10., 10.), (-10., 10.)],
            &cons,
            (),
            200,
            Some(stop_tol),
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

    fn xsinx(x: &[f64], gradient: Option<&mut [f64]>, _user_data: &mut ()) -> f64 {
        let r = (x[0] - 3.5) / std::f64::consts::PI;
        if let Some(g) = gradient {
            g[0] = f64::sin(r) + (x[0] - 3.5) * f64::cos(r) / std::f64::consts::PI;
        }
        (x[0] - 3.5) * f64::sin(r)
    }

    #[test]
    fn test_xsinx() {
        let xinit = vec![10.];

        let cons: Vec<&dyn Func<()>> = vec![];

        // x_opt = [18.9352]
        match minimize(xsinx, &xinit, &[(0., 25.)], &cons, (), 200, None) {
            Ok((_, x, _)) => {
                let exp = [18.9352];
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
