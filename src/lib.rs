mod slsqp;

use crate::slsqp::{
    nlopt_constraint, nlopt_constraint_raw_callback, nlopt_function_raw_callback, nlopt_stopping,
    NLoptConstraintCfg, NLoptFunctionCfg, NLoptObjFn,
};
use std::os::raw::c_void;

/// Minimizes a function using the SLSQP method.
/// This implementation is a translation of [NLopt](https://github.com/stevengj/nlopt) 2.7.1
///
/// See [NLopt SLSQP](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#slsqp) documentation.
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn nlopt_slsqp<'a, F: NLoptObjFn<U>, G: NLoptObjFn<U>, U: Clone>(
    func: F,
    x0: &'a mut [f64],
    cons: &[G],
    args: U,
    rhobeg: f64,
    rhoend: f64,
    maxfun: i32,
    _iprint: i32,
    bounds: (f64, f64),
) -> (i32, &'a [f64]) {
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
                constraint_fn: c as &dyn NLoptObjFn<U>,
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

    let x = x0;
    let n = x.len() as u32;
    let m = cons.len() as u32;

    let lb = vec![bounds.0; n as usize];
    let ub = vec![bounds.1; n as usize];
    let xtol_abs = vec![0.; n as usize];
    let x_weights = vec![0.; n as usize];
    let mut minf = f64::INFINITY;
    let mut nevals_p = 0;
    let mut force_stop = 0;
    let mut stop = nlopt_stopping {
        n,
        minf_max: -f64::INFINITY,
        ftol_rel: 1e-4,
        ftol_abs: 0.0,
        xtol_rel: rhoend / rhobeg,
        xtol_abs: xtol_abs.as_ptr(),
        x_weights: x_weights.as_ptr(),
        nevals_p: &mut nevals_p,
        maxeval: maxfun.into(),
        maxtime: 0.0,
        start: 0.0,
        force_stop: &mut force_stop,
        stop_msg: "".to_string(),
    };

    // XXX: Weird bug. Can not pass nlopt_constraint_raw_callback
    // Work around is to patch nlopt_eval_constraint to use nlopt_constraint_raw_callback directly
    // if !cons.is_empty() {
    //     let mut xtest = vec![1., 1.];
    //     unsafe {
    //         let mut result = -666.;
    //         let nc = cstr_cfg[0];
    //         let fc = &nc as *const nlopt_constraint;

    //         // It works: cstr1 is called
    //         let _res = nlopt_constraint_raw_callback::<&dyn NLoptObjFn<()>, ()>(
    //             2,
    //             xtest.as_mut_ptr(),
    //             std::ptr::null_mut::<libc::c_double>(),
    //             (nc).f_data,
    //         );
    //         // println!(
    //         //     "###################################### JUST direct nlopt_constraint_raw_callback is OK = {}",
    //         //     res,
    //         // );

    //         // XXX: Weird bug!
    //         // If fails : (*fc).f does call nlopt_constraint_raw_callback but
    //         // when unpacking (*fc).f_data with unsafe f = { &mut *(params as *mut NLoptConstraintCfg<F, T>) };
    //         // (*f).constraint_fn(...) calls the objective function instead of the constraint function!!!
    //         let _res = ((*fc).f.expect("func"))(
    //             2,
    //             xtest.as_mut_ptr(),
    //             std::ptr::null_mut::<libc::c_double>(),
    //             (*fc).f_data,
    //         );
    //         // println!(
    //         //     "###################################### JUST stored nlopt_constraint_raw_callback is NOT OK = {}",
    //         //     res,
    //         // );

    //         // It works: cstr1 is called
    //         // we use directly a copy of specialized nlopt_constraint_raw_callback
    //         nlopt_eval_constraint::<()>(
    //             &mut result,
    //             std::ptr::null_mut::<libc::c_double>(),
    //             fc,
    //             n,
    //             x.as_mut_ptr(),
    //         );
    //         // println!(
    //         //     "############################### TEST nlopt_eval_constraint (OK if opposite previous OK result) = {}",
    //         //     result
    //         // );
    //     }
    // }

    let status = unsafe {
        slsqp::nlopt_slsqp::<U>(
            n.into(),
            Some(nlopt_function_raw_callback::<F, U>),
            fn_cfg_ptr,
            m.into(),
            cstr_cfg.as_mut_ptr(),
            0,
            std::ptr::null_mut(),
            lb.as_ptr(),
            ub.as_ptr(),
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
    (status, x)
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

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
            slsqp::nlopt_slsqp::<()>(
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

    fn nlopt_paraboloid(x: &[f64], gradient: Option<&mut [f64]>, _data: &mut ()) -> f64 {
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
    fn test_nlopt_paraboloid() {
        let mut x = vec![1., 1.];

        let mut cons: Vec<&dyn NLoptObjFn<()>> = vec![];
        let cstr1 = |x: &[f64], gradient: Option<&mut [f64]>, _user_data: &mut ()| {
            if let Some(g) = gradient {
                g[0] = -1.;
                g[1] = 0.;
            }
            -x[0]
        };
        cons.push(&cstr1 as &dyn NLoptObjFn<()>);

        // x_opt = [0, 0]
        let (status, x_opt) = nlopt_slsqp(
            nlopt_paraboloid,
            &mut x,
            &cons,
            (),
            0.5,
            0.0,
            200,
            1,
            (-10., 10.),
        );
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        assert_abs_diff_eq!(x.as_slice(), [0., 0.].as_slice(), epsilon = 1e-3);
    }
}
