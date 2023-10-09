#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut,
    clippy::needless_return,
    clippy::zero_ptr,
    clippy::toplevel_ref_arg,
    clippy::nonminimal_bool,
    clippy::assign_op_pattern,
    clippy::collapsible_if,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::single_match,
    clippy::unnecessary_cast,
    clippy::excessive_precision,
    clippy::too_many_arguments
)]

use std::slice;
use std::time::{SystemTime, UNIX_EPOCH};

pub fn nlopt_function_raw_callback<F: NLoptObjFn<T>, T>(
    n: libc::c_uint,
    x: *const f64,
    g: *mut f64,
    params: *mut libc::c_void,
) -> f64 {
    // prepare args
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    let gradient = if g.is_null() {
        None
    } else {
        Some(unsafe { slice::from_raw_parts_mut(g, n as usize) })
    };

    // recover FunctionCfg object from supplied params and call
    let f = unsafe { &mut *(params as *mut NLoptFunctionCfg<F, T>) };
    let res = (f.objective_fn)(argument, gradient, &mut f.user_data);
    #[allow(forgetting_references)]
    std::mem::forget(f);
    res
}

pub fn nlopt_constraint_raw_callback<F: NLoptObjFn<T>, T>(
    n: libc::c_uint,
    x: *const f64,
    g: *mut f64,
    params: *mut libc::c_void,
) -> f64 {
    let f = unsafe { &mut *(params as *mut NLoptConstraintCfg<F, T>) };
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    let gradient = if g.is_null() {
        None
    } else {
        Some(unsafe { slice::from_raw_parts_mut(g, n as usize) })
    };
    (f.constraint_fn)(argument, gradient, &mut f.user_data)
}

/// Packs an objective function with a user defined parameter set of type `T`.
pub struct NLoptFunctionCfg<F: NLoptObjFn<T>, T> {
    pub objective_fn: F,
    pub user_data: T,
}

pub struct NLoptConstraintCfg<F: NLoptObjFn<T>, T> {
    pub constraint_fn: F,
    pub user_data: T,
}

/// A trait representing an objective function.
///
/// An objective function takes the form of a closure `f(x: &[f64], gradient: Option<&mut [f64], user_data: &mut U) -> f64`
///
/// * `x` - `n`-dimensional array
/// * `gradient` - `n`-dimensional array to store the gradient `grad f(x)`. If `gradient` matches
/// `Some(x)`, the user is required to provide a gradient, otherwise the optimization will
/// probabely fail.
/// * `user_data` - user defined data
pub trait NLoptObjFn<U>: Fn(&[f64], Option<&mut [f64]>, &mut U) -> f64 {}
impl<T, U> NLoptObjFn<U> for T where T: Fn(&[f64], Option<&mut [f64]>, &mut U) -> f64 {}

// #![register_tool(c2rust)]
// #![feature(c_variadic, register_tool)]
extern "C" {
    // fn sqrt(_: libc::c_double) -> libc::c_double;
    // fn fabs(_: libc::c_double) -> libc::c_double;
    // fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    // fn realloc(_: *mut libc::c_void, _: libc::c_ulong) -> *mut libc::c_void;
    // fn free(__ptr: *mut libc::c_void);
    // fn abort() -> !;
    // fn memcpy(_: *mut libc::c_void, _: *const libc::c_void, _: libc::c_ulong) -> *mut libc::c_void;
    // fn strlen(_: *const libc::c_char) -> libc::c_ulong;
    // fn gettimeofday(__tv: *mut timeval, __tz: *mut libc::c_void) -> libc::c_int;
    // fn vsnprintf(
    //     _: *mut libc::c_char,
    //     _: libc::c_ulong,
    //     _: *const libc::c_char,
    //     _: ::std::ffi::VaList,
    // ) -> libc::c_int;
}

unsafe fn memcpy(dst: *mut libc::c_void, src: *const libc::c_void, n: libc::c_ulong) {
    std::ptr::copy_nonoverlapping(src, dst, n as usize);
}

// pub type __builtin_va_list = [__va_list_tag; 1];
// #[derive(Copy, Clone)]
// #[repr(C)]
// pub struct __va_list_tag {
//     pub gp_offset: libc::c_uint,
//     pub fp_offset: libc::c_uint,
//     pub overflow_arg_area: *mut libc::c_void,
//     pub reg_save_area: *mut libc::c_void,
// }
pub type __time_t = libc::c_long;
pub type __suseconds_t = libc::c_long;
pub type size_t = libc::c_ulong;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct timeval {
    pub tv_sec: __time_t,
    pub tv_usec: __suseconds_t,
}
pub type nlopt_func = Option<
    unsafe fn(
        libc::c_uint,
        *const libc::c_double,
        *mut libc::c_double,
        *mut libc::c_void,
    ) -> libc::c_double,
>;
pub type nlopt_mfunc = Option<
    unsafe fn(
        libc::c_uint,
        *mut libc::c_double,
        libc::c_uint,
        *const libc::c_double,
        *mut libc::c_double,
        *mut libc::c_void,
    ) -> (),
>;
pub type nlopt_precond = Option<
    unsafe fn(
        libc::c_uint,
        *const libc::c_double,
        *const libc::c_double,
        *mut libc::c_double,
        *mut libc::c_void,
    ) -> (),
>;
pub type nlopt_result = libc::c_int;
pub const NLOPT_NUM_RESULTS: nlopt_result = 7;
pub const NLOPT_MAXTIME_REACHED: nlopt_result = 6;
pub const NLOPT_MAXEVAL_REACHED: nlopt_result = 5;
pub const NLOPT_XTOL_REACHED: nlopt_result = 4;
pub const NLOPT_FTOL_REACHED: nlopt_result = 3;
pub const NLOPT_STOPVAL_REACHED: nlopt_result = 2;
pub const NLOPT_SUCCESS: nlopt_result = 1;
pub const NLOPT_NUM_FAILURES: nlopt_result = -6;
pub const NLOPT_FORCED_STOP: nlopt_result = -5;
pub const NLOPT_ROUNDOFF_LIMITED: nlopt_result = -4;
pub const NLOPT_OUT_OF_MEMORY: nlopt_result = -3;
pub const NLOPT_INVALID_ARGS: nlopt_result = -2;
pub const NLOPT_FAILURE: nlopt_result = -1;
// pub type va_list = __builtin_va_list;
// #[derive(Copy, Clone)]
#[derive(Clone)]
#[repr(C)]
pub struct nlopt_stopping {
    pub n: libc::c_uint,
    pub minf_max: libc::c_double,
    pub ftol_rel: libc::c_double,
    pub ftol_abs: libc::c_double,
    pub xtol_rel: libc::c_double,
    pub xtol_abs: *const libc::c_double,
    pub x_weights: *const libc::c_double,
    pub nevals_p: *mut libc::c_int,
    pub maxeval: libc::c_int,
    pub maxtime: libc::c_double,
    pub start: libc::c_double,
    pub force_stop: *mut libc::c_int,
    pub stop_msg: String,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct nlopt_constraint {
    pub m: libc::c_uint,
    pub f: nlopt_func,
    pub mf: nlopt_mfunc,
    pub pre: nlopt_precond,
    pub f_data: *mut libc::c_void,
    pub tol: *mut libc::c_double,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct slsqpb_state {
    pub t: libc::c_double,
    pub f0: libc::c_double,
    pub h1: libc::c_double,
    pub h2: libc::c_double,
    pub h3: libc::c_double,
    pub h4: libc::c_double,
    pub n1: libc::c_int,
    pub n2: libc::c_int,
    pub n3: libc::c_int,
    pub t0: libc::c_double,
    pub gs: libc::c_double,
    pub tol: libc::c_double,
    pub line: libc::c_int,
    pub alpha: libc::c_double,
    pub iexact: libc::c_int,
    pub incons: libc::c_int,
    pub ireset: libc::c_int,
    pub itermx: libc::c_int,
    pub x0: *mut libc::c_double,
}
#[no_mangle]
pub unsafe fn nlopt_time_seed() -> libc::c_ulong {
    // let mut tv = libc::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // libc::gettimeofday(&mut tv, 0 as *mut libc::timezone);
    //return (tv.tv_sec ^ tv.tv_usec) as libc::c_ulong;
    let start = SystemTime::now();
    let since_the_epoch = start.duration_since(UNIX_EPOCH).expect("Time flies");
    since_the_epoch.as_millis() as libc::c_ulong
}
#[no_mangle]
pub unsafe fn nlopt_seconds() -> libc::c_double {
    // static mut start_inited: libc::c_int = 0 as libc::c_int;
    // static mut start: libc::timeval = libc::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // let mut tv: libc::timeval = libc::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // if start_inited == 0 {
    //     start_inited = 1 as libc::c_int;
    //     libc::gettimeofday(&mut start, 0 as *mut libc::timezone);
    // }
    // libc::gettimeofday(&mut tv, 0 as *mut libc::timezone);
    // return (tv.tv_sec - start.tv_sec) as libc::c_double
    //     + 1.0e-6f64 * (tv.tv_usec - start.tv_usec) as libc::c_double;
    static mut start_inited: bool = false;
    static mut start: SystemTime = UNIX_EPOCH;
    if !start_inited {
        start_inited = true;
        start = SystemTime::now();
    }
    start
        .duration_since(UNIX_EPOCH)
        .expect("Time flies")
        .as_secs_f64()
}
unsafe fn sc(
    mut x: libc::c_double,
    mut smin: libc::c_double,
    mut smax: libc::c_double,
) -> libc::c_double {
    return smin + x * (smax - smin);
}
unsafe fn vector_norm(
    mut n: libc::c_uint,
    mut vec: *const libc::c_double,
    mut w: *const libc::c_double,
    mut scale_min: *const libc::c_double,
    mut scale_max: *const libc::c_double,
) -> libc::c_double {
    let mut i: libc::c_uint = 0;
    let mut ret: libc::c_double = 0 as libc::c_int as libc::c_double;
    if !scale_min.is_null() && !scale_max.is_null() {
        if !w.is_null() {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += *w.offset(i as isize)
                    * (sc(
                        *vec.offset(i as isize),
                        *scale_min.offset(i as isize),
                        *scale_max.offset(i as isize),
                    ))
                    .abs();
                i = i.wrapping_add(1);
            }
        } else {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += (sc(
                    *vec.offset(i as isize),
                    *scale_min.offset(i as isize),
                    *scale_max.offset(i as isize),
                ))
                .abs();
                i = i.wrapping_add(1);
            }
        }
    } else if !w.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += *w.offset(i as isize) * (*vec.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += (*vec.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    }
    return ret;
}
unsafe fn diff_norm(
    mut n: libc::c_uint,
    mut x: *const libc::c_double,
    mut oldx: *const libc::c_double,
    mut w: *const libc::c_double,
    mut scale_min: *const libc::c_double,
    mut scale_max: *const libc::c_double,
) -> libc::c_double {
    let mut i: libc::c_uint = 0;
    let mut ret: libc::c_double = 0 as libc::c_int as libc::c_double;
    if !scale_min.is_null() && !scale_max.is_null() {
        if !w.is_null() {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += *w.offset(i as isize)
                    * (sc(
                        *x.offset(i as isize),
                        *scale_min.offset(i as isize),
                        *scale_max.offset(i as isize),
                    ) - sc(
                        *oldx.offset(i as isize),
                        *scale_min.offset(i as isize),
                        *scale_max.offset(i as isize),
                    ))
                    .abs();
                i = i.wrapping_add(1);
            }
        } else {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += (sc(
                    *x.offset(i as isize),
                    *scale_min.offset(i as isize),
                    *scale_max.offset(i as isize),
                ) - sc(
                    *oldx.offset(i as isize),
                    *scale_min.offset(i as isize),
                    *scale_max.offset(i as isize),
                ))
                .abs();
                i = i.wrapping_add(1);
            }
        }
    } else if !w.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += *w.offset(i as isize) * (*x.offset(i as isize) - *oldx.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += (*x.offset(i as isize) - *oldx.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    }
    return ret;
}
unsafe fn relstop(
    mut vold: libc::c_double,
    mut vnew: libc::c_double,
    mut reltol: libc::c_double,
    mut abstol: libc::c_double,
) -> libc::c_int {
    if nlopt_isinf(vold) != 0 {
        return 0 as libc::c_int;
    }
    return ((vnew - vold).abs() < abstol
        || (vnew - vold).abs() < reltol * ((vnew).abs() + (vold).abs()) * 0.5f64
        || reltol > 0 as libc::c_int as libc::c_double && vnew == vold) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_ftol(
    mut s: *const nlopt_stopping,
    mut f: libc::c_double,
    mut oldf: libc::c_double,
) -> libc::c_int {
    return relstop(oldf, f, (*s).ftol_rel, (*s).ftol_abs);
}
#[no_mangle]
pub unsafe fn nlopt_stop_f(
    mut s: *const nlopt_stopping,
    mut f: libc::c_double,
    mut oldf: libc::c_double,
) -> libc::c_int {
    return (f <= (*s).minf_max || nlopt_stop_ftol(s, f, oldf) != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_x(
    mut s: *const nlopt_stopping,
    mut x: *const libc::c_double,
    mut oldx: *const libc::c_double,
) -> libc::c_int {
    let mut i: libc::c_uint = 0;
    if diff_norm(
        (*s).n,
        x,
        oldx,
        (*s).x_weights,
        0 as *const libc::c_double,
        0 as *const libc::c_double,
    ) < (*s).xtol_rel
        * vector_norm(
            (*s).n,
            x,
            (*s).x_weights,
            0 as *const libc::c_double,
            0 as *const libc::c_double,
        )
    {
        return 1 as libc::c_int;
    }
    if ((*s).xtol_abs).is_null() {
        return 0 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < (*s).n {
        if (*x.offset(i as isize) - *oldx.offset(i as isize)).abs()
            >= *((*s).xtol_abs).offset(i as isize)
        {
            return 0 as libc::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_dx(
    mut s: *const nlopt_stopping,
    mut x: *const libc::c_double,
    mut dx: *const libc::c_double,
) -> libc::c_int {
    let mut i: libc::c_uint = 0;
    if vector_norm(
        (*s).n,
        dx,
        (*s).x_weights,
        0 as *const libc::c_double,
        0 as *const libc::c_double,
    ) < (*s).xtol_rel
        * vector_norm(
            (*s).n,
            x,
            (*s).x_weights,
            0 as *const libc::c_double,
            0 as *const libc::c_double,
        )
    {
        return 1 as libc::c_int;
    }
    if ((*s).xtol_abs).is_null() {
        return 0 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < (*s).n {
        if (*dx.offset(i as isize)).abs() >= *((*s).xtol_abs).offset(i as isize) {
            return 0 as libc::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_xs(
    mut s: *const nlopt_stopping,
    mut xs: *const libc::c_double,
    mut oldxs: *const libc::c_double,
    mut scale_min: *const libc::c_double,
    mut scale_max: *const libc::c_double,
) -> libc::c_int {
    let mut i: libc::c_uint = 0;
    if diff_norm((*s).n, xs, oldxs, (*s).x_weights, scale_min, scale_max)
        < (*s).xtol_rel * vector_norm((*s).n, xs, (*s).x_weights, scale_min, scale_max)
    {
        return 1 as libc::c_int;
    }
    if ((*s).xtol_abs).is_null() {
        return 0 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < (*s).n {
        if (sc(
            *xs.offset(i as isize),
            *scale_min.offset(i as isize),
            *scale_max.offset(i as isize),
        ) - sc(
            *oldxs.offset(i as isize),
            *scale_min.offset(i as isize),
            *scale_max.offset(i as isize),
        ))
        .abs()
            >= *((*s).xtol_abs).offset(i as isize)
        {
            return 0 as libc::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_isfinite(mut x: libc::c_double) -> libc::c_int {
    return ((x).abs() <= 1.7976931348623157e+308f64) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_istiny(mut x: libc::c_double) -> libc::c_int {
    if x == 0.0f64 {
        return 1 as libc::c_int;
    } else {
        return ((x).abs() < 2.2250738585072014e-308f64) as libc::c_int;
    };
}
#[no_mangle]
pub unsafe fn nlopt_isnan(mut x: libc::c_double) -> libc::c_int {
    return x.is_nan() as i32;
}
#[no_mangle]
pub unsafe fn nlopt_stop_evals(mut s: *const nlopt_stopping) -> libc::c_int {
    return ((*s).maxeval > 0 as libc::c_int && *(*s).nevals_p >= (*s).maxeval) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_time_(
    mut start: libc::c_double,
    mut maxtime: libc::c_double,
) -> libc::c_int {
    return (maxtime > 0 as libc::c_int as libc::c_double && nlopt_seconds() - start >= maxtime)
        as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_time(mut s: *const nlopt_stopping) -> libc::c_int {
    return nlopt_stop_time_((*s).start, (*s).maxtime);
}
#[no_mangle]
pub unsafe fn nlopt_stop_evalstime(mut stop: *const nlopt_stopping) -> libc::c_int {
    return (nlopt_stop_evals(stop) != 0 || nlopt_stop_time(stop) != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_forced(mut stop: *const nlopt_stopping) -> libc::c_int {
    return (!((*stop).force_stop).is_null() && *(*stop).force_stop != 0) as libc::c_int;
}
// #[no_mangle]
// pub unsafe fn nlopt_vsprintf(
//     mut p: *mut libc::c_char,
//     mut format: *const libc::c_char,
//     mut ap: ::std::ffi::VaList,
// ) -> *mut libc::c_char {
//     let mut len: size_t = (strlen(format))
//         .wrapping_add(128 as libc::c_int as libc::c_ulong);
//     let mut ret: libc::c_int = 0;
//     p = realloc(p as *mut libc::c_void, len) as *mut libc::c_char;
//     if p.is_null() {
//         abort();
//     }
//     loop {
//         ret = vsnprintf(p, len, format, ap.as_va_list());
//         if !(ret < 0 as libc::c_int || ret as size_t >= len) {
//             break;
//         }
//         len = if ret >= 0 as libc::c_int {
//             (ret + 1 as libc::c_int) as size_t
//         } else {
//             len.wrapping_mul(3 as libc::c_int as libc::c_ulong) >> 1 as libc::c_int
//         };
//         p = realloc(p as *mut libc::c_void, len) as *mut libc::c_char;
//         if p.is_null() {
//             abort();
//         }
//     }
//     return p;
// }
#[no_mangle]
pub unsafe fn nlopt_stop_msg(mut s: *mut nlopt_stopping, msg: &str) {
    (*s).stop_msg = msg.to_string();
}
#[no_mangle]
pub unsafe fn nlopt_count_constraints(
    mut p: libc::c_uint,
    mut c: *const nlopt_constraint,
) -> libc::c_uint {
    let mut i: libc::c_uint = 0;
    let mut count: libc::c_uint = 0 as libc::c_int as libc::c_uint;
    i = 0 as libc::c_int as libc::c_uint;
    while i < p {
        count = count.wrapping_add((*c.offset(i as isize)).m);
        i = i.wrapping_add(1);
    }
    return count;
}
#[no_mangle]
pub unsafe fn nlopt_max_constraint_dim(
    mut p: libc::c_uint,
    mut c: *const nlopt_constraint,
) -> libc::c_uint {
    let mut i: libc::c_uint = 0;
    let mut max_dim: libc::c_uint = 0 as libc::c_int as libc::c_uint;
    i = 0 as libc::c_int as libc::c_uint;
    while i < p {
        if (*c.offset(i as isize)).m > max_dim {
            max_dim = (*c.offset(i as isize)).m;
        }
        i = i.wrapping_add(1);
    }
    return max_dim;
}

pub unsafe fn nlopt_eval_constraint<U>(
    mut result: *mut libc::c_double,
    mut grad: *mut libc::c_double,
    mut c: *const nlopt_constraint,
    mut n: libc::c_uint,
    mut x: *const libc::c_double,
) {
    if ((*c).f).is_some() {
        *result.offset(0 as libc::c_int as isize) =
        // PATCH Weird bug ((*c).f).expect("non-null function pointer") calls the objective function!!!
        // even if (*c), nlopt_constraint object was correctly built with a nlopt_constraint_raw_callback!!! 
        //    ((*c).f).expect("non-null function pointer")(n, x, grad, (*c).f_data);
        // Maybe the U generic parameter required explains it cannot work like with C ???
        nlopt_constraint_raw_callback::<&dyn NLoptObjFn<U>, U>(n, x, grad, (*c).f_data);
    } else {
        ((*c).mf).expect("non-null function pointer")((*c).m, result, n, x, grad, (*c).f_data);
    };
}
#[no_mangle]
pub unsafe fn nlopt_isinf(mut x: libc::c_double) -> libc::c_int {
    return ((x).abs() >= ::std::f64::INFINITY * 0.99f64
        || if x.is_infinite() {
            if x.is_sign_positive() {
                1
            } else {
                -1
            }
        } else {
            0
        } != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_compute_rescaling(
    mut n: libc::c_uint,
    mut dx: *const libc::c_double,
) -> *mut libc::c_double {
    // let mut s: *mut libc::c_double = malloc(
    //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong).wrapping_mul(n as libc::c_ulong),
    // ) as *mut libc::c_double;

    let mut space: Box<Vec<libc::c_double>> = Box::new(vec![0.; usize::try_from(n).unwrap()]);
    let s = space.as_mut_ptr() as *mut libc::c_double;
    std::mem::forget(space);

    let mut i: libc::c_uint = 0;
    if s.is_null() {
        return 0 as *mut libc::c_double;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < n {
        *s.offset(i as isize) = 1.0f64;
        i = i.wrapping_add(1);
    }
    if n == 1 as libc::c_int as libc::c_uint {
        return s;
    }
    i = 1 as libc::c_int as libc::c_uint;
    while i < n
        && *dx.offset(i as isize)
            == *dx.offset(i.wrapping_sub(1 as libc::c_int as libc::c_uint) as isize)
    {
        i = i.wrapping_add(1);
    }
    if i < n {
        i = 1 as libc::c_int as libc::c_uint;
        while i < n {
            *s.offset(i as isize) = *dx.offset(i as isize) / *dx.offset(0 as libc::c_int as isize);
            i = i.wrapping_add(1);
        }
    }
    return s;
}
#[no_mangle]
pub unsafe fn nlopt_rescale(
    mut n: libc::c_uint,
    mut s: *const libc::c_double,
    mut x: *const libc::c_double,
    mut xs: *mut libc::c_double,
) {
    let mut i: libc::c_uint = 0;
    if s.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize);
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize) / *s.offset(i as isize);
            i = i.wrapping_add(1);
        }
    };
}
#[no_mangle]
pub unsafe fn nlopt_unscale(
    mut n: libc::c_uint,
    mut s: *const libc::c_double,
    mut x: *const libc::c_double,
    mut xs: *mut libc::c_double,
) {
    let mut i: libc::c_uint = 0;
    if s.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize);
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize) * *s.offset(i as isize);
            i = i.wrapping_add(1);
        }
    };
}
#[no_mangle]
pub unsafe fn nlopt_new_rescaled(
    mut n: libc::c_uint,
    mut s: *const libc::c_double,
    mut x: *const libc::c_double,
) -> *mut libc::c_double {
    // let mut xs: *mut libc::c_double = malloc(
    //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong).wrapping_mul(n as libc::c_ulong),
    // ) as *mut libc::c_double;

    let mut space: Box<Vec<libc::c_double>> = Box::new(vec![0.; usize::try_from(n).unwrap()]);
    let xs = space.as_mut_ptr() as *mut libc::c_double;
    std::mem::forget(space);

    if xs.is_null() {
        return 0 as *mut libc::c_double;
    }
    nlopt_rescale(n, s, x, xs);
    return xs;
}
#[no_mangle]
pub unsafe fn nlopt_reorder_bounds(
    mut n: libc::c_uint,
    mut lb: *mut libc::c_double,
    mut ub: *mut libc::c_double,
) {
    let mut i: libc::c_uint = 0;
    i = 0 as libc::c_int as libc::c_uint;
    while i < n {
        if *lb.offset(i as isize) > *ub.offset(i as isize) {
            let mut t: libc::c_double = *lb.offset(i as isize);
            *lb.offset(i as isize) = *ub.offset(i as isize);
            *ub.offset(i as isize) = t;
        }
        i = i.wrapping_add(1);
    }
}
unsafe fn dcopy___(
    mut n_: *mut libc::c_int,
    mut dx: *const libc::c_double,
    mut incx: libc::c_int,
    mut dy: *mut libc::c_double,
    mut incy: libc::c_int,
) {
    let mut i: libc::c_int = 0;
    let mut n: libc::c_int = *n_;
    if n <= 0 as libc::c_int {
        return;
    }
    if incx == 1 as libc::c_int && incy == 1 as libc::c_int {
        memcpy(
            dy as *mut libc::c_void,
            dx as *const libc::c_void,
            (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                .wrapping_mul(n as libc::c_uint as libc::c_ulong),
        );
    } else if incx == 0 as libc::c_int && incy == 1 as libc::c_int {
        let mut x: libc::c_double = *dx.offset(0 as libc::c_int as isize);
        i = 0 as libc::c_int;
        while i < n {
            *dy.offset(i as isize) = x;
            i += 1;
        }
    } else {
        i = 0 as libc::c_int;
        while i < n {
            *dy.offset((i * incy) as isize) = *dx.offset((i * incx) as isize);
            i += 1;
        }
    };
}
unsafe fn daxpy_sl__(
    mut n_: *mut libc::c_int,
    mut da_: *const libc::c_double,
    mut dx: *const libc::c_double,
    mut incx: libc::c_int,
    mut dy: *mut libc::c_double,
    mut incy: libc::c_int,
) {
    let mut n: libc::c_int = *n_;
    let mut i: libc::c_int = 0;
    let mut da: libc::c_double = *da_;
    if n <= 0 as libc::c_int || da == 0 as libc::c_int as libc::c_double {
        return;
    }
    i = 0 as libc::c_int;
    while i < n {
        *dy.offset((i * incy) as isize) += da * *dx.offset((i * incx) as isize);
        i += 1;
    }
}
unsafe fn ddot_sl__(
    mut n_: *mut libc::c_int,
    mut dx: *mut libc::c_double,
    mut incx: libc::c_int,
    mut dy: *mut libc::c_double,
    mut incy: libc::c_int,
) -> libc::c_double {
    let mut n: libc::c_int = *n_;
    let mut i: libc::c_int = 0;
    let mut sum: libc::c_double = 0 as libc::c_int as libc::c_double;
    if n <= 0 as libc::c_int {
        return 0 as libc::c_int as libc::c_double;
    }
    i = 0 as libc::c_int;
    while i < n {
        sum += *dx.offset((i * incx) as isize) * *dy.offset((i * incy) as isize);
        i += 1;
    }
    return sum;
}
unsafe fn dnrm2___(
    mut n_: *mut libc::c_int,
    mut dx: *mut libc::c_double,
    mut incx: libc::c_int,
) -> libc::c_double {
    let mut i: libc::c_int = 0;
    let mut n: libc::c_int = *n_;
    let mut xmax: libc::c_double = 0 as libc::c_int as libc::c_double;
    let mut scale: libc::c_double = 0.;
    let mut sum: libc::c_double = 0 as libc::c_int as libc::c_double;
    i = 0 as libc::c_int;
    while i < n {
        let mut xabs: libc::c_double = (*dx.offset((incx * i) as isize)).abs();
        if xmax < xabs {
            xmax = xabs;
        }
        i += 1;
    }
    if xmax == 0 as libc::c_int as libc::c_double {
        return 0 as libc::c_int as libc::c_double;
    }
    scale = 1.0f64 / xmax;
    i = 0 as libc::c_int;
    while i < n {
        let mut xs: libc::c_double = scale * *dx.offset((incx * i) as isize);
        sum += xs * xs;
        i += 1;
    }
    return xmax * sum.sqrt();
}
unsafe fn dsrot_(
    mut n: libc::c_int,
    mut dx: *mut libc::c_double,
    mut incx: libc::c_int,
    mut dy: *mut libc::c_double,
    mut incy: libc::c_int,
    mut c__: *mut libc::c_double,
    mut s_: *mut libc::c_double,
) {
    let mut i: libc::c_int = 0;
    let mut c: libc::c_double = *c__;
    let mut s: libc::c_double = *s_;
    i = 0 as libc::c_int;
    while i < n {
        let mut x: libc::c_double = *dx.offset((incx * i) as isize);
        let mut y: libc::c_double = *dy.offset((incy * i) as isize);
        *dx.offset((incx * i) as isize) = c * x + s * y;
        *dy.offset((incy * i) as isize) = c * y - s * x;
        i += 1;
    }
}
unsafe fn dsrotg_(
    mut da: *mut libc::c_double,
    mut db: *mut libc::c_double,
    mut c: *mut libc::c_double,
    mut s: *mut libc::c_double,
) {
    let mut absa: libc::c_double = 0.;
    let mut absb: libc::c_double = 0.;
    let mut roe: libc::c_double = 0.;
    let mut scale: libc::c_double = 0.;
    absa = (*da).abs();
    absb = (*db).abs();
    if absa > absb {
        roe = *da;
        scale = absa;
    } else {
        roe = *db;
        scale = absb;
    }
    if scale != 0 as libc::c_int as libc::c_double {
        let mut r: libc::c_double = 0.;
        let mut iscale: libc::c_double = 1 as libc::c_int as libc::c_double / scale;
        let mut tmpa: libc::c_double = *da * iscale;
        let mut tmpb: libc::c_double = *db * iscale;
        r = (if roe < 0 as libc::c_int as libc::c_double {
            -scale
        } else {
            scale
        }) * (tmpa * tmpa + tmpb * tmpb).sqrt();
        *c = *da / r;
        *s = *db / r;
        *da = r;
        if *c != 0 as libc::c_int as libc::c_double && (*c).abs() <= *s {
            *db = 1 as libc::c_int as libc::c_double / *c;
        } else {
            *db = *s;
        }
    } else {
        *c = 1 as libc::c_int as libc::c_double;
        *db = 0 as libc::c_int as libc::c_double;
        *da = *db;
        *s = *da;
    };
}
unsafe fn dscal_sl__(
    mut n_: *mut libc::c_int,
    mut da: *const libc::c_double,
    mut dx: *mut libc::c_double,
    mut incx: libc::c_int,
) {
    let mut i: libc::c_int = 0;
    let mut n: libc::c_int = *n_;
    let mut alpha: libc::c_double = *da;
    i = 0 as libc::c_int;
    while i < n {
        *dx.offset((i * incx) as isize) *= alpha;
        i += 1;
    }
}
static mut c__0: libc::c_int = 0 as libc::c_int;
static mut c__1: libc::c_int = 1 as libc::c_int;
static mut c__2: libc::c_int = 2 as libc::c_int;
unsafe fn h12_(
    mut mode: *const libc::c_int,
    mut lpivot: *mut libc::c_int,
    mut l1: *mut libc::c_int,
    mut m: *mut libc::c_int,
    mut u: *mut libc::c_double,
    mut iue: *const libc::c_int,
    mut up: *mut libc::c_double,
    mut c__: *mut libc::c_double,
    mut ice: *const libc::c_int,
    mut icv: *const libc::c_int,
    mut ncv: *const libc::c_int,
) {
    let mut current_block: u64;
    let one: libc::c_double = 1.0f64;
    let mut u_dim1: libc::c_int = 0;
    let mut u_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut b: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut i2: libc::c_int = 0;
    let mut i3: libc::c_int = 0;
    let mut i4: libc::c_int = 0;
    let mut cl: libc::c_double = 0.;
    let mut sm: libc::c_double = 0.;
    let mut incr: libc::c_int = 0;
    let mut clinv: libc::c_double = 0.;
    u_dim1 = *iue;
    u_offset = 1 as libc::c_int + u_dim1;
    u = u.offset(-(u_offset as isize));
    c__ = c__.offset(-1);
    if !(0 as libc::c_int >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
        d__1 = *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize);
        cl = (d__1).abs();
        if *mode == 2 as libc::c_int {
            if cl <= 0.0f64 {
                current_block = 12017783607140213051;
            } else {
                current_block = 16897903984108266201;
            }
        } else {
            i__1 = *m;
            j = *l1;
            while j <= i__1 {
                d__1 = *u.offset((j * u_dim1 + 1 as libc::c_int) as isize);
                sm = (d__1).abs();
                cl = if sm >= cl { sm } else { cl };
                j += 1;
            }
            if cl <= 0.0f64 {
                current_block = 12017783607140213051;
            } else {
                clinv = one / cl;
                d__1 = *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) * clinv;
                sm = d__1 * d__1;
                i__1 = *m;
                j = *l1;
                while j <= i__1 {
                    d__1 = *u.offset((j * u_dim1 + 1 as libc::c_int) as isize) * clinv;
                    sm += d__1 * d__1;
                    j += 1;
                }
                cl *= (sm).sqrt();
                if *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) > 0.0f64 {
                    cl = -cl;
                }
                *up = *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) - cl;
                *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) = cl;
                current_block = 16897903984108266201;
            }
        }
        match current_block {
            12017783607140213051 => {}
            _ => {
                if !(*ncv <= 0 as libc::c_int) {
                    b = *up * *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize);
                    if !(b >= 0.0f64) {
                        b = one / b;
                        i2 = 1 as libc::c_int - *icv + *ice * (*lpivot - 1 as libc::c_int);
                        incr = *ice * (*l1 - *lpivot);
                        i__1 = *ncv;
                        j = 1 as libc::c_int;
                        while j <= i__1 {
                            i2 += *icv;
                            i3 = i2 + incr;
                            i4 = i3;
                            sm = *c__.offset(i2 as isize) * *up;
                            i__2 = *m;
                            i__ = *l1;
                            while i__ <= i__2 {
                                sm += *c__.offset(i3 as isize)
                                    * *u.offset((i__ * u_dim1 + 1 as libc::c_int) as isize);
                                i3 += *ice;
                                i__ += 1;
                            }
                            if !(sm == 0.0f64) {
                                sm *= b;
                                *c__.offset(i2 as isize) += sm * *up;
                                i__2 = *m;
                                i__ = *l1;
                                while i__ <= i__2 {
                                    *c__.offset(i4 as isize) +=
                                        sm * *u.offset((i__ * u_dim1 + 1 as libc::c_int) as isize);
                                    i4 += *ice;
                                    i__ += 1;
                                }
                            }
                            j += 1;
                        }
                    }
                }
            }
        }
    }
}
unsafe fn nnls_(
    mut a: *mut libc::c_double,
    mut mda: *mut libc::c_int,
    mut m: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut b: *mut libc::c_double,
    mut x: *mut libc::c_double,
    mut rnorm: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut z__: *mut libc::c_double,
    mut indx: *mut libc::c_int,
    mut mode: *mut libc::c_int,
) {
    let mut current_block: u64;
    let one: libc::c_double = 1.0f64;
    let factor: libc::c_double = 0.01f64;
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut c__: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    let mut l: libc::c_int = 0;
    let mut s: libc::c_double = 0.;
    let mut t: libc::c_double = 0.;
    let mut ii: libc::c_int = 0;
    let mut jj: libc::c_int = 0;
    let mut ip: libc::c_int = 0;
    let mut iz: libc::c_int = 0;
    let mut jz: libc::c_int = 0;
    let mut up: libc::c_double = 0.;
    let mut iz1: libc::c_int = 0;
    let mut iz2: libc::c_int = 0;
    let mut npp1: libc::c_int = 0;
    let mut iter: libc::c_int = 0;
    let mut wmax: libc::c_double = 0.;
    let mut alpha: libc::c_double = 0.;
    let mut asave: libc::c_double = 0.;
    let mut itmax: libc::c_int = 0;
    let mut izmax: libc::c_int = 0 as libc::c_int;
    let mut nsetp: libc::c_int = 0;
    let mut unorm: libc::c_double = 0.;
    z__ = z__.offset(-1);
    b = b.offset(-1);
    indx = indx.offset(-1);
    w = w.offset(-1);
    x = x.offset(-1);
    a_dim1 = *mda;
    a_offset = 1 as libc::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    *mode = 2 as libc::c_int;
    if !(*m <= 0 as libc::c_int || *n <= 0 as libc::c_int) {
        *mode = 1 as libc::c_int;
        iter = 0 as libc::c_int;
        itmax = *n * 3 as libc::c_int;
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            *indx.offset(i__ as isize) = i__;
            i__ += 1;
        }
        iz1 = 1 as libc::c_int;
        iz2 = *n;
        nsetp = 0 as libc::c_int;
        npp1 = 1 as libc::c_int;
        *x.offset(1 as libc::c_int as isize) = 0.0f64;
        dcopy___(
            n,
            &*x.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *x.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        'c_11282: loop {
            if iz1 > iz2 || nsetp >= *m {
                current_block = 1795062611287366171;
                break;
            }
            i__1 = iz2;
            iz = iz1;
            while iz <= i__1 {
                j = *indx.offset(iz as isize);
                i__2 = *m - nsetp;
                *w.offset(j as isize) = ddot_sl__(
                    &mut i__2,
                    &mut *a.offset((npp1 + j * a_dim1) as isize),
                    1 as libc::c_int,
                    &mut *b.offset(npp1 as isize),
                    1 as libc::c_int,
                );
                iz += 1;
            }
            loop {
                wmax = 0.0f64;
                i__2 = iz2;
                iz = iz1;
                while iz <= i__2 {
                    j = *indx.offset(iz as isize);
                    if !(*w.offset(j as isize) <= wmax) {
                        wmax = *w.offset(j as isize);
                        izmax = iz;
                    }
                    iz += 1;
                }
                if wmax <= 0.0f64 {
                    current_block = 1795062611287366171;
                    break 'c_11282;
                }
                iz = izmax;
                j = *indx.offset(iz as isize);
                asave = *a.offset((npp1 + j * a_dim1) as isize);
                i__2 = npp1 + 1 as libc::c_int;
                h12_(
                    &c__1,
                    &mut npp1,
                    &mut i__2,
                    m,
                    &mut *a.offset((j * a_dim1 + 1 as libc::c_int) as isize),
                    &c__1,
                    &mut up,
                    &mut *z__.offset(1 as libc::c_int as isize),
                    &c__1,
                    &c__1,
                    &c__0,
                );
                unorm = dnrm2___(
                    &mut nsetp,
                    &mut *a.offset((j * a_dim1 + 1 as libc::c_int) as isize),
                    1 as libc::c_int,
                );
                d__1 = *a.offset((npp1 + j * a_dim1) as isize);
                t = factor * (d__1).abs();
                d__1 = unorm + t;
                if !(d__1 - unorm <= 0.0f64) {
                    dcopy___(
                        m,
                        &*b.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *z__.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    i__2 = npp1 + 1 as libc::c_int;
                    h12_(
                        &c__2,
                        &mut npp1,
                        &mut i__2,
                        m,
                        &mut *a.offset((j * a_dim1 + 1 as libc::c_int) as isize),
                        &c__1,
                        &mut up,
                        &mut *z__.offset(1 as libc::c_int as isize),
                        &c__1,
                        &c__1,
                        &c__1,
                    );
                    if *z__.offset(npp1 as isize) / *a.offset((npp1 + j * a_dim1) as isize) > 0.0f64
                    {
                        break;
                    }
                }
                *a.offset((npp1 + j * a_dim1) as isize) = asave;
                *w.offset(j as isize) = 0.0f64;
            }
            dcopy___(
                m,
                z__.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
                &mut *b.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
            );
            *indx.offset(iz as isize) = *indx.offset(iz1 as isize);
            *indx.offset(iz1 as isize) = j;
            iz1 += 1;
            nsetp = npp1;
            npp1 += 1;
            i__2 = iz2;
            jz = iz1;
            while jz <= i__2 {
                jj = *indx.offset(jz as isize);
                h12_(
                    &c__2,
                    &mut nsetp,
                    &mut npp1,
                    m,
                    &mut *a.offset((j * a_dim1 + 1 as libc::c_int) as isize),
                    &c__1,
                    &mut up,
                    &mut *a.offset((jj * a_dim1 + 1 as libc::c_int) as isize),
                    &c__1,
                    mda,
                    &c__1,
                );
                jz += 1;
            }
            k = if npp1 <= *mda { npp1 } else { *mda };
            *w.offset(j as isize) = 0.0f64;
            i__2 = *m - nsetp;
            dcopy___(
                &mut i__2,
                w.offset(j as isize),
                0 as libc::c_int,
                &mut *a.offset((k + j * a_dim1) as isize),
                1 as libc::c_int,
            );
            loop {
                ip = nsetp;
                while ip >= 1 as libc::c_int {
                    if !(ip == nsetp) {
                        d__1 = -*z__.offset((ip + 1 as libc::c_int) as isize);
                        daxpy_sl__(
                            &mut ip,
                            &d__1,
                            a.offset((jj * a_dim1 + 1 as libc::c_int) as isize),
                            1 as libc::c_int,
                            &mut *z__.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        );
                    }
                    jj = *indx.offset(ip as isize);
                    *z__.offset(ip as isize) /= *a.offset((ip + jj * a_dim1) as isize);
                    ip -= 1;
                }
                iter += 1;
                if !(iter <= itmax) {
                    current_block = 7781325733407116356;
                    break 'c_11282;
                }
                alpha = one;
                jj = 0 as libc::c_int;
                i__2 = nsetp;
                ip = 1 as libc::c_int;
                while ip <= i__2 {
                    if !(*z__.offset(ip as isize) > 0.0f64) {
                        l = *indx.offset(ip as isize);
                        t = -*x.offset(l as isize)
                            / (*z__.offset(ip as isize) - *x.offset(l as isize));
                        if !(alpha < t) {
                            alpha = t;
                            jj = ip;
                        }
                    }
                    ip += 1;
                }
                i__2 = nsetp;
                ip = 1 as libc::c_int;
                while ip <= i__2 {
                    l = *indx.offset(ip as isize);
                    *x.offset(l as isize) =
                        (one - alpha) * *x.offset(l as isize) + alpha * *z__.offset(ip as isize);
                    ip += 1;
                }
                if jj == 0 as libc::c_int {
                    break;
                }
                i__ = *indx.offset(jj as isize);
                'c_11330: loop {
                    *x.offset(i__ as isize) = 0.0f64;
                    jj += 1;
                    i__2 = nsetp;
                    j = jj;
                    while j <= i__2 {
                        ii = *indx.offset(j as isize);
                        *indx.offset((j - 1 as libc::c_int) as isize) = ii;
                        dsrotg_(
                            &mut *a.offset((j - 1 as libc::c_int + ii * a_dim1) as isize),
                            &mut *a.offset((j + ii * a_dim1) as isize),
                            &mut c__,
                            &mut s,
                        );
                        t = *a.offset((j - 1 as libc::c_int + ii * a_dim1) as isize);
                        dsrot_(
                            *n,
                            &mut *a.offset((j - 1 as libc::c_int + a_dim1) as isize),
                            *mda,
                            &mut *a.offset((j + a_dim1) as isize),
                            *mda,
                            &mut c__,
                            &mut s,
                        );
                        *a.offset((j - 1 as libc::c_int + ii * a_dim1) as isize) = t;
                        *a.offset((j + ii * a_dim1) as isize) = 0.0f64;
                        dsrot_(
                            1 as libc::c_int,
                            &mut *b.offset((j - 1 as libc::c_int) as isize),
                            1 as libc::c_int,
                            &mut *b.offset(j as isize),
                            1 as libc::c_int,
                            &mut c__,
                            &mut s,
                        );
                        j += 1;
                    }
                    npp1 = nsetp;
                    nsetp -= 1;
                    iz1 -= 1;
                    *indx.offset(iz1 as isize) = i__;
                    if nsetp <= 0 as libc::c_int {
                        current_block = 7781325733407116356;
                        break 'c_11282;
                    }
                    i__2 = nsetp;
                    jj = 1 as libc::c_int;
                    loop {
                        if !(jj <= i__2) {
                            break 'c_11330;
                        }
                        i__ = *indx.offset(jj as isize);
                        if *x.offset(i__ as isize) <= 0.0f64 {
                            break;
                        }
                        jj += 1;
                    }
                }
                dcopy___(
                    m,
                    b.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *z__.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
            }
        }
        match current_block {
            7781325733407116356 => {
                *mode = 3 as libc::c_int;
            }
            _ => {}
        }
        k = if npp1 <= *m { npp1 } else { *m };
        i__2 = *m - nsetp;
        *rnorm = dnrm2___(&mut i__2, &mut *b.offset(k as isize), 1 as libc::c_int);
        if npp1 > *m {
            *w.offset(1 as libc::c_int as isize) = 0.0f64;
            dcopy___(
                n,
                w.offset(1 as libc::c_int as isize),
                0 as libc::c_int,
                &mut *w.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
            );
        }
    }
}
unsafe fn ldp_(
    mut g: *mut libc::c_double,
    mut mg: *mut libc::c_int,
    mut m: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut h__: *mut libc::c_double,
    mut x: *mut libc::c_double,
    mut xnorm: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut indx: *mut libc::c_int,
    mut mode: *mut libc::c_int,
) {
    let one: libc::c_double = 1.0f64;
    let mut g_dim1: libc::c_int = 0;
    let mut g_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut n1: libc::c_int = 0;
    let mut if__: libc::c_int = 0;
    let mut iw: libc::c_int = 0;
    let mut iy: libc::c_int = 0;
    let mut iz: libc::c_int = 0;
    let mut fac: libc::c_double = 0.;
    let mut rnorm: libc::c_double = 0.;
    let mut iwdual: libc::c_int = 0;
    indx = indx.offset(-1);
    h__ = h__.offset(-1);
    x = x.offset(-1);
    g_dim1 = *mg;
    g_offset = 1 as libc::c_int + g_dim1;
    g = g.offset(-(g_offset as isize));
    w = w.offset(-1);
    *mode = 2 as libc::c_int;
    if !(*n <= 0 as libc::c_int) {
        *mode = 1 as libc::c_int;
        *x.offset(1 as libc::c_int as isize) = 0.0f64;
        dcopy___(
            n,
            x.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *x.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        *xnorm = 0.0f64;
        if !(*m == 0 as libc::c_int) {
            iw = 0 as libc::c_int;
            i__1 = *m;
            j = 1 as libc::c_int;
            while j <= i__1 {
                i__2 = *n;
                i__ = 1 as libc::c_int;
                while i__ <= i__2 {
                    iw += 1;
                    *w.offset(iw as isize) = *g.offset((j + i__ * g_dim1) as isize);
                    i__ += 1;
                }
                iw += 1;
                *w.offset(iw as isize) = *h__.offset(j as isize);
                j += 1;
            }
            if__ = iw + 1 as libc::c_int;
            i__1 = *n;
            i__ = 1 as libc::c_int;
            while i__ <= i__1 {
                iw += 1;
                *w.offset(iw as isize) = 0.0f64;
                i__ += 1;
            }
            *w.offset((iw + 1 as libc::c_int) as isize) = one;
            n1 = *n + 1 as libc::c_int;
            iz = iw + 2 as libc::c_int;
            iy = iz + n1;
            iwdual = iy + *m;
            nnls_(
                &mut *w.offset(1 as libc::c_int as isize),
                &mut n1,
                &mut n1,
                m,
                &mut *w.offset(if__ as isize),
                &mut *w.offset(iy as isize),
                &mut rnorm,
                &mut *w.offset(iwdual as isize),
                &mut *w.offset(iz as isize),
                &mut *indx.offset(1 as libc::c_int as isize),
                mode,
            );
            if !(*mode != 1 as libc::c_int) {
                *mode = 4 as libc::c_int;
                if !(rnorm <= 0.0f64) {
                    fac = one
                        - ddot_sl__(
                            m,
                            &mut *h__.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                            &mut *w.offset(iy as isize),
                            1 as libc::c_int,
                        );
                    d__1 = one + fac;
                    if !(d__1 - one <= 0.0f64) {
                        *mode = 1 as libc::c_int;
                        fac = one / fac;
                        i__1 = *n;
                        j = 1 as libc::c_int;
                        while j <= i__1 {
                            *x.offset(j as isize) = fac
                                * ddot_sl__(
                                    m,
                                    &mut *g.offset((j * g_dim1 + 1 as libc::c_int) as isize),
                                    1 as libc::c_int,
                                    &mut *w.offset(iy as isize),
                                    1 as libc::c_int,
                                );
                            j += 1;
                        }
                        *xnorm = dnrm2___(
                            n,
                            &mut *x.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        );
                        *w.offset(1 as libc::c_int as isize) = 0.0f64;
                        dcopy___(
                            m,
                            w.offset(1 as libc::c_int as isize),
                            0 as libc::c_int,
                            &mut *w.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        );
                        daxpy_sl__(
                            m,
                            &fac,
                            w.offset(iy as isize),
                            1 as libc::c_int,
                            &mut *w.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        );
                    }
                }
            }
        }
    }
}
unsafe fn lsi_(
    mut e: *mut libc::c_double,
    mut f: *mut libc::c_double,
    mut g: *mut libc::c_double,
    mut h__: *mut libc::c_double,
    mut le: *mut libc::c_int,
    mut me: *mut libc::c_int,
    mut lg: *mut libc::c_int,
    mut mg: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut x: *mut libc::c_double,
    mut xnorm: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut jw: *mut libc::c_int,
    mut mode: *mut libc::c_int,
) {
    let mut current_block: u64;
    let epmach: libc::c_double = 2.22e-16f64;
    let one: libc::c_double = 1.0f64;
    let mut e_dim1: libc::c_int = 0;
    let mut e_offset: libc::c_int = 0;
    let mut g_dim1: libc::c_int = 0;
    let mut g_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut i__3: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut t: libc::c_double = 0.;
    f = f.offset(-1);
    jw = jw.offset(-1);
    h__ = h__.offset(-1);
    x = x.offset(-1);
    g_dim1 = *lg;
    g_offset = 1 as libc::c_int + g_dim1;
    g = g.offset(-(g_offset as isize));
    e_dim1 = *le;
    e_offset = 1 as libc::c_int + e_dim1;
    e = e.offset(-(e_offset as isize));
    w = w.offset(-1);
    i__1 = *n;
    i__ = 1 as libc::c_int;
    while i__ <= i__1 {
        i__2 = i__ + 1 as libc::c_int;
        j = if i__2 <= *n { i__2 } else { *n };
        i__2 = i__ + 1 as libc::c_int;
        i__3 = *n - i__;
        h12_(
            &c__1,
            &mut i__,
            &mut i__2,
            me,
            &mut *e.offset((i__ * e_dim1 + 1 as libc::c_int) as isize),
            &c__1,
            &mut t,
            &mut *e.offset((j * e_dim1 + 1 as libc::c_int) as isize),
            &c__1,
            le,
            &i__3,
        );
        i__2 = i__ + 1 as libc::c_int;
        h12_(
            &c__2,
            &mut i__,
            &mut i__2,
            me,
            &mut *e.offset((i__ * e_dim1 + 1 as libc::c_int) as isize),
            &c__1,
            &mut t,
            &mut *f.offset(1 as libc::c_int as isize),
            &c__1,
            &c__1,
            &c__1,
        );
        i__ += 1;
    }
    *mode = 5 as libc::c_int;
    i__2 = *mg;
    i__ = 1 as libc::c_int;
    's_121: loop {
        if !(i__ <= i__2) {
            current_block = 14434620278749266018;
            break;
        }
        i__1 = *n;
        j = 1 as libc::c_int;
        while j <= i__1 {
            d__1 = *e.offset((j + j * e_dim1) as isize);
            if (d__1).abs() < epmach {
                current_block = 1714560199797219869;
                break 's_121;
            }
            i__3 = j - 1 as libc::c_int;
            *g.offset((i__ + j * g_dim1) as isize) = (*g.offset((i__ + j * g_dim1) as isize)
                - ddot_sl__(
                    &mut i__3,
                    &mut *g.offset((i__ + g_dim1) as isize),
                    *lg,
                    &mut *e.offset((j * e_dim1 + 1 as libc::c_int) as isize),
                    1 as libc::c_int,
                ))
                / *e.offset((j + j * e_dim1) as isize);
            j += 1;
        }
        *h__.offset(i__ as isize) -= ddot_sl__(
            n,
            &mut *g.offset((i__ + g_dim1) as isize),
            *lg,
            &mut *f.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        i__ += 1;
    }
    match current_block {
        14434620278749266018 => {
            ldp_(
                &mut *g.offset(g_offset as isize),
                lg,
                mg,
                n,
                &mut *h__.offset(1 as libc::c_int as isize),
                &mut *x.offset(1 as libc::c_int as isize),
                xnorm,
                &mut *w.offset(1 as libc::c_int as isize),
                &mut *jw.offset(1 as libc::c_int as isize),
                mode,
            );
            if !(*mode != 1 as libc::c_int) {
                daxpy_sl__(
                    n,
                    &one,
                    f.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                i__ = *n;
                while i__ >= 1 as libc::c_int {
                    i__2 = i__ + 1 as libc::c_int;
                    j = if i__2 <= *n { i__2 } else { *n };
                    i__2 = *n - i__;
                    *x.offset(i__ as isize) = (*x.offset(i__ as isize)
                        - ddot_sl__(
                            &mut i__2,
                            &mut *e.offset((i__ + j * e_dim1) as isize),
                            *le,
                            &mut *x.offset(j as isize),
                            1 as libc::c_int,
                        ))
                        / *e.offset((i__ + i__ * e_dim1) as isize);
                    i__ -= 1;
                }
                i__2 = *n + 1 as libc::c_int;
                j = if i__2 <= *me { i__2 } else { *me };
                i__2 = *me - *n;
                t = dnrm2___(&mut i__2, &mut *f.offset(j as isize), 1 as libc::c_int);
                *xnorm = (*xnorm * *xnorm + t * t).sqrt();
            }
        }
        _ => {}
    };
}
unsafe fn hfti_(
    mut a: *mut libc::c_double,
    mut mda: *mut libc::c_int,
    mut m: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut b: *mut libc::c_double,
    mut mdb: *mut libc::c_int,
    mut nb: *const libc::c_int,
    mut tau: *mut libc::c_double,
    mut krank: *mut libc::c_int,
    mut rnorm: *mut libc::c_double,
    mut h__: *mut libc::c_double,
    mut g: *mut libc::c_double,
    mut ip: *mut libc::c_int,
) {
    let mut current_block: u64;
    let factor: libc::c_double = 0.001f64;
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut b_dim1: libc::c_int = 0;
    let mut b_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut i__3: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    let mut l: libc::c_int = 0;
    let mut jb: libc::c_int = 0;
    let mut kp1: libc::c_int = 0;
    let mut tmp: libc::c_double = 0.;
    let mut hmax: libc::c_double = 0.;
    let mut lmax: libc::c_int = 0;
    let mut ldiag: libc::c_int = 0;
    ip = ip.offset(-1);
    g = g.offset(-1);
    h__ = h__.offset(-1);
    a_dim1 = *mda;
    a_offset = 1 as libc::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    rnorm = rnorm.offset(-1);
    b_dim1 = *mdb;
    b_offset = 1 as libc::c_int + b_dim1;
    b = b.offset(-(b_offset as isize));
    k = 0 as libc::c_int;
    ldiag = if *m <= *n { *m } else { *n };
    if !(ldiag <= 0 as libc::c_int) {
        i__1 = ldiag;
        j = 1 as libc::c_int;
        while j <= i__1 {
            let mut current_block_56: u64;
            if j == 1 as libc::c_int {
                current_block_56 = 18323326546081147283;
            } else {
                lmax = j;
                i__2 = *n;
                l = j;
                while l <= i__2 {
                    d__1 = *a.offset((j - 1 as libc::c_int + l * a_dim1) as isize);
                    *h__.offset(l as isize) -= d__1 * d__1;
                    if *h__.offset(l as isize) > *h__.offset(lmax as isize) {
                        lmax = l;
                    }
                    l += 1;
                }
                d__1 = hmax + factor * *h__.offset(lmax as isize);
                if d__1 - hmax > 0.0f64 {
                    current_block_56 = 17342314209269636488;
                } else {
                    current_block_56 = 18323326546081147283;
                }
            }
            match current_block_56 {
                18323326546081147283 => {
                    lmax = j;
                    i__2 = *n;
                    l = j;
                    while l <= i__2 {
                        *h__.offset(l as isize) = 0.0f64;
                        i__3 = *m;
                        i__ = j;
                        while i__ <= i__3 {
                            d__1 = *a.offset((i__ + l * a_dim1) as isize);
                            *h__.offset(l as isize) += d__1 * d__1;
                            i__ += 1;
                        }
                        if *h__.offset(l as isize) > *h__.offset(lmax as isize) {
                            lmax = l;
                        }
                        l += 1;
                    }
                    hmax = *h__.offset(lmax as isize);
                }
                _ => {}
            }
            *ip.offset(j as isize) = lmax;
            if !(*ip.offset(j as isize) == j) {
                i__2 = *m;
                i__ = 1 as libc::c_int;
                while i__ <= i__2 {
                    tmp = *a.offset((i__ + j * a_dim1) as isize);
                    *a.offset((i__ + j * a_dim1) as isize) =
                        *a.offset((i__ + lmax * a_dim1) as isize);
                    *a.offset((i__ + lmax * a_dim1) as isize) = tmp;
                    i__ += 1;
                }
                *h__.offset(lmax as isize) = *h__.offset(j as isize);
            }
            i__2 = j + 1 as libc::c_int;
            i__ = if i__2 <= *n { i__2 } else { *n };
            i__2 = j + 1 as libc::c_int;
            i__3 = *n - j;
            h12_(
                &c__1,
                &mut j,
                &mut i__2,
                m,
                &mut *a.offset((j * a_dim1 + 1 as libc::c_int) as isize),
                &c__1,
                &mut *h__.offset(j as isize),
                &mut *a.offset((i__ * a_dim1 + 1 as libc::c_int) as isize),
                &c__1,
                mda,
                &i__3,
            );
            i__2 = j + 1 as libc::c_int;
            h12_(
                &c__2,
                &mut j,
                &mut i__2,
                m,
                &mut *a.offset((j * a_dim1 + 1 as libc::c_int) as isize),
                &c__1,
                &mut *h__.offset(j as isize),
                &mut *b.offset(b_offset as isize),
                &c__1,
                mdb,
                nb,
            );
            j += 1;
        }
        i__2 = ldiag;
        j = 1 as libc::c_int;
        loop {
            if !(j <= i__2) {
                current_block = 18038362259723567392;
                break;
            }
            d__1 = *a.offset((j + j * a_dim1) as isize);
            if (d__1).abs() <= *tau {
                current_block = 11399179239179359665;
                break;
            }
            j += 1;
        }
        match current_block {
            18038362259723567392 => {
                k = ldiag;
            }
            _ => {
                k = j - 1 as libc::c_int;
            }
        }
        kp1 = k + 1 as libc::c_int;
        i__2 = *nb;
        jb = 1 as libc::c_int;
        while jb <= i__2 {
            i__1 = *m - k;
            *rnorm.offset(jb as isize) = dnrm2___(
                &mut i__1,
                &mut *b.offset((kp1 + jb * b_dim1) as isize),
                1 as libc::c_int,
            );
            jb += 1;
        }
        if k > 0 as libc::c_int {
            if !(k == *n) {
                i__ = k;
                while i__ >= 1 as libc::c_int {
                    i__2 = i__ - 1 as libc::c_int;
                    h12_(
                        &c__1,
                        &mut i__,
                        &mut kp1,
                        n,
                        &mut *a.offset((i__ + a_dim1) as isize),
                        mda,
                        &mut *g.offset(i__ as isize),
                        &mut *a.offset(a_offset as isize),
                        mda,
                        &c__1,
                        &i__2,
                    );
                    i__ -= 1;
                }
            }
            i__2 = *nb;
            jb = 1 as libc::c_int;
            while jb <= i__2 {
                i__ = k;
                while i__ >= 1 as libc::c_int {
                    i__1 = i__ + 1 as libc::c_int;
                    j = if i__1 <= *n { i__1 } else { *n };
                    i__1 = k - i__;
                    *b.offset((i__ + jb * b_dim1) as isize) = (*b
                        .offset((i__ + jb * b_dim1) as isize)
                        - ddot_sl__(
                            &mut i__1,
                            &mut *a.offset((i__ + j * a_dim1) as isize),
                            *mda,
                            &mut *b.offset((j + jb * b_dim1) as isize),
                            1 as libc::c_int,
                        ))
                        / *a.offset((i__ + i__ * a_dim1) as isize);
                    i__ -= 1;
                }
                if !(k == *n) {
                    i__1 = *n;
                    j = kp1;
                    while j <= i__1 {
                        *b.offset((j + jb * b_dim1) as isize) = 0.0f64;
                        j += 1;
                    }
                    i__1 = k;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        h12_(
                            &c__2,
                            &mut i__,
                            &mut kp1,
                            n,
                            &mut *a.offset((i__ + a_dim1) as isize),
                            mda,
                            &mut *g.offset(i__ as isize),
                            &mut *b.offset((jb * b_dim1 + 1 as libc::c_int) as isize),
                            &c__1,
                            mdb,
                            &c__1,
                        );
                        i__ += 1;
                    }
                }
                j = ldiag;
                while j >= 1 as libc::c_int {
                    if !(*ip.offset(j as isize) == j) {
                        l = *ip.offset(j as isize);
                        tmp = *b.offset((l + jb * b_dim1) as isize);
                        *b.offset((l + jb * b_dim1) as isize) =
                            *b.offset((j + jb * b_dim1) as isize);
                        *b.offset((j + jb * b_dim1) as isize) = tmp;
                    }
                    j -= 1;
                }
                jb += 1;
            }
        } else {
            i__1 = *nb;
            jb = 1 as libc::c_int;
            while jb <= i__1 {
                i__2 = *n;
                i__ = 1 as libc::c_int;
                while i__ <= i__2 {
                    *b.offset((i__ + jb * b_dim1) as isize) = 0.0f64;
                    i__ += 1;
                }
                jb += 1;
            }
        }
    }
    *krank = k;
}
unsafe fn lsei_(
    mut c__: *mut libc::c_double,
    mut d__: *mut libc::c_double,
    mut e: *mut libc::c_double,
    mut f: *mut libc::c_double,
    mut g: *mut libc::c_double,
    mut h__: *mut libc::c_double,
    mut lc: *mut libc::c_int,
    mut mc: *mut libc::c_int,
    mut le: *mut libc::c_int,
    mut me: *mut libc::c_int,
    mut lg: *mut libc::c_int,
    mut mg: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut x: *mut libc::c_double,
    mut xnrm: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut jw: *mut libc::c_int,
    mut mode: *mut libc::c_int,
) {
    let mut current_block: u64;
    let epmach: libc::c_double = 2.22e-16f64;
    let mut c_dim1: libc::c_int = 0;
    let mut c_offset: libc::c_int = 0;
    let mut e_dim1: libc::c_int = 0;
    let mut e_offset: libc::c_int = 0;
    let mut g_dim1: libc::c_int = 0;
    let mut g_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut i__3: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    let mut l: libc::c_int = 0;
    let mut t: libc::c_double = 0.;
    let mut ie: libc::c_int = 0;
    let mut if__: libc::c_int = 0;
    let mut ig: libc::c_int = 0;
    let mut iw: libc::c_int = 0;
    let mut mc1: libc::c_int = 0;
    let mut krank: libc::c_int = 0;
    d__ = d__.offset(-1);
    f = f.offset(-1);
    h__ = h__.offset(-1);
    x = x.offset(-1);
    g_dim1 = *lg;
    g_offset = 1 as libc::c_int + g_dim1;
    g = g.offset(-(g_offset as isize));
    e_dim1 = *le;
    e_offset = 1 as libc::c_int + e_dim1;
    e = e.offset(-(e_offset as isize));
    c_dim1 = *lc;
    c_offset = 1 as libc::c_int + c_dim1;
    c__ = c__.offset(-(c_offset as isize));
    w = w.offset(-1);
    jw = jw.offset(-1);
    *mode = 2 as libc::c_int;
    if !(*mc > *n) {
        l = *n - *mc;
        mc1 = *mc + 1 as libc::c_int;
        iw = (l + 1 as libc::c_int) * (*mg + 2 as libc::c_int) + (*mg << 1 as libc::c_int) + *mc;
        ie = iw + *mc + 1 as libc::c_int;
        if__ = ie + *me * l;
        ig = if__ + *me;
        i__1 = *mc;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            i__2 = i__ + 1 as libc::c_int;
            j = if i__2 <= *lc { i__2 } else { *lc };
            i__2 = i__ + 1 as libc::c_int;
            i__3 = *mc - i__;
            h12_(
                &c__1,
                &mut i__,
                &mut i__2,
                n,
                &mut *c__.offset((i__ + c_dim1) as isize),
                lc,
                &mut *w.offset((iw + i__) as isize),
                &mut *c__.offset((j + c_dim1) as isize),
                lc,
                &c__1,
                &i__3,
            );
            i__2 = i__ + 1 as libc::c_int;
            h12_(
                &c__2,
                &mut i__,
                &mut i__2,
                n,
                &mut *c__.offset((i__ + c_dim1) as isize),
                lc,
                &mut *w.offset((iw + i__) as isize),
                &mut *e.offset(e_offset as isize),
                le,
                &c__1,
                me,
            );
            i__2 = i__ + 1 as libc::c_int;
            h12_(
                &c__2,
                &mut i__,
                &mut i__2,
                n,
                &mut *c__.offset((i__ + c_dim1) as isize),
                lc,
                &mut *w.offset((iw + i__) as isize),
                &mut *g.offset(g_offset as isize),
                lg,
                &c__1,
                mg,
            );
            i__ += 1;
        }
        *mode = 6 as libc::c_int;
        i__2 = *mc;
        i__ = 1 as libc::c_int;
        loop {
            if !(i__ <= i__2) {
                current_block = 3222590281903869779;
                break;
            }
            d__1 = *c__.offset((i__ + i__ * c_dim1) as isize);
            if (d__1).abs() < epmach {
                current_block = 1640519170444703904;
                break;
            }
            i__1 = i__ - 1 as libc::c_int;
            *x.offset(i__ as isize) = (*d__.offset(i__ as isize)
                - ddot_sl__(
                    &mut i__1,
                    &mut *c__.offset((i__ + c_dim1) as isize),
                    *lc,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                ))
                / *c__.offset((i__ + i__ * c_dim1) as isize);
            i__ += 1;
        }
        match current_block {
            1640519170444703904 => {}
            _ => {
                *mode = 1 as libc::c_int;
                *w.offset(mc1 as isize) = 0.0f64;
                i__2 = *mg;
                dcopy___(
                    &mut i__2,
                    w.offset(mc1 as isize),
                    0 as libc::c_int,
                    &mut *w.offset(mc1 as isize),
                    1 as libc::c_int,
                );
                if *mc == *n {
                    current_block = 10769700236867903670;
                } else {
                    i__2 = *me;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__2 {
                        *w.offset((if__ - 1 as libc::c_int + i__) as isize) = *f
                            .offset(i__ as isize)
                            - ddot_sl__(
                                mc,
                                &mut *e.offset((i__ + e_dim1) as isize),
                                *le,
                                &mut *x.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            );
                        i__ += 1;
                    }
                    i__2 = *me;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__2 {
                        dcopy___(
                            &mut l,
                            e.offset((i__ + mc1 * e_dim1) as isize),
                            *le,
                            &mut *w.offset((ie - 1 as libc::c_int + i__) as isize),
                            *me,
                        );
                        i__ += 1;
                    }
                    i__2 = *mg;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__2 {
                        dcopy___(
                            &mut l,
                            g.offset((i__ + mc1 * g_dim1) as isize),
                            *lg,
                            &mut *w.offset((ig - 1 as libc::c_int + i__) as isize),
                            *mg,
                        );
                        i__ += 1;
                    }
                    if *mg > 0 as libc::c_int {
                        i__2 = *mg;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__2 {
                            *h__.offset(i__ as isize) -= ddot_sl__(
                                mc,
                                &mut *g.offset((i__ + g_dim1) as isize),
                                *lg,
                                &mut *x.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            );
                            i__ += 1;
                        }
                        lsi_(
                            &mut *w.offset(ie as isize),
                            &mut *w.offset(if__ as isize),
                            &mut *w.offset(ig as isize),
                            &mut *h__.offset(1 as libc::c_int as isize),
                            me,
                            me,
                            mg,
                            mg,
                            &mut l,
                            &mut *x.offset(mc1 as isize),
                            xnrm,
                            &mut *w.offset(mc1 as isize),
                            &mut *jw.offset(1 as libc::c_int as isize),
                            mode,
                        );
                        if *mc == 0 as libc::c_int {
                            current_block = 1640519170444703904;
                        } else {
                            t = dnrm2___(
                                mc,
                                &mut *x.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            );
                            *xnrm = (*xnrm * *xnrm + t * t).sqrt();
                            if *mode != 1 as libc::c_int {
                                current_block = 1640519170444703904;
                            } else {
                                current_block = 10769700236867903670;
                            }
                        }
                    } else {
                        *mode = 7 as libc::c_int;
                        k = if *le >= *n { *le } else { *n };
                        t = (epmach).sqrt();
                        hfti_(
                            &mut *w.offset(ie as isize),
                            me,
                            me,
                            &mut l,
                            &mut *w.offset(if__ as isize),
                            &mut k,
                            &c__1,
                            &mut t,
                            &mut krank,
                            xnrm,
                            &mut *w.offset(1 as libc::c_int as isize),
                            &mut *w.offset((l + 1 as libc::c_int) as isize),
                            &mut *jw.offset(1 as libc::c_int as isize),
                        );
                        dcopy___(
                            &mut l,
                            w.offset(if__ as isize),
                            1 as libc::c_int,
                            &mut *x.offset(mc1 as isize),
                            1 as libc::c_int,
                        );
                        if krank != l {
                            current_block = 1640519170444703904;
                        } else {
                            *mode = 1 as libc::c_int;
                            current_block = 10769700236867903670;
                        }
                    }
                }
                match current_block {
                    1640519170444703904 => {}
                    _ => {
                        i__2 = *me;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__2 {
                            *f.offset(i__ as isize) = ddot_sl__(
                                n,
                                &mut *e.offset((i__ + e_dim1) as isize),
                                *le,
                                &mut *x.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            ) - *f.offset(i__ as isize);
                            i__ += 1;
                        }
                        i__2 = *mc;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__2 {
                            *d__.offset(i__ as isize) = ddot_sl__(
                                me,
                                &mut *e.offset((i__ * e_dim1 + 1 as libc::c_int) as isize),
                                1 as libc::c_int,
                                &mut *f.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            ) - ddot_sl__(
                                mg,
                                &mut *g.offset((i__ * g_dim1 + 1 as libc::c_int) as isize),
                                1 as libc::c_int,
                                &mut *w.offset(mc1 as isize),
                                1 as libc::c_int,
                            );
                            i__ += 1;
                        }
                        i__ = *mc;
                        while i__ >= 1 as libc::c_int {
                            i__2 = i__ + 1 as libc::c_int;
                            h12_(
                                &c__2,
                                &mut i__,
                                &mut i__2,
                                n,
                                &mut *c__.offset((i__ + c_dim1) as isize),
                                lc,
                                &mut *w.offset((iw + i__) as isize),
                                &mut *x.offset(1 as libc::c_int as isize),
                                &c__1,
                                &c__1,
                                &c__1,
                            );
                            i__ -= 1;
                        }
                        i__ = *mc;
                        while i__ >= 1 as libc::c_int {
                            i__2 = i__ + 1 as libc::c_int;
                            j = if i__2 <= *lc { i__2 } else { *lc };
                            i__2 = *mc - i__;
                            *w.offset(i__ as isize) = (*d__.offset(i__ as isize)
                                - ddot_sl__(
                                    &mut i__2,
                                    &mut *c__.offset((j + i__ * c_dim1) as isize),
                                    1 as libc::c_int,
                                    &mut *w.offset(j as isize),
                                    1 as libc::c_int,
                                ))
                                / *c__.offset((i__ + i__ * c_dim1) as isize);
                            i__ -= 1;
                        }
                    }
                }
            }
        }
    }
}
unsafe fn lsq_(
    mut m: *mut libc::c_int,
    mut meq: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut nl: *mut libc::c_int,
    mut la: *mut libc::c_int,
    mut l: *mut libc::c_double,
    mut g: *mut libc::c_double,
    mut a: *mut libc::c_double,
    mut b: *mut libc::c_double,
    mut xl: *const libc::c_double,
    mut xu: *const libc::c_double,
    mut x: *mut libc::c_double,
    mut y: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut jw: *mut libc::c_int,
    mut mode: *mut libc::c_int,
) {
    let one: libc::c_double = 1.0f64;
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut i1: libc::c_int = 0;
    let mut i2: libc::c_int = 0;
    let mut i3: libc::c_int = 0;
    let mut i4: libc::c_int = 0;
    let mut m1: libc::c_int = 0;
    let mut n1: libc::c_int = 0;
    let mut n2: libc::c_int = 0;
    let mut n3: libc::c_int = 0;
    let mut ic: libc::c_int = 0;
    let mut id: libc::c_int = 0;
    let mut ie: libc::c_int = 0;
    let mut if__: libc::c_int = 0;
    let mut ig: libc::c_int = 0;
    let mut ih: libc::c_int = 0;
    let mut il: libc::c_int = 0;
    let mut im: libc::c_int = 0;
    let mut ip: libc::c_int = 0;
    let mut iu: libc::c_int = 0;
    let mut iw: libc::c_int = 0;
    let mut diag: libc::c_double = 0.;
    let mut mineq: libc::c_int = 0;
    let mut xnorm: libc::c_double = 0.;
    y = y.offset(-1);
    x = x.offset(-1);
    xu = xu.offset(-1);
    xl = xl.offset(-1);
    g = g.offset(-1);
    l = l.offset(-1);
    b = b.offset(-1);
    a_dim1 = *la;
    a_offset = 1 as libc::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    w = w.offset(-1);
    jw = jw.offset(-1);
    n1 = *n + 1 as libc::c_int;
    mineq = *m - *meq;
    m1 = mineq + *n + *n;
    n2 = n1 * *n / 2 as libc::c_int + 1 as libc::c_int;
    if n2 == *nl {
        n2 = 0 as libc::c_int;
    } else {
        n2 = 1 as libc::c_int;
    }
    n3 = *n - n2;
    i2 = 1 as libc::c_int;
    i3 = 1 as libc::c_int;
    i4 = 1 as libc::c_int;
    ie = 1 as libc::c_int;
    if__ = *n * *n + 1 as libc::c_int;
    i__1 = n3;
    i__ = 1 as libc::c_int;
    while i__ <= i__1 {
        i1 = n1 - i__;
        diag = (*l.offset(i2 as isize)).sqrt();
        *w.offset(i3 as isize) = 0.0f64;
        dcopy___(
            &mut i1,
            w.offset(i3 as isize),
            0 as libc::c_int,
            &mut *w.offset(i3 as isize),
            1 as libc::c_int,
        );
        i__2 = i1 - n2;
        dcopy___(
            &mut i__2,
            l.offset(i2 as isize),
            1 as libc::c_int,
            &mut *w.offset(i3 as isize),
            *n,
        );
        i__2 = i1 - n2;
        dscal_sl__(&mut i__2, &diag, &mut *w.offset(i3 as isize), *n);
        *w.offset(i3 as isize) = diag;
        i__2 = i__ - 1 as libc::c_int;
        *w.offset((if__ - 1 as libc::c_int + i__) as isize) = (*g.offset(i__ as isize)
            - ddot_sl__(
                &mut i__2,
                &mut *w.offset(i4 as isize),
                1 as libc::c_int,
                &mut *w.offset(if__ as isize),
                1 as libc::c_int,
            ))
            / diag;
        i2 = i2 + i1 - n2;
        i3 += n1;
        i4 += *n;
        i__ += 1;
    }
    if n2 == 1 as libc::c_int {
        *w.offset(i3 as isize) = *l.offset(*nl as isize);
        *w.offset(i4 as isize) = 0.0f64;
        dcopy___(
            &mut n3,
            w.offset(i4 as isize),
            0 as libc::c_int,
            &mut *w.offset(i4 as isize),
            1 as libc::c_int,
        );
        *w.offset((if__ - 1 as libc::c_int + *n) as isize) = 0.0f64;
    }
    d__1 = -one;
    dscal_sl__(n, &d__1, &mut *w.offset(if__ as isize), 1 as libc::c_int);
    ic = if__ + *n;
    id = ic + *meq * *n;
    if *meq > 0 as libc::c_int {
        i__1 = *meq;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            dcopy___(
                n,
                a.offset((i__ + a_dim1) as isize),
                *la,
                &mut *w.offset((ic - 1 as libc::c_int + i__) as isize),
                *meq,
            );
            i__ += 1;
        }
        dcopy___(
            meq,
            b.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
            &mut *w.offset(id as isize),
            1 as libc::c_int,
        );
        d__1 = -one;
        dscal_sl__(meq, &d__1, &mut *w.offset(id as isize), 1 as libc::c_int);
    }
    ig = id + *meq;
    if mineq > 0 as libc::c_int {
        i__1 = mineq;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            dcopy___(
                n,
                a.offset((*meq + i__ + a_dim1) as isize),
                *la,
                &mut *w.offset((ig - 1 as libc::c_int + i__) as isize),
                m1,
            );
            i__ += 1;
        }
    }
    ip = ig + mineq;
    i__1 = *n;
    i__ = 1 as libc::c_int;
    while i__ <= i__1 {
        *w.offset((ip - 1 as libc::c_int + i__) as isize) = 0.0f64;
        dcopy___(
            n,
            w.offset((ip - 1 as libc::c_int + i__) as isize),
            0 as libc::c_int,
            &mut *w.offset((ip - 1 as libc::c_int + i__) as isize),
            m1,
        );
        i__ += 1;
    }
    i__1 = m1 + 1 as libc::c_int;
    i__ = 1 as libc::c_int;
    while i__ <= *n {
        if nlopt_isinf(*xl.offset(i__ as isize)) == 0 {
            *w.offset((ip - i__1 + i__ * i__1) as isize) = 1.0f64;
        }
        i__ += 1;
    }
    im = ip + *n;
    i__1 = *n;
    i__ = 1 as libc::c_int;
    while i__ <= i__1 {
        *w.offset((im - 1 as libc::c_int + i__) as isize) = 0.0f64;
        dcopy___(
            n,
            w.offset((im - 1 as libc::c_int + i__) as isize),
            0 as libc::c_int,
            &mut *w.offset((im - 1 as libc::c_int + i__) as isize),
            m1,
        );
        i__ += 1;
    }
    i__1 = m1 + 1 as libc::c_int;
    i__ = 1 as libc::c_int;
    while i__ <= *n {
        if nlopt_isinf(*xu.offset(i__ as isize)) == 0 {
            *w.offset((im - i__1 + i__ * i__1) as isize) = -1.0f64;
        }
        i__ += 1;
    }
    ih = ig + m1 * *n;
    if mineq > 0 as libc::c_int {
        dcopy___(
            &mut mineq,
            b.offset((*meq + 1 as libc::c_int) as isize),
            1 as libc::c_int,
            &mut *w.offset(ih as isize),
            1 as libc::c_int,
        );
        d__1 = -one;
        dscal_sl__(
            &mut mineq,
            &d__1,
            &mut *w.offset(ih as isize),
            1 as libc::c_int,
        );
    }
    il = ih + mineq;
    iu = il + *n;
    i__ = 1 as libc::c_int;
    while i__ <= *n {
        *w.offset((il - 1 as libc::c_int + i__) as isize) =
            if nlopt_isinf(*xl.offset(i__ as isize)) != 0 {
                0 as libc::c_int as libc::c_double
            } else {
                *xl.offset(i__ as isize)
            };
        *w.offset((iu - 1 as libc::c_int + i__) as isize) =
            if nlopt_isinf(*xu.offset(i__ as isize)) != 0 {
                0 as libc::c_int as libc::c_double
            } else {
                -*xu.offset(i__ as isize)
            };
        i__ += 1;
    }
    iw = iu + *n;
    i__1 = if 1 as libc::c_int >= *meq {
        1 as libc::c_int
    } else {
        *meq
    };
    lsei_(
        &mut *w.offset(ic as isize),
        &mut *w.offset(id as isize),
        &mut *w.offset(ie as isize),
        &mut *w.offset(if__ as isize),
        &mut *w.offset(ig as isize),
        &mut *w.offset(ih as isize),
        &mut i__1,
        meq,
        n,
        n,
        &mut m1,
        &mut m1,
        n,
        &mut *x.offset(1 as libc::c_int as isize),
        &mut xnorm,
        &mut *w.offset(iw as isize),
        &mut *jw.offset(1 as libc::c_int as isize),
        mode,
    );
    if *mode == 1 as libc::c_int {
        dcopy___(
            m,
            w.offset(iw as isize),
            1 as libc::c_int,
            &mut *y.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        dcopy___(
            &mut n3,
            w.offset((iw + *m) as isize),
            1 as libc::c_int,
            &mut *y.offset((*m + 1 as libc::c_int) as isize),
            1 as libc::c_int,
        );
        dcopy___(
            &mut n3,
            w.offset((iw + *m + *n) as isize),
            1 as libc::c_int,
            &mut *y.offset((*m + n3 + 1 as libc::c_int) as isize),
            1 as libc::c_int,
        );
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            if *x.offset(i__ as isize) < *xl.offset(i__ as isize) {
                *x.offset(i__ as isize) = *xl.offset(i__ as isize);
            } else if *x.offset(i__ as isize) > *xu.offset(i__ as isize) {
                *x.offset(i__ as isize) = *xu.offset(i__ as isize);
            }
            i__ += 1;
        }
    }
}
unsafe fn ldl_(
    mut n: *mut libc::c_int,
    mut a: *mut libc::c_double,
    mut z__: *mut libc::c_double,
    mut sigma: *mut libc::c_double,
    mut w: *mut libc::c_double,
) {
    let one: libc::c_double = 1.0f64;
    let four: libc::c_double = 4.0f64;
    let epmach: libc::c_double = 2.22e-16f64;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut t: libc::c_double = 0.;
    let mut u: libc::c_double = 0.;
    let mut v: libc::c_double = 0.;
    let mut ij: libc::c_int = 0;
    let mut tp: libc::c_double = 0.;
    let mut beta: libc::c_double = 0.;
    let mut gamma_: libc::c_double = 0.;
    let mut alpha: libc::c_double = 0.;
    let mut delta: libc::c_double = 0.;
    w = w.offset(-1);
    z__ = z__.offset(-1);
    a = a.offset(-1);
    if !(*sigma == 0.0f64) {
        ij = 1 as libc::c_int;
        t = one / *sigma;
        if !(*sigma > 0.0f64) {
            i__1 = *n;
            i__ = 1 as libc::c_int;
            while i__ <= i__1 {
                *w.offset(i__ as isize) = *z__.offset(i__ as isize);
                i__ += 1;
            }
            i__1 = *n;
            i__ = 1 as libc::c_int;
            while i__ <= i__1 {
                v = *w.offset(i__ as isize);
                t += v * v / *a.offset(ij as isize);
                i__2 = *n;
                j = i__ + 1 as libc::c_int;
                while j <= i__2 {
                    ij += 1;
                    *w.offset(j as isize) -= v * *a.offset(ij as isize);
                    j += 1;
                }
                ij += 1;
                i__ += 1;
            }
            if t >= 0.0f64 {
                t = epmach / *sigma;
            }
            i__1 = *n;
            i__ = 1 as libc::c_int;
            while i__ <= i__1 {
                j = *n + 1 as libc::c_int - i__;
                ij -= i__;
                u = *w.offset(j as isize);
                *w.offset(j as isize) = t;
                t -= u * u / *a.offset(ij as isize);
                i__ += 1;
            }
        }
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            v = *z__.offset(i__ as isize);
            delta = v / *a.offset(ij as isize);
            if *sigma < 0.0f64 {
                tp = *w.offset(i__ as isize);
            } else {
                tp = t + delta * v;
            }
            alpha = tp / t;
            *a.offset(ij as isize) = alpha * *a.offset(ij as isize);
            if i__ == *n {
                break;
            }
            beta = delta / tp;
            if alpha > four {
                gamma_ = t / tp;
                i__2 = *n;
                j = i__ + 1 as libc::c_int;
                while j <= i__2 {
                    ij += 1;
                    u = *a.offset(ij as isize);
                    *a.offset(ij as isize) = gamma_ * u + beta * *z__.offset(j as isize);
                    *z__.offset(j as isize) -= v * u;
                    j += 1;
                }
            } else {
                i__2 = *n;
                j = i__ + 1 as libc::c_int;
                while j <= i__2 {
                    ij += 1;
                    *z__.offset(j as isize) -= v * *a.offset(ij as isize);
                    *a.offset(ij as isize) += beta * *z__.offset(j as isize);
                    j += 1;
                }
            }
            ij += 1;
            t = tp;
            i__ += 1;
        }
    }
}
unsafe fn slsqpb_(
    mut m: *mut libc::c_int,
    mut meq: *mut libc::c_int,
    mut la: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut x: *mut libc::c_double,
    mut xl: *const libc::c_double,
    mut xu: *const libc::c_double,
    mut f: *mut libc::c_double,
    mut c__: *mut libc::c_double,
    mut g: *mut libc::c_double,
    mut a: *mut libc::c_double,
    mut acc: *mut libc::c_double,
    mut iter: *mut libc::c_int,
    mut mode: *mut libc::c_int,
    mut r__: *mut libc::c_double,
    mut l: *mut libc::c_double,
    mut x0: *mut libc::c_double,
    mut mu: *mut libc::c_double,
    mut s: *mut libc::c_double,
    mut u: *mut libc::c_double,
    mut v: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut iw: *mut libc::c_int,
    mut state: *mut slsqpb_state,
) {
    let mut current_block: u64;
    let one: libc::c_double = 1.0f64;
    let alfmin: libc::c_double = 0.1f64;
    let hun: libc::c_double = 100.0f64;
    let ten: libc::c_double = 10.0f64;
    let two: libc::c_double = 2.0f64;
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut d__2: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    let mut t: libc::c_double = 0.;
    let mut f0: libc::c_double = 0.;
    let mut h1: libc::c_double = 0.;
    let mut h2: libc::c_double = 0.;
    let mut h3: libc::c_double = 0.;
    let mut h4: libc::c_double = 0.;
    let mut n1: libc::c_int = 0;
    let mut n2: libc::c_int = 0;
    let mut n3: libc::c_int = 0;
    let mut t0: libc::c_double = 0.;
    let mut gs: libc::c_double = 0.;
    let mut tol: libc::c_double = 0.;
    let mut line: libc::c_int = 0;
    let mut alpha: libc::c_double = 0.;
    let mut iexact: libc::c_int = 0;
    let mut incons: libc::c_int = 0;
    let mut ireset: libc::c_int = 0;
    let mut itermx: libc::c_int = 0;
    t = (*state).t;
    f0 = (*state).f0;
    h1 = (*state).h1;
    h2 = (*state).h2;
    h3 = (*state).h3;
    h4 = (*state).h4;
    n1 = (*state).n1;
    n2 = (*state).n2;
    n3 = (*state).n3;
    t0 = (*state).t0;
    gs = (*state).gs;
    tol = (*state).tol;
    line = (*state).line;
    alpha = (*state).alpha;
    iexact = (*state).iexact;
    incons = (*state).incons;
    ireset = (*state).ireset;
    itermx = (*state).itermx;
    mu = mu.offset(-1);
    c__ = c__.offset(-1);
    v = v.offset(-1);
    u = u.offset(-1);
    s = s.offset(-1);
    x0 = x0.offset(-1);
    l = l.offset(-1);
    r__ = r__.offset(-1);
    a_dim1 = *la;
    a_offset = 1 as libc::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    g = g.offset(-1);
    xu = xu.offset(-1);
    xl = xl.offset(-1);
    x = x.offset(-1);
    w = w.offset(-1);
    iw = iw.offset(-1);
    if *mode == -(1 as libc::c_int) {
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            *u.offset(i__ as isize) = *g.offset(i__ as isize)
                - ddot_sl__(
                    m,
                    &mut *a.offset((i__ * a_dim1 + 1 as libc::c_int) as isize),
                    1 as libc::c_int,
                    &mut *r__.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                )
                - *v.offset(i__ as isize);
            i__ += 1;
        }
        k = 0 as libc::c_int;
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            h1 = 0.0f64;
            k += 1;
            i__2 = *n;
            j = i__ + 1 as libc::c_int;
            while j <= i__2 {
                k += 1;
                h1 += *l.offset(k as isize) * *s.offset(j as isize);
                j += 1;
            }
            *v.offset(i__ as isize) = *s.offset(i__ as isize) + h1;
            i__ += 1;
        }
        k = 1 as libc::c_int;
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            *v.offset(i__ as isize) = *l.offset(k as isize) * *v.offset(i__ as isize);
            k = k + n1 - i__;
            i__ += 1;
        }
        i__ = *n;
        while i__ >= 1 as libc::c_int {
            h1 = 0.0f64;
            k = i__;
            i__1 = i__ - 1 as libc::c_int;
            j = 1 as libc::c_int;
            while j <= i__1 {
                h1 += *l.offset(k as isize) * *v.offset(j as isize);
                k = k + *n - j;
                j += 1;
            }
            *v.offset(i__ as isize) += h1;
            i__ -= 1;
        }
        h1 = ddot_sl__(
            n,
            &mut *s.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
            &mut *u.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        h2 = ddot_sl__(
            n,
            &mut *s.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
            &mut *v.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        h3 = h2 * 0.2f64;
        if h1 < h3 {
            h4 = (h2 - h3) / (h2 - h1);
            h1 = h3;
            dscal_sl__(
                n,
                &h4,
                &mut *u.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
            );
            d__1 = one - h4;
            daxpy_sl__(
                n,
                &d__1,
                v.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
                &mut *u.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
            );
        }
        d__1 = one / h1;
        ldl_(
            n,
            &mut *l.offset(1 as libc::c_int as isize),
            &mut *u.offset(1 as libc::c_int as isize),
            &mut d__1,
            &mut *v.offset(1 as libc::c_int as isize),
        );
        d__1 = -one / h2;
        ldl_(
            n,
            &mut *l.offset(1 as libc::c_int as isize),
            &mut *v.offset(1 as libc::c_int as isize),
            &mut d__1,
            &mut *u.offset(1 as libc::c_int as isize),
        );
        current_block = 6888897117580844842;
    } else if *mode == 0 as libc::c_int {
        itermx = *iter;
        if *acc >= 0.0f64 {
            iexact = 0 as libc::c_int;
        } else {
            iexact = 1 as libc::c_int;
        }
        *acc = (*acc).abs();
        tol = ten * *acc;
        *iter = 0 as libc::c_int;
        ireset = 0 as libc::c_int;
        n1 = *n + 1 as libc::c_int;
        n2 = n1 * *n / 2 as libc::c_int;
        n3 = n2 + 1 as libc::c_int;
        *s.offset(1 as libc::c_int as isize) = 0.0f64;
        *mu.offset(1 as libc::c_int as isize) = 0.0f64;
        dcopy___(
            n,
            s.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *s.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        dcopy___(
            m,
            mu.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *mu.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        current_block = 10732253873254351846;
    } else {
        t = *f;
        i__1 = *m;
        j = 1 as libc::c_int;
        while j <= i__1 {
            if j <= *meq {
                h1 = *c__.offset(j as isize);
            } else {
                h1 = 0.0f64;
            }
            d__1 = -*c__.offset(j as isize);
            t += *mu.offset(j as isize) * (if d__1 >= h1 { d__1 } else { h1 });
            j += 1;
        }
        h1 = t - t0;
        match iexact + 1 as libc::c_int {
            1 => {
                current_block = 11231496966810938197;
                match current_block {
                    11231496966810938197 => {
                        if nlopt_isfinite(h1) != 0 {
                            if h1 <= h3 / ten || line > 10 as libc::c_int {
                                current_block = 2367719738530698383;
                            } else {
                                d__1 = h3 / (two * (h3 - h1));
                                alpha = if d__1 >= alfmin { d__1 } else { alfmin };
                                current_block = 18018047807550034048;
                            }
                        } else {
                            alpha = if alpha * 0.5f64 >= alfmin {
                                alpha * 0.5f64
                            } else {
                                alfmin
                            };
                            current_block = 18018047807550034048;
                        }
                    }
                    _ => {}
                }
                match current_block {
                    18018047807550034048 => {}
                    _ => {
                        h3 = 0.0f64;
                        i__1 = *m;
                        j = 1 as libc::c_int;
                        while j <= i__1 {
                            if j <= *meq {
                                h1 = *c__.offset(j as isize);
                            } else {
                                h1 = 0.0f64;
                            }
                            d__1 = -*c__.offset(j as isize);
                            h3 += if d__1 >= h1 { d__1 } else { h1 };
                            j += 1;
                        }
                        d__1 = *f - f0;
                        if ((d__1).abs() < *acc
                            || dnrm2___(
                                n,
                                &mut *s.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            ) < *acc)
                            && h3 < *acc
                        {
                            *mode = 0 as libc::c_int;
                        } else {
                            *mode = -(1 as libc::c_int);
                        }
                        current_block = 10145886053577303292;
                    }
                }
            }
            2 => {
                current_block = 18080407663732281061;
            }
            _ => {
                current_block = 2367719738530698383;
                match current_block {
                    11231496966810938197 => {
                        if nlopt_isfinite(h1) != 0 {
                            if h1 <= h3 / ten || line > 10 as libc::c_int {
                                current_block = 2367719738530698383;
                            } else {
                                d__1 = h3 / (two * (h3 - h1));
                                alpha = if d__1 >= alfmin { d__1 } else { alfmin };
                                current_block = 18018047807550034048;
                            }
                        } else {
                            alpha = if alpha * 0.5f64 >= alfmin {
                                alpha * 0.5f64
                            } else {
                                alfmin
                            };
                            current_block = 18018047807550034048;
                        }
                    }
                    _ => {}
                }
                match current_block {
                    18018047807550034048 => {}
                    _ => {
                        h3 = 0.0f64;
                        i__1 = *m;
                        j = 1 as libc::c_int;
                        while j <= i__1 {
                            if j <= *meq {
                                h1 = *c__.offset(j as isize);
                            } else {
                                h1 = 0.0f64;
                            }
                            d__1 = -*c__.offset(j as isize);
                            h3 += if d__1 >= h1 { d__1 } else { h1 };
                            j += 1;
                        }
                        d__1 = *f - f0;
                        if ((d__1).abs() < *acc
                            || dnrm2___(
                                n,
                                &mut *s.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            ) < *acc)
                            && h3 < *acc
                        {
                            *mode = 0 as libc::c_int;
                        } else {
                            *mode = -(1 as libc::c_int);
                        }
                        current_block = 10145886053577303292;
                    }
                }
            }
        }
    }
    'c_6480: loop {
        match current_block {
            18080407663732281061 => {
                *mode = 9 as libc::c_int;
                return;
            }
            18018047807550034048 => {
                line += 1;
                h3 = alpha * h3;
                dscal_sl__(
                    n,
                    &alpha,
                    &mut *s.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                dcopy___(
                    n,
                    x0.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                daxpy_sl__(
                    n,
                    &one,
                    s.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                i__1 = *n;
                i__ = 1 as libc::c_int;
                while i__ <= i__1 {
                    if *x.offset(i__ as isize) < *xl.offset(i__ as isize) {
                        *x.offset(i__ as isize) = *xl.offset(i__ as isize);
                    } else if *x.offset(i__ as isize) > *xu.offset(i__ as isize) {
                        *x.offset(i__ as isize) = *xu.offset(i__ as isize);
                    }
                    i__ += 1;
                }
                *mode = if line == 1 as libc::c_int {
                    -(2 as libc::c_int)
                } else {
                    1 as libc::c_int
                };
                current_block = 10145886053577303292;
            }
            6888897117580844842 => {
                *iter += 1;
                *mode = 9 as libc::c_int;
                if *iter > itermx && itermx > 0 as libc::c_int {
                    current_block = 10145886053577303292;
                } else {
                    dcopy___(
                        n,
                        &*xl.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *u.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    dcopy___(
                        n,
                        &*xu.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *v.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    d__1 = -one;
                    daxpy_sl__(
                        n,
                        &d__1,
                        x.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *u.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    d__1 = -one;
                    daxpy_sl__(
                        n,
                        &d__1,
                        x.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *v.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    h4 = one;
                    lsq_(
                        m,
                        meq,
                        n,
                        &mut n3,
                        la,
                        &mut *l.offset(1 as libc::c_int as isize),
                        &mut *g.offset(1 as libc::c_int as isize),
                        &mut *a.offset(a_offset as isize),
                        &mut *c__.offset(1 as libc::c_int as isize),
                        u.offset(1 as libc::c_int as isize),
                        v.offset(1 as libc::c_int as isize),
                        &mut *s.offset(1 as libc::c_int as isize),
                        &mut *r__.offset(1 as libc::c_int as isize),
                        &mut *w.offset(1 as libc::c_int as isize),
                        &mut *iw.offset(1 as libc::c_int as isize),
                        mode,
                    );
                    if *mode == 6 as libc::c_int {
                        if *n == *meq {
                            *mode = 4 as libc::c_int;
                        }
                    }
                    if *mode == 4 as libc::c_int {
                        i__1 = *m;
                        j = 1 as libc::c_int;
                        while j <= i__1 {
                            if j <= *meq {
                                *a.offset((j + n1 * a_dim1) as isize) = -*c__.offset(j as isize);
                            } else {
                                d__1 = -*c__.offset(j as isize);
                                *a.offset((j + n1 * a_dim1) as isize) =
                                    if d__1 >= 0.0f64 { d__1 } else { 0.0f64 };
                            }
                            j += 1;
                        }
                        *s.offset(1 as libc::c_int as isize) = 0.0f64;
                        dcopy___(
                            n,
                            s.offset(1 as libc::c_int as isize),
                            0 as libc::c_int,
                            &mut *s.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        );
                        h3 = 0.0f64;
                        *g.offset(n1 as isize) = 0.0f64;
                        *l.offset(n3 as isize) = hun;
                        *s.offset(n1 as isize) = one;
                        *u.offset(n1 as isize) = 0.0f64;
                        *v.offset(n1 as isize) = one;
                        incons = 0 as libc::c_int;
                        loop {
                            lsq_(
                                m,
                                meq,
                                &mut n1,
                                &mut n3,
                                la,
                                &mut *l.offset(1 as libc::c_int as isize),
                                &mut *g.offset(1 as libc::c_int as isize),
                                &mut *a.offset(a_offset as isize),
                                &mut *c__.offset(1 as libc::c_int as isize),
                                u.offset(1 as libc::c_int as isize),
                                v.offset(1 as libc::c_int as isize),
                                &mut *s.offset(1 as libc::c_int as isize),
                                &mut *r__.offset(1 as libc::c_int as isize),
                                &mut *w.offset(1 as libc::c_int as isize),
                                &mut *iw.offset(1 as libc::c_int as isize),
                                mode,
                            );
                            h4 = one - *s.offset(n1 as isize);
                            if !(*mode == 4 as libc::c_int) {
                                break;
                            }
                            *l.offset(n3 as isize) = ten * *l.offset(n3 as isize);
                            incons += 1;
                            if incons > 5 as libc::c_int {
                                current_block = 10145886053577303292;
                                continue 'c_6480;
                            }
                        }
                        if *mode != 1 as libc::c_int {
                            current_block = 10145886053577303292;
                            continue;
                        }
                    } else if *mode != 1 as libc::c_int {
                        current_block = 10145886053577303292;
                        continue;
                    }
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        *v.offset(i__ as isize) = *g.offset(i__ as isize)
                            - ddot_sl__(
                                m,
                                &mut *a.offset((i__ * a_dim1 + 1 as libc::c_int) as isize),
                                1 as libc::c_int,
                                &mut *r__.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            );
                        i__ += 1;
                    }
                    f0 = *f;
                    dcopy___(
                        n,
                        x.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *x0.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    gs = ddot_sl__(
                        n,
                        &mut *g.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *s.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    h1 = (gs).abs();
                    h2 = 0.0f64;
                    i__1 = *m;
                    j = 1 as libc::c_int;
                    while j <= i__1 {
                        if j <= *meq {
                            h3 = *c__.offset(j as isize);
                        } else {
                            h3 = 0.0f64;
                        }
                        d__1 = -*c__.offset(j as isize);
                        h2 += if d__1 >= h3 { d__1 } else { h3 };
                        d__1 = *r__.offset(j as isize);
                        h3 = (d__1).abs();
                        d__1 = h3;
                        d__2 = (*mu.offset(j as isize) + h3) / two;
                        *mu.offset(j as isize) = if d__1 >= d__2 { d__1 } else { d__2 };
                        d__1 = *c__.offset(j as isize);
                        h1 += h3 * (d__1).abs();
                        j += 1;
                    }
                    *mode = 0 as libc::c_int;
                    if h1 < *acc && h2 < *acc {
                        current_block = 10145886053577303292;
                        continue;
                    }
                    h1 = 0.0f64;
                    i__1 = *m;
                    j = 1 as libc::c_int;
                    while j <= i__1 {
                        if j <= *meq {
                            h3 = *c__.offset(j as isize);
                        } else {
                            h3 = 0.0f64;
                        }
                        d__1 = -*c__.offset(j as isize);
                        h1 += *mu.offset(j as isize) * (if d__1 >= h3 { d__1 } else { h3 });
                        j += 1;
                    }
                    t0 = *f + h1;
                    h3 = gs - h1 * h4;
                    *mode = 8 as libc::c_int;
                    if h3 >= 0.0f64 {
                        current_block = 10732253873254351846;
                        continue;
                    }
                    line = 0 as libc::c_int;
                    alpha = one;
                    if iexact == 1 as libc::c_int {
                        current_block = 18080407663732281061;
                    } else {
                        current_block = 18018047807550034048;
                    }
                }
            }
            10732253873254351846 => {
                ireset += 1;
                if ireset > 5 as libc::c_int {
                    d__1 = *f - f0;
                    if ((d__1).abs() < tol
                        || dnrm2___(
                            n,
                            &mut *s.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        ) < tol)
                        && h3 < tol
                    {
                        *mode = 0 as libc::c_int;
                    } else {
                        *mode = 8 as libc::c_int;
                    }
                    current_block = 10145886053577303292;
                } else {
                    *l.offset(1 as libc::c_int as isize) = 0.0f64;
                    dcopy___(
                        &mut n2,
                        l.offset(1 as libc::c_int as isize),
                        0 as libc::c_int,
                        &mut *l.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    j = 1 as libc::c_int;
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        *l.offset(j as isize) = one;
                        j = j + n1 - i__;
                        i__ += 1;
                    }
                    current_block = 6888897117580844842;
                }
            }
            _ => {
                (*state).t = t;
                (*state).f0 = f0;
                (*state).h1 = h1;
                (*state).h2 = h2;
                (*state).h3 = h3;
                (*state).h4 = h4;
                (*state).n1 = n1;
                (*state).n2 = n2;
                (*state).n3 = n3;
                (*state).t0 = t0;
                (*state).gs = gs;
                (*state).tol = tol;
                (*state).line = line;
                (*state).alpha = alpha;
                (*state).iexact = iexact;
                (*state).incons = incons;
                (*state).ireset = ireset;
                (*state).itermx = itermx;
                return;
            }
        }
    }
}
unsafe fn slsqp(
    mut m: *mut libc::c_int,
    mut meq: *mut libc::c_int,
    mut la: *mut libc::c_int,
    mut n: *mut libc::c_int,
    mut x: *mut libc::c_double,
    mut xl: *const libc::c_double,
    mut xu: *const libc::c_double,
    mut f: *mut libc::c_double,
    mut c__: *mut libc::c_double,
    mut g: *mut libc::c_double,
    mut a: *mut libc::c_double,
    mut acc: *mut libc::c_double,
    mut iter: *mut libc::c_int,
    mut mode: *mut libc::c_int,
    mut w: *mut libc::c_double,
    mut l_w__: *mut libc::c_int,
    mut jw: *mut libc::c_int,
    mut l_jw__: *mut libc::c_int,
    mut state: *mut slsqpb_state,
) {
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut n1: libc::c_int = 0;
    let mut il: libc::c_int = 0;
    let mut im: libc::c_int = 0;
    let mut ir: libc::c_int = 0;
    let mut is: libc::c_int = 0;
    let mut iu: libc::c_int = 0;
    let mut iv: libc::c_int = 0;
    let mut iw: libc::c_int = 0;
    let mut ix: libc::c_int = 0;
    let mut mineq: libc::c_int = 0;
    c__ = c__.offset(-1);
    a_dim1 = *la;
    a_offset = 1 as libc::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    g = g.offset(-1);
    xu = xu.offset(-1);
    xl = xl.offset(-1);
    x = x.offset(-1);
    w = w.offset(-1);
    jw = jw.offset(-1);
    n1 = *n + 1 as libc::c_int;
    mineq = *m - *meq + n1 + n1;
    il = (n1 * 3 as libc::c_int + *m) * (n1 + 1 as libc::c_int)
        + (n1 - *meq + 1 as libc::c_int) * (mineq + 2 as libc::c_int)
        + (mineq << 1 as libc::c_int)
        + (n1 + mineq) * (n1 - *meq)
        + (*meq << 1 as libc::c_int)
        + n1 * *n / 2 as libc::c_int
        + (*m << 1 as libc::c_int)
        + *n * 3 as libc::c_int
        + (n1 << 2 as libc::c_int)
        + 1 as libc::c_int;
    i__1 = mineq;
    i__2 = n1 - *meq;
    im = if i__1 >= i__2 { i__1 } else { i__2 };
    if *l_w__ < il || *l_jw__ < im {
        *mode = (if 10 as libc::c_int >= il {
            10 as libc::c_int
        } else {
            il
        }) * 1000 as libc::c_int;
        *mode += if 10 as libc::c_int >= im {
            10 as libc::c_int
        } else {
            im
        };
        return;
    }
    im = 1 as libc::c_int;
    il = im
        + (if 1 as libc::c_int >= *m {
            1 as libc::c_int
        } else {
            *m
        });
    il = im + *la;
    ix = il + n1 * *n / 2 as libc::c_int + 1 as libc::c_int;
    ir = ix + *n;
    is = ir
        + *n
        + *n
        + (if 1 as libc::c_int >= *m {
            1 as libc::c_int
        } else {
            *m
        });
    is = ir + *n + *n + *la;
    iu = is + n1;
    iv = iu + n1;
    iw = iv + n1;
    slsqpb_(
        m,
        meq,
        la,
        n,
        &mut *x.offset(1 as libc::c_int as isize),
        &*xl.offset(1 as libc::c_int as isize),
        &*xu.offset(1 as libc::c_int as isize),
        f,
        &mut *c__.offset(1 as libc::c_int as isize),
        &mut *g.offset(1 as libc::c_int as isize),
        &mut *a.offset(a_offset as isize),
        acc,
        iter,
        mode,
        &mut *w.offset(ir as isize),
        &mut *w.offset(il as isize),
        &mut *w.offset(ix as isize),
        &mut *w.offset(im as isize),
        &mut *w.offset(is as isize),
        &mut *w.offset(iu as isize),
        &mut *w.offset(iv as isize),
        &mut *w.offset(iw as isize),
        &mut *jw.offset(1 as libc::c_int as isize),
        state,
    );
    let ref mut fresh1 = (*state).x0;
    *fresh1 = &mut *w.offset(ix as isize) as *mut libc::c_double;
}
unsafe fn length_work(
    mut LEN_W: *mut libc::c_int,
    mut LEN_JW: *mut libc::c_int,
    mut M: libc::c_int,
    mut MEQ: libc::c_int,
    mut N: libc::c_int,
) {
    let mut N1: libc::c_int = N + 1 as libc::c_int;
    let mut MINEQ: libc::c_int = M - MEQ + N1 + N1;
    *LEN_W = (3 as libc::c_int * N1 + M) * (N1 + 1 as libc::c_int)
        + (N1 - MEQ + 1 as libc::c_int) * (MINEQ + 2 as libc::c_int)
        + 2 as libc::c_int * MINEQ
        + (N1 + MINEQ) * (N1 - MEQ)
        + 2 as libc::c_int * MEQ
        + N1
        + (N + 1 as libc::c_int) * N / 2 as libc::c_int
        + 2 as libc::c_int * M
        + 3 as libc::c_int * N
        + 3 as libc::c_int * N1
        + 1 as libc::c_int;
    *LEN_JW = MINEQ;
}

pub unsafe fn nlopt_slsqp<U>(
    mut n: libc::c_uint,
    mut f: nlopt_func,
    mut f_data: *mut libc::c_void,
    mut m: libc::c_uint,
    mut fc: *mut nlopt_constraint,
    mut p: libc::c_uint,
    mut h: *mut nlopt_constraint,
    mut lb: *const libc::c_double,
    mut ub: *const libc::c_double,
    mut x: *mut libc::c_double,
    mut minf: *mut libc::c_double,
    mut stop: *mut nlopt_stopping,
) -> nlopt_result {
    let mut current_block: u64;
    let mut state: slsqpb_state = slsqpb_state {
        t: 0 as libc::c_int as libc::c_double,
        f0: 0 as libc::c_int as libc::c_double,
        h1: 0 as libc::c_int as libc::c_double,
        h2: 0 as libc::c_int as libc::c_double,
        h3: 0 as libc::c_int as libc::c_double,
        h4: 0 as libc::c_int as libc::c_double,
        n1: 0 as libc::c_int,
        n2: 0 as libc::c_int,
        n3: 0 as libc::c_int,
        t0: 0 as libc::c_int as libc::c_double,
        gs: 0 as libc::c_int as libc::c_double,
        tol: 0 as libc::c_int as libc::c_double,
        line: 0 as libc::c_int,
        alpha: 0 as libc::c_int as libc::c_double,
        iexact: 0 as libc::c_int,
        incons: 0 as libc::c_int,
        ireset: 0 as libc::c_int,
        itermx: 0 as libc::c_int,
        x0: 0 as *mut libc::c_double,
    };

    let mut mtot: libc::c_uint = nlopt_count_constraints(m, fc);
    let mut ptot: libc::c_uint = nlopt_count_constraints(p, h);
    let mut work: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut cgrad: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut c: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut grad: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut w: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut fcur: libc::c_double = 0.;
    let mut xcur: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut fprev: libc::c_double = 0.;
    let mut xprev: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut cgradtmp: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut mpi: libc::c_int = mtot.wrapping_add(ptot) as libc::c_int;
    let mut pi: libc::c_int = ptot as libc::c_int;
    let mut ni: libc::c_int = n as libc::c_int;
    let mut mpi1: libc::c_int = if mpi > 0 as libc::c_int {
        mpi
    } else {
        1 as libc::c_int
    };
    let mut len_w: libc::c_int = 0;
    let mut len_jw: libc::c_int = 0;
    let mut jw: *mut libc::c_int = 0 as *mut libc::c_int;
    let mut mode: libc::c_int = 0 as libc::c_int;
    let mut prev_mode: libc::c_int = 0 as libc::c_int;
    let mut acc: libc::c_double = 0 as libc::c_int as libc::c_double;
    let mut iter: libc::c_int = 0 as libc::c_int;
    let mut i: libc::c_uint = 0;
    let mut ii: libc::c_uint = 0;
    let mut ret: nlopt_result = NLOPT_SUCCESS;
    let mut feasible: libc::c_int = 0;
    let mut feasible_cur: libc::c_int = 0;
    let mut infeasibility: libc::c_double = ::std::f64::INFINITY;
    let mut infeasibility_cur: libc::c_double = ::std::f64::INFINITY;
    let mut max_cdim: libc::c_uint = 0;
    let mut want_grad: libc::c_int = 1 as libc::c_int;
    max_cdim = if nlopt_max_constraint_dim(m, fc) >= nlopt_max_constraint_dim(p, h) {
        nlopt_max_constraint_dim(m, fc)
    } else {
        nlopt_max_constraint_dim(p, h)
    };
    length_work(&mut len_w, &mut len_jw, mpi, pi, ni);
    // work = malloc(
    //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
    //         .wrapping_mul(
    //             (mpi1 as libc::c_uint)
    //                 .wrapping_mul(n.wrapping_add(1 as libc::c_int as libc::c_uint))
    //                 .wrapping_add(mpi as libc::c_uint)
    //                 .wrapping_add(n)
    //                 .wrapping_add(1 as libc::c_int as libc::c_uint)
    //                 .wrapping_add(n)
    //                 .wrapping_add(n)
    //                 .wrapping_add(max_cdim.wrapping_mul(n))
    //                 .wrapping_add(len_w as libc::c_uint) as libc::c_ulong,
    //         )
    //         .wrapping_add(
    //             (::std::mem::size_of::<libc::c_int>() as libc::c_ulong)
    //                 .wrapping_mul(len_jw as libc::c_uint as libc::c_ulong),
    //         ),
    // ) as *mut libc::c_double;

    let space_size = (mpi1 as libc::c_uint)
        .wrapping_mul(n.wrapping_add(1 as libc::c_int as libc::c_uint))
        .wrapping_add(mpi as libc::c_uint)
        .wrapping_add(n)
        .wrapping_add(1 as libc::c_int as libc::c_uint)
        .wrapping_add(n)
        .wrapping_add(n)
        .wrapping_add(max_cdim.wrapping_mul(n))
        .wrapping_add(len_w as libc::c_uint);
    let space_size = space_size.wrapping_add(
        (::std::mem::size_of::<libc::c_int>() as libc::c_ulong)
            .wrapping_mul(len_jw as libc::c_uint as libc::c_ulong),
    );
    let mut space: Box<Vec<libc::c_double>> =
        Box::new(vec![0.; usize::try_from(space_size).unwrap()]);
    work = space.as_mut_ptr() as *mut libc::c_double;
    std::mem::forget(space);

    if work.is_null() {
        return NLOPT_OUT_OF_MEMORY;
    }
    cgrad = work;
    c = cgrad.offset(
        (mpi1 as libc::c_uint).wrapping_mul(n.wrapping_add(1 as libc::c_int as libc::c_uint))
            as isize,
    );
    grad = c.offset(mpi as isize);
    xcur = grad.offset(n as isize).offset(1 as libc::c_int as isize);
    xprev = xcur.offset(n as isize);
    cgradtmp = xprev.offset(n as isize);
    w = cgradtmp.offset(max_cdim.wrapping_mul(n) as isize);
    jw = w.offset(len_w as isize) as *mut libc::c_int;
    memcpy(
        xcur as *mut libc::c_void,
        x as *const libc::c_void,
        (::std::mem::size_of::<libc::c_double>() as libc::c_ulong).wrapping_mul(n as libc::c_ulong),
    );
    memcpy(
        xprev as *mut libc::c_void,
        x as *const libc::c_void,
        (::std::mem::size_of::<libc::c_double>() as libc::c_ulong).wrapping_mul(n as libc::c_ulong),
    );
    *minf = ::std::f64::INFINITY;
    fcur = *minf;
    fprev = fcur;
    feasible_cur = 0 as libc::c_int;
    feasible = feasible_cur;
    'c_6042: loop {
        want_grad = 1 as libc::c_int;
        's_146: loop {
            let mut newgrad: *mut libc::c_double = 0 as *mut libc::c_double;
            let mut newcgrad: *mut libc::c_double = 0 as *mut libc::c_double;
            if want_grad != 0 {
                newgrad = grad;
                newcgrad = cgradtmp;
            }
            feasible_cur = 1 as libc::c_int;
            infeasibility_cur = 0 as libc::c_int as libc::c_double;
            fcur = f.expect("non-null function pointer")(n, xcur, newgrad, f_data);
            let ref mut fresh2 = *(*stop).nevals_p;
            *fresh2 += 1;
            if nlopt_stop_forced(stop) != 0 {
                fcur = ::std::f64::INFINITY;
                ret = NLOPT_FORCED_STOP;
                break 'c_6042;
            } else {
                if nlopt_isfinite(fcur) != 0 {
                    want_grad = 0 as libc::c_int;
                    ii = 0 as libc::c_int as libc::c_uint;
                    i = 0 as libc::c_int as libc::c_uint;
                    while i < p {
                        let mut j: libc::c_uint = 0;
                        let mut k: libc::c_uint = 0;
                        nlopt_eval_constraint::<U>(
                            c.offset(ii as isize),
                            newcgrad,
                            h.offset(i as isize),
                            n,
                            xcur,
                        );
                        if nlopt_stop_forced(stop) != 0 {
                            ret = NLOPT_FORCED_STOP;
                            break 'c_6042;
                        } else {
                            k = 0 as libc::c_int as libc::c_uint;
                            while k < (*h.offset(i as isize)).m {
                                infeasibility_cur =
                                    if infeasibility_cur >= (*c.offset(ii as isize)).abs() {
                                        infeasibility_cur
                                    } else {
                                        (*c.offset(ii as isize)).abs()
                                    };
                                feasible_cur = (feasible_cur != 0
                                    && (*c.offset(ii as isize)).abs()
                                        <= *((*h.offset(i as isize)).tol).offset(k as isize))
                                    as libc::c_int;
                                if !newcgrad.is_null() {
                                    j = 0 as libc::c_int as libc::c_uint;
                                    while j < n {
                                        *cgrad.offset(
                                            j.wrapping_mul(mpi1 as libc::c_uint).wrapping_add(ii)
                                                as isize,
                                        ) = *cgradtmp
                                            .offset(k.wrapping_mul(n).wrapping_add(j) as isize);
                                        j = j.wrapping_add(1);
                                    }
                                }
                                k = k.wrapping_add(1);
                                ii = ii.wrapping_add(1);
                            }
                            i = i.wrapping_add(1);
                        }
                    }
                    i = 0 as libc::c_int as libc::c_uint;
                    while i < m {
                        let mut j_0: libc::c_uint = 0;
                        let mut k_0: libc::c_uint = 0;
                        nlopt_eval_constraint::<U>(
                            c.offset(ii as isize),
                            newcgrad,
                            fc.offset(i as isize),
                            n,
                            xcur,
                        );
                        if nlopt_stop_forced(stop) != 0 {
                            ret = NLOPT_FORCED_STOP;
                            break 'c_6042;
                        } else {
                            k_0 = 0 as libc::c_int as libc::c_uint;
                            while k_0 < (*fc.offset(i as isize)).m {
                                infeasibility_cur = if infeasibility_cur >= *c.offset(ii as isize) {
                                    infeasibility_cur
                                } else {
                                    *c.offset(ii as isize)
                                };
                                feasible_cur = (feasible_cur != 0
                                    && *c.offset(ii as isize)
                                        <= *((*fc.offset(i as isize)).tol).offset(k_0 as isize))
                                    as libc::c_int;
                                if !newcgrad.is_null() {
                                    j_0 = 0 as libc::c_int as libc::c_uint;
                                    while j_0 < n {
                                        *cgrad.offset(
                                            j_0.wrapping_mul(mpi1 as libc::c_uint).wrapping_add(ii)
                                                as isize,
                                        ) = -*cgradtmp
                                            .offset(k_0.wrapping_mul(n).wrapping_add(j_0) as isize);
                                        j_0 = j_0.wrapping_add(1);
                                    }
                                }
                                *c.offset(ii as isize) = -*c.offset(ii as isize);
                                k_0 = k_0.wrapping_add(1);
                                ii = ii.wrapping_add(1);
                            }
                            i = i.wrapping_add(1);
                        }
                    }
                }
                loop {
                    prev_mode = mode;
                    if nlopt_isfinite(fcur) != 0
                        && (fcur < *minf && (feasible_cur != 0 || feasible == 0)
                            || feasible == 0 && infeasibility_cur < infeasibility)
                    {
                        *minf = fcur;
                        feasible = feasible_cur;
                        infeasibility = infeasibility_cur;
                        memcpy(
                            x as *mut libc::c_void,
                            xcur as *const libc::c_void,
                            (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                                .wrapping_mul(n as libc::c_ulong),
                        );
                    }
                    if mode == -(1 as libc::c_int) {
                        if nlopt_isinf(fprev) == 0 && feasible != 0 {
                            if nlopt_stop_ftol(stop, fcur, fprev) != 0 {
                                ret = NLOPT_FTOL_REACHED;
                            } else if nlopt_stop_x(stop, xcur, xprev) != 0 {
                                ret = NLOPT_XTOL_REACHED;
                            }
                        }
                        fprev = fcur;
                        memcpy(
                            xprev as *mut libc::c_void,
                            xcur as *const libc::c_void,
                            (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                                .wrapping_mul(n as libc::c_ulong),
                        );
                    }
                    if nlopt_stop_evals(stop) != 0 {
                        ret = NLOPT_MAXEVAL_REACHED;
                    } else if nlopt_stop_time(stop) != 0 {
                        ret = NLOPT_MAXTIME_REACHED;
                    } else if feasible != 0 && *minf < (*stop).minf_max {
                        ret = NLOPT_STOPVAL_REACHED;
                    }
                    if !(ret as libc::c_int == NLOPT_SUCCESS as libc::c_int) {
                        break 'c_6042;
                    }
                    slsqp(
                        &mut mpi,
                        &mut pi,
                        &mut mpi1,
                        &mut ni,
                        xcur,
                        lb,
                        ub,
                        &mut fcur,
                        c,
                        grad,
                        cgrad,
                        &mut acc,
                        &mut iter,
                        &mut mode,
                        w,
                        &mut len_w,
                        jw,
                        &mut len_jw,
                        &mut state,
                    );
                    match mode {
                        -1 => {
                            if !(prev_mode == -(2 as libc::c_int) && want_grad == 0) {
                                continue 'c_6042;
                            }
                        }
                        -2 => {
                            continue 'c_6042;
                        }
                        1 => {
                            break;
                        }
                        0 => {
                            break 'c_6042;
                        }
                        8 => {
                            ret = NLOPT_ROUNDOFF_LIMITED;
                            if feasible_cur != 0 {
                                let mut save_ftol_rel: libc::c_double = (*stop).ftol_rel;
                                let mut save_xtol_rel: libc::c_double = (*stop).xtol_rel;
                                let mut save_ftol_abs: libc::c_double = (*stop).ftol_abs;
                                (*stop).ftol_rel *= 10 as libc::c_int as libc::c_double;
                                (*stop).ftol_abs *= 10 as libc::c_int as libc::c_double;
                                (*stop).xtol_rel *= 10 as libc::c_int as libc::c_double;
                                if nlopt_stop_ftol(stop, fcur, state.f0) != 0 {
                                    ret = NLOPT_FTOL_REACHED;
                                } else if nlopt_stop_x(stop, xcur, state.x0) != 0 {
                                    ret = NLOPT_XTOL_REACHED;
                                }
                                (*stop).ftol_rel = save_ftol_rel;
                                (*stop).ftol_abs = save_ftol_abs;
                                (*stop).xtol_rel = save_xtol_rel;
                            }
                            break 'c_6042;
                        }
                        5 => {
                            current_block = 9714366428826605325;
                            break 's_146;
                        }
                        6 | 7 => {
                            current_block = 9714366428826605325;
                            break 's_146;
                        }
                        4 => {
                            current_block = 8507021683849773288;
                            break 's_146;
                        }
                        3 | 9 => {
                            current_block = 8507021683849773288;
                            break 's_146;
                        }
                        // 2 |
                        _ => {
                            nlopt_stop_msg(stop, "bug: workspace is too small");
                            ret = NLOPT_INVALID_ARGS;
                            break 'c_6042;
                        }
                    }
                }
            }
        }
        match current_block {
            9714366428826605325 => {
                ret = NLOPT_ROUNDOFF_LIMITED;
                break;
            }
            _ => {
                nlopt_stop_msg(stop, "bug: more than iter SQP iterations");
                ret = NLOPT_FAILURE;
                break;
            }
        }
    }
    if nlopt_isinf(*minf) != 0 {
        if nlopt_isinf(fcur) != 0 {
            *minf = fprev;
            memcpy(
                x as *mut libc::c_void,
                xprev as *const libc::c_void,
                (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                    .wrapping_mul(n as libc::c_ulong),
            );
        } else {
            *minf = fcur;
            memcpy(
                x as *mut libc::c_void,
                xcur as *const libc::c_void,
                (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                    .wrapping_mul(n as libc::c_ulong),
            );
        }
    }
    // free(work as *mut libc::c_void);
    let _ = Box::from_raw(work);
    return ret;
}
