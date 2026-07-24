#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut,
    static_mut_refs,
    unsafe_op_in_unsafe_fn,
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

pub(crate) fn nlopt_function_raw_callback<F: Func<T>, T>(
    n: ::core::ffi::c_uint,
    x: *const f64,
    g: *mut f64,
    params: *mut ::core::ffi::c_void,
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
    res
}

pub(crate) fn nlopt_constraint_raw_callback<F: Func<T>, T>(
    n: ::core::ffi::c_uint,
    x: *const f64,
    g: *mut f64,
    params: *mut ::core::ffi::c_void,
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
pub(crate) struct NLoptFunctionCfg<F: Func<T>, T> {
    pub objective_fn: F,
    pub user_data: T,
}

pub(crate) struct NLoptConstraintCfg<F: Func<T>, T> {
    pub constraint_fn: F,
    pub user_data: T,
}

/// A trait representing a function. It is used for objective and constraint functions.
///
/// A function takes the form of a closure `f(x: &[f64], gradient: Option<&mut [f64], user_data: &mut U) -> f64`
///
/// * `x` - `n`-dimensional array
/// * `gradient` - `n`-dimensional array to store the gradient `grad f(x)`. If `gradient` matches
///   `Some(x)`, the user is required to provide a gradient, otherwise the optimization will fail.
/// * `user_data` - user defined data used to make objective and constraints computations.
pub trait Func<U>: Fn(&[f64], Option<&mut [f64]>, &mut U) -> f64 {}
impl<T, U> Func<U> for T where T: Fn(&[f64], Option<&mut [f64]>, &mut U) -> f64 {}

// #![register_tool(c2rust)]
// #![feature(c_variadic, register_tool)]
//extern "C" {
// fn sqrt(_: ::core::ffi::c_double) -> ::core::ffi::c_double;
// fn fabs(_: ::core::ffi::c_double) -> ::core::ffi::c_double;
// fn malloc(_: ::core::ffi::c_ulong) -> *mut ::core::ffi::c_void;
// fn realloc(_: *mut ::core::ffi::c_void, _: ::core::ffi::c_ulong) -> *mut ::core::ffi::c_void;
// fn free(__ptr: *mut ::core::ffi::c_void);
// fn abort() -> !;
// fn memcpy(_: *mut ::core::ffi::c_void, _: *const ::core::ffi::c_void, _: ::core::ffi::c_ulong) -> *mut ::core::ffi::c_void;
// fn strlen(_: *const ::core::ffi::c_char) -> ::core::ffi::c_ulong;
// fn gettimeofday(__tv: *mut timeval, __tz: *mut ::core::ffi::c_void) -> ::core::ffi::c_int;
// fn vsnprintf(
//     _: *mut ::core::ffi::c_char,
//     _: ::core::ffi::c_ulong,
//     _: *const ::core::ffi::c_char,
//     _: ::std::ffi::VaList,
// ) -> ::core::ffi::c_int;
//}

unsafe fn memcpy(dst: *mut ::core::ffi::c_void, src: *const ::core::ffi::c_void, n: ::core::ffi::c_ulong) { unsafe {
    std::ptr::copy_nonoverlapping(src, dst, n as usize);
}}

// pub type __builtin_va_list = [__va_list_tag; 1];
// #[derive(Copy, Clone)]
// #[repr(C)]
// pub struct __va_list_tag {
//     pub gp_offset: ::core::ffi::c_uint,
//     pub fp_offset: ::core::ffi::c_uint,
//     pub overflow_arg_area: *mut ::core::ffi::c_void,
//     pub reg_save_area: *mut ::core::ffi::c_void,
// }
type __time_t = ::core::ffi::c_long;
type __suseconds_t = ::core::ffi::c_long;
type size_t = ::core::ffi::c_ulong;
#[derive(Copy, Clone)]
#[repr(C)]
struct timeval {
    pub tv_sec: __time_t,
    pub tv_usec: __suseconds_t,
}
type nlopt_func = Option<
    unsafe fn(
        ::core::ffi::c_uint,
        *const ::core::ffi::c_double,
        *mut ::core::ffi::c_double,
        *mut ::core::ffi::c_void,
    ) -> ::core::ffi::c_double,
>;
type nlopt_mfunc = Option<
    unsafe fn(
        ::core::ffi::c_uint,
        *mut ::core::ffi::c_double,
        ::core::ffi::c_uint,
        *const ::core::ffi::c_double,
        *mut ::core::ffi::c_double,
        *mut ::core::ffi::c_void,
    ) -> (),
>;
type nlopt_precond = Option<
    unsafe fn(
        ::core::ffi::c_uint,
        *const ::core::ffi::c_double,
        *const ::core::ffi::c_double,
        *mut ::core::ffi::c_double,
        *mut ::core::ffi::c_void,
    ) -> (),
>;
type nlopt_result = ::core::ffi::c_int;
const NLOPT_NUM_RESULTS: nlopt_result = 7;
const NLOPT_MAXTIME_REACHED: nlopt_result = 6;
const NLOPT_MAXEVAL_REACHED: nlopt_result = 5;
const NLOPT_XTOL_REACHED: nlopt_result = 4;
const NLOPT_FTOL_REACHED: nlopt_result = 3;
const NLOPT_STOPVAL_REACHED: nlopt_result = 2;
const NLOPT_SUCCESS: nlopt_result = 1;
const NLOPT_NUM_FAILURES: nlopt_result = -6;
const NLOPT_FORCED_STOP: nlopt_result = -5;
const NLOPT_ROUNDOFF_LIMITED: nlopt_result = -4;
const NLOPT_OUT_OF_MEMORY: nlopt_result = -3;
const NLOPT_INVALID_ARGS: nlopt_result = -2;
const NLOPT_FAILURE: nlopt_result = -1;
// pub type va_list = __builtin_va_list;
// #[derive(Copy, Clone)]
#[derive(Clone)]
#[repr(C)]
pub(crate) struct nlopt_stopping {
    pub n: ::core::ffi::c_uint,
    pub minf_max: ::core::ffi::c_double,
    pub ftol_rel: ::core::ffi::c_double,
    pub ftol_abs: ::core::ffi::c_double,
    pub xtol_rel: ::core::ffi::c_double,
    pub xtol_abs: *const ::core::ffi::c_double,
    pub x_weights: *const ::core::ffi::c_double,
    pub nevals_p: *mut ::core::ffi::c_int,
    pub maxeval: ::core::ffi::c_int,
    pub maxtime: ::core::ffi::c_double,
    pub start: ::core::ffi::c_double,
    pub force_stop: *mut ::core::ffi::c_int,
    pub stop_msg: String,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub(crate) struct nlopt_constraint {
    pub m: ::core::ffi::c_uint,
    pub f: nlopt_func,
    pub mf: nlopt_mfunc,
    pub pre: nlopt_precond,
    pub f_data: *mut ::core::ffi::c_void,
    pub tol: *mut ::core::ffi::c_double,
}
#[derive(Copy, Clone)]
#[repr(C)]
struct slsqpb_state {
    pub t: ::core::ffi::c_double,
    pub f0: ::core::ffi::c_double,
    pub h1: ::core::ffi::c_double,
    pub h2: ::core::ffi::c_double,
    pub h3: ::core::ffi::c_double,
    pub h4: ::core::ffi::c_double,
    pub n1: ::core::ffi::c_int,
    pub n2: ::core::ffi::c_int,
    pub n3: ::core::ffi::c_int,
    pub t0: ::core::ffi::c_double,
    pub gs: ::core::ffi::c_double,
    pub tol: ::core::ffi::c_double,
    pub line: ::core::ffi::c_int,
    pub alpha: ::core::ffi::c_double,
    pub iexact: ::core::ffi::c_int,
    pub incons: ::core::ffi::c_int,
    pub ireset: ::core::ffi::c_int,
    pub itermx: ::core::ffi::c_int,
    pub x0: *mut ::core::ffi::c_double,
}

unsafe fn nlopt_time_seed() -> ::core::ffi::c_ulong {
    // let mut tv = ::core::ffi::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // ::core::ffi::gettimeofday(&mut tv, 0 as *mut ::core::ffi::timezone);
    //return (tv.tv_sec ^ tv.tv_usec) as ::core::ffi::c_ulong;
    let start = SystemTime::now();
    let since_the_epoch = start.duration_since(UNIX_EPOCH).expect("Time flies");
    since_the_epoch.as_millis() as ::core::ffi::c_ulong
}

unsafe fn nlopt_seconds() -> ::core::ffi::c_double { unsafe {
    // static mut start_inited: ::core::ffi::c_int = 0 as ::core::ffi::c_int;
    // static mut start: ::core::ffi::timeval = ::core::ffi::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // let mut tv: ::core::ffi::timeval = ::core::ffi::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // if start_inited == 0 {
    //     start_inited = 1 as ::core::ffi::c_int;
    //     ::core::ffi::gettimeofday(&mut start, 0 as *mut ::core::ffi::timezone);
    // }
    // ::core::ffi::gettimeofday(&mut tv, 0 as *mut ::core::ffi::timezone);
    // return (tv.tv_sec - start.tv_sec) as ::core::ffi::c_double
    //     + 1.0e-6f64 * (tv.tv_usec - start.tv_usec) as ::core::ffi::c_double;
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
}}
unsafe fn sc(
    mut x: ::core::ffi::c_double,
    mut smin: ::core::ffi::c_double,
    mut smax: ::core::ffi::c_double,
) -> ::core::ffi::c_double {
    return smin + x * (smax - smin);
}
unsafe fn vector_norm(
    mut n: ::core::ffi::c_uint,
    mut vec: *const ::core::ffi::c_double,
    mut w: *const ::core::ffi::c_double,
    mut scale_min: *const ::core::ffi::c_double,
    mut scale_max: *const ::core::ffi::c_double,
) -> ::core::ffi::c_double {
    let mut i: ::core::ffi::c_uint = 0;
    let mut ret: ::core::ffi::c_double = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    if !scale_min.is_null() && !scale_max.is_null() {
        if !w.is_null() {
            i = 0 as ::core::ffi::c_uint;
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
            i = 0 as ::core::ffi::c_uint;
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
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            ret += *w.offset(i as isize) * (*vec.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            ret += (*vec.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    }
    return ret;
}
unsafe fn diff_norm(
    mut n: ::core::ffi::c_uint,
    mut x: *const ::core::ffi::c_double,
    mut oldx: *const ::core::ffi::c_double,
    mut w: *const ::core::ffi::c_double,
    mut scale_min: *const ::core::ffi::c_double,
    mut scale_max: *const ::core::ffi::c_double,
) -> ::core::ffi::c_double {
    let mut i: ::core::ffi::c_uint = 0;
    let mut ret: ::core::ffi::c_double = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    if !scale_min.is_null() && !scale_max.is_null() {
        if !w.is_null() {
            i = 0 as ::core::ffi::c_uint;
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
            i = 0 as ::core::ffi::c_uint;
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
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            ret += *w.offset(i as isize) * (*x.offset(i as isize) - *oldx.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            ret += (*x.offset(i as isize) - *oldx.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    }
    return ret;
}
unsafe fn relstop(
    mut vold: ::core::ffi::c_double,
    mut vnew: ::core::ffi::c_double,
    mut reltol: ::core::ffi::c_double,
    mut abstol: ::core::ffi::c_double,
) -> ::core::ffi::c_int {
    if nlopt_isinf(vold) != 0 {
        return 0 as ::core::ffi::c_int;
    }
    return ((vnew - vold).abs() < abstol
        || (vnew - vold).abs() < reltol * ((vnew).abs() + (vold).abs()) * 0.5f64
        || reltol > 0 as ::core::ffi::c_int as ::core::ffi::c_double && vnew == vold) as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_ftol(
    mut s: *const nlopt_stopping,
    mut f: ::core::ffi::c_double,
    mut oldf: ::core::ffi::c_double,
) -> ::core::ffi::c_int {
    return relstop(oldf, f, (*s).ftol_rel, (*s).ftol_abs);
}

unsafe fn nlopt_stop_f(
    mut s: *const nlopt_stopping,
    mut f: ::core::ffi::c_double,
    mut oldf: ::core::ffi::c_double,
) -> ::core::ffi::c_int {
    return (f <= (*s).minf_max || nlopt_stop_ftol(s, f, oldf) != 0) as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_x(
    mut s: *const nlopt_stopping,
    mut x: *const ::core::ffi::c_double,
    mut oldx: *const ::core::ffi::c_double,
) -> ::core::ffi::c_int {
    let mut i: ::core::ffi::c_uint = 0;
    if diff_norm(
        (*s).n,
        x,
        oldx,
        (*s).x_weights,
        ::core::ptr::null::<::core::ffi::c_double>(),
        ::core::ptr::null::<::core::ffi::c_double>(),
    ) < (*s).xtol_rel
        * vector_norm(
            (*s).n,
            x,
            (*s).x_weights,
            ::core::ptr::null::<::core::ffi::c_double>(),
            ::core::ptr::null::<::core::ffi::c_double>(),
        )
    {
        return 1 as ::core::ffi::c_int;
    }
    if (*s).xtol_abs.is_null() {
        return 0 as ::core::ffi::c_int;
    }
    i = 0 as ::core::ffi::c_uint;
    while i < (*s).n {
        if (*x.offset(i as isize) - *oldx.offset(i as isize)).abs()
            >= *((*s).xtol_abs).offset(i as isize)
        {
            return 0 as ::core::ffi::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_dx(
    mut s: *const nlopt_stopping,
    mut x: *const ::core::ffi::c_double,
    mut dx: *const ::core::ffi::c_double,
) -> ::core::ffi::c_int {
    let mut i: ::core::ffi::c_uint = 0;
    if vector_norm(
        (*s).n,
        dx,
        (*s).x_weights,
        ::core::ptr::null::<::core::ffi::c_double>(),
        ::core::ptr::null::<::core::ffi::c_double>(),
    ) < (*s).xtol_rel
        * vector_norm(
            (*s).n,
            x,
            (*s).x_weights,
            ::core::ptr::null::<::core::ffi::c_double>(),
            ::core::ptr::null::<::core::ffi::c_double>(),
        )
    {
        return 1 as ::core::ffi::c_int;
    }
    if (*s).xtol_abs.is_null() {
        return 0 as ::core::ffi::c_int;
    }
    i = 0 as ::core::ffi::c_uint;
    while i < (*s).n {
        if (*dx.offset(i as isize)).abs() >= *((*s).xtol_abs).offset(i as isize) {
            return 0 as ::core::ffi::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_xs(
    mut s: *const nlopt_stopping,
    mut xs: *const ::core::ffi::c_double,
    mut oldxs: *const ::core::ffi::c_double,
    mut scale_min: *const ::core::ffi::c_double,
    mut scale_max: *const ::core::ffi::c_double,
) -> ::core::ffi::c_int {
    let mut i: ::core::ffi::c_uint = 0;
    if diff_norm((*s).n, xs, oldxs, (*s).x_weights, scale_min, scale_max)
        < (*s).xtol_rel * vector_norm((*s).n, xs, (*s).x_weights, scale_min, scale_max)
    {
        return 1 as ::core::ffi::c_int;
    }
    if (*s).xtol_abs.is_null() {
        return 0 as ::core::ffi::c_int;
    }
    i = 0 as ::core::ffi::c_uint;
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
            return 0 as ::core::ffi::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as ::core::ffi::c_int;
}

unsafe fn nlopt_isfinite(mut x: ::core::ffi::c_double) -> ::core::ffi::c_int {
    return ((x).abs() <= 1.7976931348623157e+308f64) as ::core::ffi::c_int;
}

unsafe fn nlopt_istiny(mut x: ::core::ffi::c_double) -> ::core::ffi::c_int {
    if x == 0.0f64 {
        return 1 as ::core::ffi::c_int;
    } else {
        return ((x).abs() < 2.2250738585072014e-308f64) as ::core::ffi::c_int;
    };
}

unsafe fn nlopt_isnan(mut x: ::core::ffi::c_double) -> ::core::ffi::c_int {
    return x.is_nan() as i32;
}

unsafe fn nlopt_stop_evals(mut s: *const nlopt_stopping) -> ::core::ffi::c_int {
    return ((*s).maxeval > 0 as ::core::ffi::c_int && *(*s).nevals_p >= (*s).maxeval) as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_time_(mut start: ::core::ffi::c_double, mut maxtime: ::core::ffi::c_double) -> ::core::ffi::c_int {
    return (maxtime > 0 as ::core::ffi::c_int as ::core::ffi::c_double && nlopt_seconds() - start >= maxtime)
        as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_time(mut s: *const nlopt_stopping) -> ::core::ffi::c_int {
    return nlopt_stop_time_((*s).start, (*s).maxtime);
}

unsafe fn nlopt_stop_evalstime(mut stop: *const nlopt_stopping) -> ::core::ffi::c_int {
    return (nlopt_stop_evals(stop) != 0 || nlopt_stop_time(stop) != 0) as ::core::ffi::c_int;
}

unsafe fn nlopt_stop_forced(mut stop: *const nlopt_stopping) -> ::core::ffi::c_int {
    return (!((*stop).force_stop).is_null() && *(*stop).force_stop != 0) as ::core::ffi::c_int;
}
//
// pub unsafe fn nlopt_vsprintf(
//     mut p: *mut ::core::ffi::c_char,
//     mut format: *const ::core::ffi::c_char,
//     mut ap: ::std::ffi::VaList,
// ) -> *mut ::core::ffi::c_char {
//     let mut len: size_t = (strlen(format))
//         .wrapping_add(128 as ::core::ffi::c_int as ::core::ffi::c_ulong);
//     let mut ret: ::core::ffi::c_int = 0;
//     p = realloc(p as *mut ::core::ffi::c_void, len) as *mut ::core::ffi::c_char;
//     if p.is_null() {
//         abort();
//     }
//     loop {
//         ret = vsnprintf(p, len, format, ap.as_va_list());
//         if !(ret < 0 as ::core::ffi::c_int || ret as size_t >= len) {
//             break;
//         }
//         len = if ret >= 0 as ::core::ffi::c_int {
//             (ret + 1 as ::core::ffi::c_int) as size_t
//         } else {
//             len.wrapping_mul(3 as ::core::ffi::c_int as ::core::ffi::c_ulong) >> 1 as ::core::ffi::c_int
//         };
//         p = realloc(p as *mut ::core::ffi::c_void, len) as *mut ::core::ffi::c_char;
//         if p.is_null() {
//             abort();
//         }
//     }
//     return p;
// }

unsafe fn nlopt_stop_msg(mut s: *mut nlopt_stopping, msg: &str) {
    (*s).stop_msg = msg.to_string();
}

unsafe fn nlopt_count_constraints(
    mut p: ::core::ffi::c_uint,
    mut c: *const nlopt_constraint,
) -> ::core::ffi::c_uint {
    let mut i: ::core::ffi::c_uint = 0;
    let mut count: ::core::ffi::c_uint = 0 as ::core::ffi::c_uint;
    i = 0 as ::core::ffi::c_uint;
    while i < p {
        count = count.wrapping_add((*c.offset(i as isize)).m);
        i = i.wrapping_add(1);
    }
    return count;
}

unsafe fn nlopt_max_constraint_dim(
    mut p: ::core::ffi::c_uint,
    mut c: *const nlopt_constraint,
) -> ::core::ffi::c_uint {
    let mut i: ::core::ffi::c_uint = 0;
    let mut max_dim: ::core::ffi::c_uint = 0 as ::core::ffi::c_uint;
    i = 0 as ::core::ffi::c_uint;
    while i < p {
        if (*c.offset(i as isize)).m > max_dim {
            max_dim = (*c.offset(i as isize)).m;
        }
        i = i.wrapping_add(1);
    }
    return max_dim;
}

unsafe fn nlopt_eval_constraint<U>(
    mut result: *mut ::core::ffi::c_double,
    mut grad: *mut ::core::ffi::c_double,
    mut c: *const nlopt_constraint,
    mut n: ::core::ffi::c_uint,
    mut x: *const ::core::ffi::c_double,
) {
    if (*c).f.is_some() {
        *result.offset(0 as ::core::ffi::c_int as isize) =
        // PATCH Weird bug ((*c).f).expect("non-null function pointer") calls the objective function!!!
        // even if (*c), nlopt_constraint object was correctly built with a nlopt_constraint_raw_callback!!! 
        //    ((*c).f).expect("non-null function pointer")(n, x, grad, (*c).f_data);
        // Maybe the U generic parameter required explains it cannot work like with C ???
        nlopt_constraint_raw_callback::<&dyn Func<U>, U>(n, x, grad, (*c).f_data);
    } else {
        (*c).mf.expect("non-null function pointer")((*c).m, result, n, x, grad, (*c).f_data);
    };
}

unsafe fn nlopt_isinf(mut x: ::core::ffi::c_double) -> ::core::ffi::c_int {
    return ((x).abs() >= f64::INFINITY * 0.99f64
        || if x.is_infinite() {
            if x.is_sign_positive() {
                1
            } else {
                -1
            }
        } else {
            0
        } != 0) as ::core::ffi::c_int;
}

unsafe fn nlopt_compute_rescaling(
    mut n: ::core::ffi::c_uint,
    mut dx: *const ::core::ffi::c_double,
) -> *mut ::core::ffi::c_double {
    // let mut s: *mut ::core::ffi::c_double = malloc(
    //     (::std::mem::size_of::<::core::ffi::c_double>() as ::core::ffi::c_ulong).wrapping_mul(n as ::core::ffi::c_ulong),
    // ) as *mut ::core::ffi::c_double;

    let mut space: Box<Vec<::core::ffi::c_double>> = Box::new(vec![0.; usize::try_from(n).unwrap()]);
    let s = space.as_mut_ptr() as *mut ::core::ffi::c_double;
    std::mem::forget(space);

    let mut i: ::core::ffi::c_uint = 0;
    if s.is_null() {
        return ::core::ptr::null_mut::<::core::ffi::c_double>();
    }
    i = 0 as ::core::ffi::c_uint;
    while i < n {
        *s.offset(i as isize) = 1.0f64;
        i = i.wrapping_add(1);
    }
    if n == 1 as ::core::ffi::c_uint {
        return s;
    }
    i = 1 as ::core::ffi::c_uint;
    while i < n
        && *dx.offset(i as isize) == *dx.offset(i.wrapping_sub(1 as ::core::ffi::c_uint) as isize)
    {
        i = i.wrapping_add(1);
    }
    if i < n {
        i = 1 as ::core::ffi::c_uint;
        while i < n {
            *s.offset(i as isize) =
                *dx.offset(i as isize) / *dx.offset(0 as ::core::ffi::c_int as isize);
            i = i.wrapping_add(1);
        }
    }
    return s;
}

unsafe fn nlopt_rescale(
    mut n: ::core::ffi::c_uint,
    mut s: *const ::core::ffi::c_double,
    mut x: *const ::core::ffi::c_double,
    mut xs: *mut ::core::ffi::c_double,
) {
    let mut i: ::core::ffi::c_uint = 0;
    if s.is_null() {
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize);
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize) / *s.offset(i as isize);
            i = i.wrapping_add(1);
        }
    };
}

unsafe fn nlopt_unscale(
    mut n: ::core::ffi::c_uint,
    mut s: *const ::core::ffi::c_double,
    mut x: *const ::core::ffi::c_double,
    mut xs: *mut ::core::ffi::c_double,
) {
    let mut i: ::core::ffi::c_uint = 0;
    if s.is_null() {
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize);
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as ::core::ffi::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize) * *s.offset(i as isize);
            i = i.wrapping_add(1);
        }
    };
}

unsafe fn nlopt_new_rescaled(
    mut n: ::core::ffi::c_uint,
    mut s: *const ::core::ffi::c_double,
    mut x: *const ::core::ffi::c_double,
) -> *mut ::core::ffi::c_double {
    // let mut xs: *mut ::core::ffi::c_double = malloc(
    //     (::std::mem::size_of::<::core::ffi::c_double>() as ::core::ffi::c_ulong).wrapping_mul(n as ::core::ffi::c_ulong),
    // ) as *mut ::core::ffi::c_double;

    let mut space: Box<Vec<::core::ffi::c_double>> = Box::new(vec![0.; usize::try_from(n).unwrap()]);
    let xs = space.as_mut_ptr() as *mut ::core::ffi::c_double;
    std::mem::forget(space);

    if xs.is_null() {
        return ::core::ptr::null_mut::<::core::ffi::c_double>();
    }
    nlopt_rescale(n, s, x, xs);
    return xs;
}

unsafe fn nlopt_reorder_bounds(
    mut n: ::core::ffi::c_uint,
    mut lb: *mut ::core::ffi::c_double,
    mut ub: *mut ::core::ffi::c_double,
) {
    let mut i: ::core::ffi::c_uint = 0;
    i = 0 as ::core::ffi::c_uint;
    while i < n {
        if *lb.offset(i as isize) > *ub.offset(i as isize) {
            let mut t: ::core::ffi::c_double = *lb.offset(i as isize);
            *lb.offset(i as isize) = *ub.offset(i as isize);
            *ub.offset(i as isize) = t;
        }
        i = i.wrapping_add(1);
    }
}
unsafe fn dcopy___(
    mut n_: *mut ::core::ffi::c_int,
    mut dx: *const ::core::ffi::c_double,
    mut incx: ::core::ffi::c_int,
    mut dy: *mut ::core::ffi::c_double,
    mut incy: ::core::ffi::c_int,
) {
    let mut i: ::core::ffi::c_int = 0;
    let mut n: ::core::ffi::c_int = *n_;
    if n <= 0 as ::core::ffi::c_int {
        return;
    }
    if incx == 1 as ::core::ffi::c_int && incy == 1 as ::core::ffi::c_int {
        memcpy(
            dy as *mut ::core::ffi::c_void,
            dx as *const ::core::ffi::c_void,
            (::core::mem::size_of::<::core::ffi::c_double>() as size_t)
                .wrapping_mul(n as ::core::ffi::c_uint as size_t),
        );
    } else if incx == 0 as ::core::ffi::c_int && incy == 1 as ::core::ffi::c_int {
        let mut x: ::core::ffi::c_double = *dx.offset(0 as ::core::ffi::c_int as isize);
        i = 0 as ::core::ffi::c_int;
        while i < n {
            *dy.offset(i as isize) = x;
            i += 1;
        }
    } else {
        i = 0 as ::core::ffi::c_int;
        while i < n {
            *dy.offset((i * incy) as isize) = *dx.offset((i * incx) as isize);
            i += 1;
        }
    };
}
unsafe fn daxpy_sl__(
    mut n_: *mut ::core::ffi::c_int,
    mut da_: *const ::core::ffi::c_double,
    mut dx: *const ::core::ffi::c_double,
    mut incx: ::core::ffi::c_int,
    mut dy: *mut ::core::ffi::c_double,
    mut incy: ::core::ffi::c_int,
) {
    let mut n: ::core::ffi::c_int = *n_;
    let mut i: ::core::ffi::c_int = 0;
    let mut da: ::core::ffi::c_double = *da_;
    if n <= 0 as ::core::ffi::c_int || da == 0 as ::core::ffi::c_int as ::core::ffi::c_double {
        return;
    }
    i = 0 as ::core::ffi::c_int;
    while i < n {
        *dy.offset((i * incy) as isize) += da * *dx.offset((i * incx) as isize);
        i += 1;
    }
}
unsafe fn ddot_sl__(
    mut n_: *mut ::core::ffi::c_int,
    mut dx: *mut ::core::ffi::c_double,
    mut incx: ::core::ffi::c_int,
    mut dy: *mut ::core::ffi::c_double,
    mut incy: ::core::ffi::c_int,
) -> ::core::ffi::c_double {
    let mut n: ::core::ffi::c_int = *n_;
    let mut i: ::core::ffi::c_int = 0;
    let mut sum: ::core::ffi::c_double = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    if n <= 0 as ::core::ffi::c_int {
        return 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    }
    i = 0 as ::core::ffi::c_int;
    while i < n {
        sum += *dx.offset((i * incx) as isize) * *dy.offset((i * incy) as isize);
        i += 1;
    }
    return sum;
}
unsafe fn dnrm2___(
    mut n_: *mut ::core::ffi::c_int,
    mut dx: *mut ::core::ffi::c_double,
    mut incx: ::core::ffi::c_int,
) -> ::core::ffi::c_double {
    let mut i: ::core::ffi::c_int = 0;
    let mut n: ::core::ffi::c_int = *n_;
    let mut xmax: ::core::ffi::c_double = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    let mut scale: ::core::ffi::c_double = 0.;
    let mut sum: ::core::ffi::c_double = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    i = 0 as ::core::ffi::c_int;
    while i < n {
        let mut xabs: ::core::ffi::c_double = (*dx.offset((incx * i) as isize)).abs();
        if xmax < xabs {
            xmax = xabs;
        }
        i += 1;
    }
    if xmax == 0 as ::core::ffi::c_int as ::core::ffi::c_double {
        return 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    }
    scale = 1.0f64 / xmax;
    i = 0 as ::core::ffi::c_int;
    while i < n {
        let mut xs: ::core::ffi::c_double = scale * *dx.offset((incx * i) as isize);
        sum += xs * xs;
        i += 1;
    }
    return xmax * sum.sqrt();
}
unsafe fn dsrot_(
    mut n: ::core::ffi::c_int,
    mut dx: *mut ::core::ffi::c_double,
    mut incx: ::core::ffi::c_int,
    mut dy: *mut ::core::ffi::c_double,
    mut incy: ::core::ffi::c_int,
    mut c__: *mut ::core::ffi::c_double,
    mut s_: *mut ::core::ffi::c_double,
) {
    let mut i: ::core::ffi::c_int = 0;
    let mut c: ::core::ffi::c_double = *c__;
    let mut s: ::core::ffi::c_double = *s_;
    i = 0 as ::core::ffi::c_int;
    while i < n {
        let mut x: ::core::ffi::c_double = *dx.offset((incx * i) as isize);
        let mut y: ::core::ffi::c_double = *dy.offset((incy * i) as isize);
        *dx.offset((incx * i) as isize) = c * x + s * y;
        *dy.offset((incy * i) as isize) = c * y - s * x;
        i += 1;
    }
}
unsafe fn dsrotg_(
    mut da: *mut ::core::ffi::c_double,
    mut db: *mut ::core::ffi::c_double,
    mut c: *mut ::core::ffi::c_double,
    mut s: *mut ::core::ffi::c_double,
) {
    let mut absa: ::core::ffi::c_double = 0.;
    let mut absb: ::core::ffi::c_double = 0.;
    let mut roe: ::core::ffi::c_double = 0.;
    let mut scale: ::core::ffi::c_double = 0.;
    absa = (*da).abs();
    absb = (*db).abs();
    if absa > absb {
        roe = *da;
        scale = absa;
    } else {
        roe = *db;
        scale = absb;
    }
    if scale != 0 as ::core::ffi::c_int as ::core::ffi::c_double {
        let mut r: ::core::ffi::c_double = 0.;
        let mut iscale: ::core::ffi::c_double = 1 as ::core::ffi::c_int as ::core::ffi::c_double / scale;
        let mut tmpa: ::core::ffi::c_double = *da * iscale;
        let mut tmpb: ::core::ffi::c_double = *db * iscale;
        r = (if roe < 0 as ::core::ffi::c_int as ::core::ffi::c_double {
            -scale
        } else {
            scale
        }) * (tmpa * tmpa + tmpb * tmpb).sqrt();
        *c = *da / r;
        *s = *db / r;
        *da = r;
        if *c != 0 as ::core::ffi::c_int as ::core::ffi::c_double && (*c).abs() <= *s {
            *db = 1 as ::core::ffi::c_int as ::core::ffi::c_double / *c;
        } else {
            *db = *s;
        }
    } else {
        *c = 1 as ::core::ffi::c_int as ::core::ffi::c_double;
        *db = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
        *da = *db;
        *s = *da;
    };
}
unsafe fn dscal_sl__(
    mut n_: *mut ::core::ffi::c_int,
    mut da: *const ::core::ffi::c_double,
    mut dx: *mut ::core::ffi::c_double,
    mut incx: ::core::ffi::c_int,
) {
    let mut i: ::core::ffi::c_int = 0;
    let mut n: ::core::ffi::c_int = *n_;
    let mut alpha: ::core::ffi::c_double = *da;
    i = 0 as ::core::ffi::c_int;
    while i < n {
        *dx.offset((i * incx) as isize) *= alpha;
        i += 1;
    }
}
static mut c__0: ::core::ffi::c_int = 0 as ::core::ffi::c_int;
static mut c__1: ::core::ffi::c_int = 1 as ::core::ffi::c_int;
static mut c__2: ::core::ffi::c_int = 2 as ::core::ffi::c_int;
unsafe fn h12_(
    mut mode: *const ::core::ffi::c_int,
    mut lpivot: *mut ::core::ffi::c_int,
    mut l1: *mut ::core::ffi::c_int,
    mut m: *mut ::core::ffi::c_int,
    mut u: *mut ::core::ffi::c_double,
    mut iue: *const ::core::ffi::c_int,
    mut up: *mut ::core::ffi::c_double,
    mut c__: *mut ::core::ffi::c_double,
    mut ice: *const ::core::ffi::c_int,
    mut icv: *const ::core::ffi::c_int,
    mut ncv: *const ::core::ffi::c_int,
) {
    let mut current_block: u64;
    let one: ::core::ffi::c_double = 1.0f64;
    let mut u_dim1: ::core::ffi::c_int = 0;
    let mut u_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut b: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut i2: ::core::ffi::c_int = 0;
    let mut i3: ::core::ffi::c_int = 0;
    let mut i4: ::core::ffi::c_int = 0;
    let mut cl: ::core::ffi::c_double = 0.;
    let mut sm: ::core::ffi::c_double = 0.;
    let mut incr: ::core::ffi::c_int = 0;
    let mut clinv: ::core::ffi::c_double = 0.;
    u_dim1 = *iue;
    u_offset = 1 as ::core::ffi::c_int + u_dim1;
    u = u.offset(-(u_offset as isize));
    c__ = c__.offset(-1);
    if !(0 as ::core::ffi::c_int >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
        d__1 = *u.offset((*lpivot * u_dim1 + 1 as ::core::ffi::c_int) as isize);
        cl = (d__1).abs();
        if *mode == 2 as ::core::ffi::c_int {
            if cl <= 0.0f64 {
                current_block = 17095853620682146675;
            } else {
                current_block = 17641539699934407521;
            }
        } else {
            i__1 = *m;
            j = *l1;
            while j <= i__1 {
                d__1 = *u.offset((j * u_dim1 + 1 as ::core::ffi::c_int) as isize);
                sm = (d__1).abs();
                cl = if sm >= cl { sm } else { cl };
                j += 1;
            }
            if cl <= 0.0f64 {
                current_block = 17095853620682146675;
            } else {
                clinv = one / cl;
                d__1 = *u.offset((*lpivot * u_dim1 + 1 as ::core::ffi::c_int) as isize) * clinv;
                sm = d__1 * d__1;
                i__1 = *m;
                j = *l1;
                while j <= i__1 {
                    d__1 = *u.offset((j * u_dim1 + 1 as ::core::ffi::c_int) as isize) * clinv;
                    sm += d__1 * d__1;
                    j += 1;
                }
                cl *= (sm).sqrt();
                if *u.offset((*lpivot * u_dim1 + 1 as ::core::ffi::c_int) as isize) > 0.0f64 {
                    cl = -cl;
                }
                *up = *u.offset((*lpivot * u_dim1 + 1 as ::core::ffi::c_int) as isize) - cl;
                *u.offset((*lpivot * u_dim1 + 1 as ::core::ffi::c_int) as isize) = cl;
                current_block = 17641539699934407521;
            }
        }
        match current_block {
            17095853620682146675 => {}
            _ => {
                if !(*ncv <= 0 as ::core::ffi::c_int) {
                    b = *up * *u.offset((*lpivot * u_dim1 + 1 as ::core::ffi::c_int) as isize);
                    if !(b >= 0.0f64) {
                        b = one / b;
                        i2 = 1 as ::core::ffi::c_int - *icv + *ice * (*lpivot - 1 as ::core::ffi::c_int);
                        incr = *ice * (*l1 - *lpivot);
                        i__1 = *ncv;
                        j = 1 as ::core::ffi::c_int;
                        while j <= i__1 {
                            i2 += *icv;
                            i3 = i2 + incr;
                            i4 = i3;
                            sm = *c__.offset(i2 as isize) * *up;
                            i__2 = *m;
                            i__ = *l1;
                            while i__ <= i__2 {
                                sm += *c__.offset(i3 as isize)
                                    * *u.offset((i__ * u_dim1 + 1 as ::core::ffi::c_int) as isize);
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
                                        sm * *u.offset((i__ * u_dim1 + 1 as ::core::ffi::c_int) as isize);
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
    mut a: *mut ::core::ffi::c_double,
    mut mda: *mut ::core::ffi::c_int,
    mut m: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut b: *mut ::core::ffi::c_double,
    mut x: *mut ::core::ffi::c_double,
    mut rnorm: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
    mut z__: *mut ::core::ffi::c_double,
    mut indx: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
) {
    let mut current_block: u64;
    let one: ::core::ffi::c_double = 1.0f64;
    let factor: ::core::ffi::c_double = 0.01f64;
    let mut a_dim1: ::core::ffi::c_int = 0;
    let mut a_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut c__: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut k: ::core::ffi::c_int = 0;
    let mut l: ::core::ffi::c_int = 0;
    let mut s: ::core::ffi::c_double = 0.;
    let mut t: ::core::ffi::c_double = 0.;
    let mut ii: ::core::ffi::c_int = 0;
    let mut jj: ::core::ffi::c_int = 0;
    let mut ip: ::core::ffi::c_int = 0;
    let mut iz: ::core::ffi::c_int = 0;
    let mut jz: ::core::ffi::c_int = 0;
    let mut up: ::core::ffi::c_double = 0.;
    let mut iz1: ::core::ffi::c_int = 0;
    let mut iz2: ::core::ffi::c_int = 0;
    let mut npp1: ::core::ffi::c_int = 0;
    let mut iter: ::core::ffi::c_int = 0;
    let mut wmax: ::core::ffi::c_double = 0.;
    let mut alpha: ::core::ffi::c_double = 0.;
    let mut asave: ::core::ffi::c_double = 0.;
    let mut itmax: ::core::ffi::c_int = 0;
    let mut izmax: ::core::ffi::c_int = 0 as ::core::ffi::c_int;
    let mut nsetp: ::core::ffi::c_int = 0;
    let mut unorm: ::core::ffi::c_double = 0.;
    z__ = z__.offset(-1);
    b = b.offset(-1);
    indx = indx.offset(-1);
    w = w.offset(-1);
    x = x.offset(-1);
    a_dim1 = *mda;
    a_offset = 1 as ::core::ffi::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    *mode = 2 as ::core::ffi::c_int;
    if !(*m <= 0 as ::core::ffi::c_int || *n <= 0 as ::core::ffi::c_int) {
        *mode = 1 as ::core::ffi::c_int;
        iter = 0 as ::core::ffi::c_int;
        itmax = *n * 3 as ::core::ffi::c_int;
        i__1 = *n;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            *indx.offset(i__ as isize) = i__;
            i__ += 1;
        }
        iz1 = 1 as ::core::ffi::c_int;
        iz2 = *n;
        nsetp = 0 as ::core::ffi::c_int;
        npp1 = 1 as ::core::ffi::c_int;
        *x.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
        dcopy___(
            n,
            x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        '_L110: loop {
            if iz1 > iz2 || nsetp >= *m {
                current_block = 1598957747739433155;
                break;
            }
            i__1 = iz2;
            iz = iz1;
            while iz <= i__1 {
                j = *indx.offset(iz as isize);
                i__2 = *m - nsetp;
                *w.offset(j as isize) = ddot_sl__(
                    &raw mut i__2,
                    a.offset((npp1 + j * a_dim1) as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    b.offset(npp1 as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
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
                    current_block = 1598957747739433155;
                    break '_L110;
                }
                iz = izmax;
                j = *indx.offset(iz as isize);
                asave = *a.offset((npp1 + j * a_dim1) as isize);
                i__2 = npp1 + 1 as ::core::ffi::c_int;
                h12_(
                    &raw const c__1,
                    &raw mut npp1,
                    &raw mut i__2,
                    m,
                    a.offset((j * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                        as *mut ::core::ffi::c_double,
                    &raw const c__1,
                    &raw mut up,
                    z__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    &raw const c__1,
                    &raw const c__1,
                    &raw const c__0,
                );
                unorm = dnrm2___(
                    &raw mut nsetp,
                    a.offset((j * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                        as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                d__1 = *a.offset((npp1 + j * a_dim1) as isize);
                t = factor * (d__1).abs();
                d__1 = unorm + t;
                if !(d__1 - unorm <= 0.0f64) {
                    dcopy___(
                        m,
                        b.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        1 as ::core::ffi::c_int,
                        z__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        1 as ::core::ffi::c_int,
                    );
                    i__2 = npp1 + 1 as ::core::ffi::c_int;
                    h12_(
                        &raw const c__2,
                        &raw mut npp1,
                        &raw mut i__2,
                        m,
                        a.offset((j * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                            as *mut ::core::ffi::c_double,
                        &raw const c__1,
                        &raw mut up,
                        z__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        &raw const c__1,
                        &raw const c__1,
                        &raw const c__1,
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
                z__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
                b.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
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
                    &raw const c__2,
                    &raw mut nsetp,
                    &raw mut npp1,
                    m,
                    a.offset((j * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                        as *mut ::core::ffi::c_double,
                    &raw const c__1,
                    &raw mut up,
                    a.offset((jj * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                        as *mut ::core::ffi::c_double,
                    &raw const c__1,
                    mda,
                    &raw const c__1,
                );
                jz += 1;
            }
            k = if npp1 <= *mda { npp1 } else { *mda };
            *w.offset(j as isize) = 0.0f64;
            i__2 = *m - nsetp;
            dcopy___(
                &raw mut i__2,
                w.offset(j as isize) as *mut ::core::ffi::c_double,
                0 as ::core::ffi::c_int,
                a.offset((k + j * a_dim1) as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
            );
            loop {
                ip = nsetp;
                while ip >= 1 as ::core::ffi::c_int {
                    if !(ip == nsetp) {
                        d__1 = -*z__.offset((ip + 1 as ::core::ffi::c_int) as isize);
                        daxpy_sl__(
                            &raw mut ip,
                            &raw mut d__1,
                            a.offset((jj * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            z__.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                    }
                    jj = *indx.offset(ip as isize);
                    *z__.offset(ip as isize) /= *a.offset((ip + jj * a_dim1) as isize);
                    ip -= 1;
                }
                iter += 1;
                if !(iter <= itmax) {
                    current_block = 15774334763445529798;
                    break '_L110;
                }
                alpha = one;
                jj = 0 as ::core::ffi::c_int;
                i__2 = nsetp;
                ip = 1 as ::core::ffi::c_int;
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
                ip = 1 as ::core::ffi::c_int;
                while ip <= i__2 {
                    l = *indx.offset(ip as isize);
                    *x.offset(l as isize) =
                        (one - alpha) * *x.offset(l as isize) + alpha * *z__.offset(ip as isize);
                    ip += 1;
                }
                if jj == 0 as ::core::ffi::c_int {
                    break;
                }
                i__ = *indx.offset(jj as isize);
                '_L250: loop {
                    *x.offset(i__ as isize) = 0.0f64;
                    jj += 1;
                    i__2 = nsetp;
                    j = jj;
                    while j <= i__2 {
                        ii = *indx.offset(j as isize);
                        *indx.offset((j - 1 as ::core::ffi::c_int) as isize) = ii;
                        dsrotg_(
                            &mut *a.offset((j - 1 as ::core::ffi::c_int + ii * a_dim1) as isize),
                            &mut *a.offset((j + ii * a_dim1) as isize),
                            &mut c__,
                            &mut s,
                        );
                        t = *a.offset((j - 1 as ::core::ffi::c_int + ii * a_dim1) as isize);
                        dsrot_(
                            *n,
                            a.offset((j - 1 as ::core::ffi::c_int + a_dim1) as isize)
                                as *mut ::core::ffi::c_double,
                            *mda,
                            a.offset((j + a_dim1) as isize) as *mut ::core::ffi::c_double,
                            *mda,
                            &raw mut c__,
                            &raw mut s,
                        );
                        *a.offset((j - 1 as ::core::ffi::c_int + ii * a_dim1) as isize) = t;
                        *a.offset((j + ii * a_dim1) as isize) = 0.0f64;
                        dsrot_(
                            1 as ::core::ffi::c_int,
                            b.offset((j - 1 as ::core::ffi::c_int) as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            b.offset(j as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            &raw mut c__,
                            &raw mut s,
                        );
                        j += 1;
                    }
                    npp1 = nsetp;
                    nsetp -= 1;
                    iz1 -= 1;
                    *indx.offset(iz1 as isize) = i__;
                    if nsetp <= 0 as ::core::ffi::c_int {
                        current_block = 15774334763445529798;
                        break '_L110;
                    }
                    i__2 = nsetp;
                    jj = 1 as ::core::ffi::c_int;
                    loop {
                        if !(jj <= i__2) {
                            break '_L250;
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
                    b.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    z__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
            }
        }
        match current_block {
            15774334763445529798 => {
                *mode = 3 as ::core::ffi::c_int;
            }
            _ => {}
        }
        k = if npp1 <= *m { npp1 } else { *m };
        i__2 = *m - nsetp;
        *rnorm = dnrm2___(
            &raw mut i__2,
            b.offset(k as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        if npp1 > *m {
            *w.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
            dcopy___(
                n,
                w.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                0 as ::core::ffi::c_int,
                w.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
            );
        }
    }
}
unsafe fn ldp_(
    mut g: *mut ::core::ffi::c_double,
    mut mg: *mut ::core::ffi::c_int,
    mut m: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut h__: *mut ::core::ffi::c_double,
    mut x: *mut ::core::ffi::c_double,
    mut xnorm: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
    mut indx: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
) {
    let one: ::core::ffi::c_double = 1.0f64;
    let mut g_dim1: ::core::ffi::c_int = 0;
    let mut g_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut n1: ::core::ffi::c_int = 0;
    let mut if__: ::core::ffi::c_int = 0;
    let mut iw: ::core::ffi::c_int = 0;
    let mut iy: ::core::ffi::c_int = 0;
    let mut iz: ::core::ffi::c_int = 0;
    let mut fac: ::core::ffi::c_double = 0.;
    let mut rnorm: ::core::ffi::c_double = 0.;
    let mut iwdual: ::core::ffi::c_int = 0;
    indx = indx.offset(-1);
    h__ = h__.offset(-1);
    x = x.offset(-1);
    g_dim1 = *mg;
    g_offset = 1 as ::core::ffi::c_int + g_dim1;
    g = g.offset(-(g_offset as isize));
    w = w.offset(-1);
    *mode = 2 as ::core::ffi::c_int;
    if !(*n <= 0 as ::core::ffi::c_int) {
        *mode = 1 as ::core::ffi::c_int;
        *x.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
        dcopy___(
            n,
            x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        *xnorm = 0.0f64;
        if !(*m == 0 as ::core::ffi::c_int) {
            iw = 0 as ::core::ffi::c_int;
            i__1 = *m;
            j = 1 as ::core::ffi::c_int;
            while j <= i__1 {
                i__2 = *n;
                i__ = 1 as ::core::ffi::c_int;
                while i__ <= i__2 {
                    iw += 1;
                    *w.offset(iw as isize) = *g.offset((j + i__ * g_dim1) as isize);
                    i__ += 1;
                }
                iw += 1;
                *w.offset(iw as isize) = *h__.offset(j as isize);
                j += 1;
            }
            if__ = iw + 1 as ::core::ffi::c_int;
            i__1 = *n;
            i__ = 1 as ::core::ffi::c_int;
            while i__ <= i__1 {
                iw += 1;
                *w.offset(iw as isize) = 0.0f64;
                i__ += 1;
            }
            *w.offset((iw + 1 as ::core::ffi::c_int) as isize) = one;
            n1 = *n + 1 as ::core::ffi::c_int;
            iz = iw + 2 as ::core::ffi::c_int;
            iy = iz + n1;
            iwdual = iy + *m;
            nnls_(
                w.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                &raw mut n1,
                &raw mut n1,
                m,
                w.offset(if__ as isize) as *mut ::core::ffi::c_double,
                w.offset(iy as isize) as *mut ::core::ffi::c_double,
                &raw mut rnorm,
                w.offset(iwdual as isize) as *mut ::core::ffi::c_double,
                w.offset(iz as isize) as *mut ::core::ffi::c_double,
                indx.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
                mode,
            );
            if !(*mode != 1 as ::core::ffi::c_int) {
                *mode = 4 as ::core::ffi::c_int;
                if !(rnorm <= 0.0f64) {
                    fac = one
                        - ddot_sl__(
                            m,
                            h__.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            w.offset(iy as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                    d__1 = one + fac;
                    if !(d__1 - one <= 0.0f64) {
                        *mode = 1 as ::core::ffi::c_int;
                        fac = one / fac;
                        i__1 = *n;
                        j = 1 as ::core::ffi::c_int;
                        while j <= i__1 {
                            *x.offset(j as isize) = fac
                                * ddot_sl__(
                                    m,
                                    g.offset((j * g_dim1 + 1 as ::core::ffi::c_int) as isize)
                                        as *mut ::core::ffi::c_double,
                                    1 as ::core::ffi::c_int,
                                    w.offset(iy as isize) as *mut ::core::ffi::c_double,
                                    1 as ::core::ffi::c_int,
                                );
                            j += 1;
                        }
                        *xnorm = dnrm2___(
                            n,
                            x.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                        *w.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
                        dcopy___(
                            m,
                            w.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            0 as ::core::ffi::c_int,
                            w.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                        daxpy_sl__(
                            m,
                            &raw mut fac,
                            w.offset(iy as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            w.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                    }
                }
            }
        }
    }
}
unsafe fn lsi_(
    mut e: *mut ::core::ffi::c_double,
    mut f: *mut ::core::ffi::c_double,
    mut g: *mut ::core::ffi::c_double,
    mut h__: *mut ::core::ffi::c_double,
    mut le: *mut ::core::ffi::c_int,
    mut me: *mut ::core::ffi::c_int,
    mut lg: *mut ::core::ffi::c_int,
    mut mg: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut x: *mut ::core::ffi::c_double,
    mut xnorm: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
    mut jw: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
) {
    let mut current_block: u64;
    let epmach: ::core::ffi::c_double = 2.22e-16f64;
    let one: ::core::ffi::c_double = 1.0f64;
    let mut e_dim1: ::core::ffi::c_int = 0;
    let mut e_offset: ::core::ffi::c_int = 0;
    let mut g_dim1: ::core::ffi::c_int = 0;
    let mut g_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut i__3: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut t: ::core::ffi::c_double = 0.;
    f = f.offset(-1);
    jw = jw.offset(-1);
    h__ = h__.offset(-1);
    x = x.offset(-1);
    g_dim1 = *lg;
    g_offset = 1 as ::core::ffi::c_int + g_dim1;
    g = g.offset(-(g_offset as isize));
    e_dim1 = *le;
    e_offset = 1 as ::core::ffi::c_int + e_dim1;
    e = e.offset(-(e_offset as isize));
    w = w.offset(-1);
    i__1 = *n;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= i__1 {
        i__2 = i__ + 1 as ::core::ffi::c_int;
        j = if i__2 <= *n { i__2 } else { *n };
        i__2 = i__ + 1 as ::core::ffi::c_int;
        i__3 = *n - i__;
        h12_(
            &raw const c__1,
            &raw mut i__,
            &raw mut i__2,
            me,
            e.offset((i__ * e_dim1 + 1 as ::core::ffi::c_int) as isize)
                as *mut ::core::ffi::c_double,
            &raw const c__1,
            &raw mut t,
            e.offset((j * e_dim1 + 1 as ::core::ffi::c_int) as isize) as *mut ::core::ffi::c_double,
            &raw const c__1,
            le,
            &raw mut i__3,
        );
        i__2 = i__ + 1 as ::core::ffi::c_int;
        h12_(
            &raw const c__2,
            &raw mut i__,
            &raw mut i__2,
            me,
            e.offset((i__ * e_dim1 + 1 as ::core::ffi::c_int) as isize)
                as *mut ::core::ffi::c_double,
            &raw const c__1,
            &raw mut t,
            f.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            &raw const c__1,
            &raw const c__1,
            &raw const c__1,
        );
        i__ += 1;
    }
    *mode = 5 as ::core::ffi::c_int;
    i__2 = *mg;
    i__ = 1 as ::core::ffi::c_int;
    's_92: loop {
        if !(i__ <= i__2) {
            current_block = 11057878835866523405;
            break;
        }
        i__1 = *n;
        j = 1 as ::core::ffi::c_int;
        while j <= i__1 {
            d__1 = *e.offset((j + j * e_dim1) as isize);
            if (d__1).abs() < epmach {
                current_block = 15870820249773205402;
                break 's_92;
            }
            i__3 = j - 1 as ::core::ffi::c_int;
            *g.offset((i__ + j * g_dim1) as isize) = (*g.offset((i__ + j * g_dim1) as isize)
                - ddot_sl__(
                    &raw mut i__3,
                    g.offset((i__ + g_dim1) as isize) as *mut ::core::ffi::c_double,
                    *lg,
                    e.offset((j * e_dim1 + 1 as ::core::ffi::c_int) as isize)
                        as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                ))
                / *e.offset((j + j * e_dim1) as isize);
            j += 1;
        }
        *h__.offset(i__ as isize) -= ddot_sl__(
            n,
            g.offset((i__ + g_dim1) as isize) as *mut ::core::ffi::c_double,
            *lg,
            f.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        i__ += 1;
    }
    match current_block {
        11057878835866523405 => {
            ldp_(
                g.offset(g_offset as isize) as *mut ::core::ffi::c_double,
                lg,
                mg,
                n,
                h__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                xnorm,
                w.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                jw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
                mode,
            );
            if !(*mode != 1 as ::core::ffi::c_int) {
                daxpy_sl__(
                    n,
                    &raw const one,
                    f.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                i__ = *n;
                while i__ >= 1 as ::core::ffi::c_int {
                    i__2 = i__ + 1 as ::core::ffi::c_int;
                    j = if i__2 <= *n { i__2 } else { *n };
                    i__2 = *n - i__;
                    *x.offset(i__ as isize) = (*x.offset(i__ as isize)
                        - ddot_sl__(
                            &raw mut i__2,
                            e.offset((i__ + j * e_dim1) as isize) as *mut ::core::ffi::c_double,
                            *le,
                            x.offset(j as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        ))
                        / *e.offset((i__ + i__ * e_dim1) as isize);
                    i__ -= 1;
                }
                i__2 = *n + 1 as ::core::ffi::c_int;
                j = if i__2 <= *me { i__2 } else { *me };
                i__2 = *me - *n;
                t = dnrm2___(
                    &raw mut i__2,
                    f.offset(j as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                *xnorm = (*xnorm * *xnorm + t * t).sqrt();
            }
        }
        _ => {}
    };
}
unsafe fn hfti_(
    mut a: *mut ::core::ffi::c_double,
    mut mda: *mut ::core::ffi::c_int,
    mut m: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut b: *mut ::core::ffi::c_double,
    mut mdb: *mut ::core::ffi::c_int,
    mut nb: *const ::core::ffi::c_int,
    mut tau: *mut ::core::ffi::c_double,
    mut krank: *mut ::core::ffi::c_int,
    mut rnorm: *mut ::core::ffi::c_double,
    mut h__: *mut ::core::ffi::c_double,
    mut g: *mut ::core::ffi::c_double,
    mut ip: *mut ::core::ffi::c_int,
) {
    let mut current_block: u64;
    let factor: ::core::ffi::c_double = 0.001f64;
    let mut a_dim1: ::core::ffi::c_int = 0;
    let mut a_offset: ::core::ffi::c_int = 0;
    let mut b_dim1: ::core::ffi::c_int = 0;
    let mut b_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut i__3: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut k: ::core::ffi::c_int = 0;
    let mut l: ::core::ffi::c_int = 0;
    let mut jb: ::core::ffi::c_int = 0;
    let mut kp1: ::core::ffi::c_int = 0;
    let mut tmp: ::core::ffi::c_double = 0.;
    let mut hmax: ::core::ffi::c_double = 0.;
    let mut lmax: ::core::ffi::c_int = 0;
    let mut ldiag: ::core::ffi::c_int = 0;
    ip = ip.offset(-1);
    g = g.offset(-1);
    h__ = h__.offset(-1);
    a_dim1 = *mda;
    a_offset = 1 as ::core::ffi::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    rnorm = rnorm.offset(-1);
    b_dim1 = *mdb;
    b_offset = 1 as ::core::ffi::c_int + b_dim1;
    b = b.offset(-(b_offset as isize));
    k = 0 as ::core::ffi::c_int;
    ldiag = if *m <= *n { *m } else { *n };
    if !(ldiag <= 0 as ::core::ffi::c_int) {
        i__1 = ldiag;
        j = 1 as ::core::ffi::c_int;
        while j <= i__1 {
            let mut current_block_56: u64;
            if j == 1 as ::core::ffi::c_int {
                current_block_56 = 14634468614531113631;
            } else {
                lmax = j;
                i__2 = *n;
                l = j;
                while l <= i__2 {
                    d__1 = *a.offset((j - 1 as ::core::ffi::c_int + l * a_dim1) as isize);
                    *h__.offset(l as isize) -= d__1 * d__1;
                    if *h__.offset(l as isize) > *h__.offset(lmax as isize) {
                        lmax = l;
                    }
                    l += 1;
                }
                d__1 = hmax + factor * *h__.offset(lmax as isize);
                if d__1 - hmax > 0.0f64 {
                    current_block_56 = 2432358212778639914;
                } else {
                    current_block_56 = 14634468614531113631;
                }
            }
            match current_block_56 {
                14634468614531113631 => {
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
                i__ = 1 as ::core::ffi::c_int;
                while i__ <= i__2 {
                    tmp = *a.offset((i__ + j * a_dim1) as isize);
                    *a.offset((i__ + j * a_dim1) as isize) =
                        *a.offset((i__ + lmax * a_dim1) as isize);
                    *a.offset((i__ + lmax * a_dim1) as isize) = tmp;
                    i__ += 1;
                }
                *h__.offset(lmax as isize) = *h__.offset(j as isize);
            }
            i__2 = j + 1 as ::core::ffi::c_int;
            i__ = if i__2 <= *n { i__2 } else { *n };
            i__2 = j + 1 as ::core::ffi::c_int;
            i__3 = *n - j;
            h12_(
                &raw const c__1,
                &raw mut j,
                &raw mut i__2,
                m,
                a.offset((j * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                    as *mut ::core::ffi::c_double,
                &raw const c__1,
                h__.offset(j as isize) as *mut ::core::ffi::c_double,
                a.offset((i__ * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                    as *mut ::core::ffi::c_double,
                &raw const c__1,
                mda,
                &raw mut i__3,
            );
            i__2 = j + 1 as ::core::ffi::c_int;
            h12_(
                &raw const c__2,
                &raw mut j,
                &raw mut i__2,
                m,
                a.offset((j * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                    as *mut ::core::ffi::c_double,
                &raw const c__1,
                h__.offset(j as isize) as *mut ::core::ffi::c_double,
                b.offset(b_offset as isize) as *mut ::core::ffi::c_double,
                &raw const c__1,
                mdb,
                nb,
            );
            j += 1;
        }
        i__2 = ldiag;
        j = 1 as ::core::ffi::c_int;
        loop {
            if !(j <= i__2) {
                current_block = 8545136480011357681;
                break;
            }
            d__1 = *a.offset((j + j * a_dim1) as isize);
            if (d__1).abs() <= *tau {
                current_block = 1748094236375284402;
                break;
            }
            j += 1;
        }
        match current_block {
            8545136480011357681 => {
                k = ldiag;
            }
            _ => {
                k = j - 1 as ::core::ffi::c_int;
            }
        }
        kp1 = k + 1 as ::core::ffi::c_int;
        i__2 = *nb;
        jb = 1 as ::core::ffi::c_int;
        while jb <= i__2 {
            i__1 = *m - k;
            *rnorm.offset(jb as isize) = dnrm2___(
                &raw mut i__1,
                b.offset((kp1 + jb * b_dim1) as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
            );
            jb += 1;
        }
        if k > 0 as ::core::ffi::c_int {
            if !(k == *n) {
                i__ = k;
                while i__ >= 1 as ::core::ffi::c_int {
                    i__2 = i__ - 1 as ::core::ffi::c_int;
                    h12_(
                        &raw const c__1,
                        &raw mut i__,
                        &raw mut kp1,
                        n,
                        a.offset((i__ + a_dim1) as isize) as *mut ::core::ffi::c_double,
                        mda,
                        g.offset(i__ as isize) as *mut ::core::ffi::c_double,
                        a.offset(a_offset as isize) as *mut ::core::ffi::c_double,
                        mda,
                        &raw const c__1,
                        &raw mut i__2,
                    );
                    i__ -= 1;
                }
            }
            i__2 = *nb;
            jb = 1 as ::core::ffi::c_int;
            while jb <= i__2 {
                i__ = k;
                while i__ >= 1 as ::core::ffi::c_int {
                    i__1 = i__ + 1 as ::core::ffi::c_int;
                    j = if i__1 <= *n { i__1 } else { *n };
                    i__1 = k - i__;
                    *b.offset((i__ + jb * b_dim1) as isize) = (*b
                        .offset((i__ + jb * b_dim1) as isize)
                        - ddot_sl__(
                            &raw mut i__1,
                            a.offset((i__ + j * a_dim1) as isize) as *mut ::core::ffi::c_double,
                            *mda,
                            b.offset((j + jb * b_dim1) as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
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
                    i__ = 1 as ::core::ffi::c_int;
                    while i__ <= i__1 {
                        h12_(
                            &raw const c__2,
                            &raw mut i__,
                            &raw mut kp1,
                            n,
                            a.offset((i__ + a_dim1) as isize) as *mut ::core::ffi::c_double,
                            mda,
                            g.offset(i__ as isize) as *mut ::core::ffi::c_double,
                            b.offset((jb * b_dim1 + 1 as ::core::ffi::c_int) as isize)
                                as *mut ::core::ffi::c_double,
                            &raw const c__1,
                            mdb,
                            &raw const c__1,
                        );
                        i__ += 1;
                    }
                }
                j = ldiag;
                while j >= 1 as ::core::ffi::c_int {
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
            jb = 1 as ::core::ffi::c_int;
            while jb <= i__1 {
                i__2 = *n;
                i__ = 1 as ::core::ffi::c_int;
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
    mut c__: *mut ::core::ffi::c_double,
    mut d__: *mut ::core::ffi::c_double,
    mut e: *mut ::core::ffi::c_double,
    mut f: *mut ::core::ffi::c_double,
    mut g: *mut ::core::ffi::c_double,
    mut h__: *mut ::core::ffi::c_double,
    mut lc: *mut ::core::ffi::c_int,
    mut mc: *mut ::core::ffi::c_int,
    mut le: *mut ::core::ffi::c_int,
    mut me: *mut ::core::ffi::c_int,
    mut lg: *mut ::core::ffi::c_int,
    mut mg: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut x: *mut ::core::ffi::c_double,
    mut xnrm: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
    mut jw: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
) {
    let mut current_block: u64;
    let epmach: ::core::ffi::c_double = 2.22e-16f64;
    let mut c_dim1: ::core::ffi::c_int = 0;
    let mut c_offset: ::core::ffi::c_int = 0;
    let mut e_dim1: ::core::ffi::c_int = 0;
    let mut e_offset: ::core::ffi::c_int = 0;
    let mut g_dim1: ::core::ffi::c_int = 0;
    let mut g_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut i__3: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut k: ::core::ffi::c_int = 0;
    let mut l: ::core::ffi::c_int = 0;
    let mut t: ::core::ffi::c_double = 0.;
    let mut ie: ::core::ffi::c_int = 0;
    let mut if__: ::core::ffi::c_int = 0;
    let mut ig: ::core::ffi::c_int = 0;
    let mut iw: ::core::ffi::c_int = 0;
    let mut mc1: ::core::ffi::c_int = 0;
    let mut krank: ::core::ffi::c_int = 0;
    d__ = d__.offset(-1);
    f = f.offset(-1);
    h__ = h__.offset(-1);
    x = x.offset(-1);
    g_dim1 = *lg;
    g_offset = 1 as ::core::ffi::c_int + g_dim1;
    g = g.offset(-(g_offset as isize));
    e_dim1 = *le;
    e_offset = 1 as ::core::ffi::c_int + e_dim1;
    e = e.offset(-(e_offset as isize));
    c_dim1 = *lc;
    c_offset = 1 as ::core::ffi::c_int + c_dim1;
    c__ = c__.offset(-(c_offset as isize));
    w = w.offset(-1);
    jw = jw.offset(-1);
    *mode = 2 as ::core::ffi::c_int;
    if !(*mc > *n) {
        l = *n - *mc;
        mc1 = *mc + 1 as ::core::ffi::c_int;
        iw = (l + 1 as ::core::ffi::c_int) * (*mg + 2 as ::core::ffi::c_int)
            + (*mg << 1 as ::core::ffi::c_int)
            + *mc;
        ie = iw + *mc + 1 as ::core::ffi::c_int;
        if__ = ie + *me * l;
        ig = if__ + *me;
        i__1 = *mc;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            i__2 = i__ + 1 as ::core::ffi::c_int;
            j = if i__2 <= *lc { i__2 } else { *lc };
            i__2 = i__ + 1 as ::core::ffi::c_int;
            i__3 = *mc - i__;
            h12_(
                &raw const c__1,
                &raw mut i__,
                &raw mut i__2,
                n,
                c__.offset((i__ + c_dim1) as isize) as *mut ::core::ffi::c_double,
                lc,
                w.offset((iw + i__) as isize) as *mut ::core::ffi::c_double,
                c__.offset((j + c_dim1) as isize) as *mut ::core::ffi::c_double,
                lc,
                &raw const c__1,
                &raw mut i__3,
            );
            i__2 = i__ + 1 as ::core::ffi::c_int;
            h12_(
                &raw const c__2,
                &raw mut i__,
                &raw mut i__2,
                n,
                c__.offset((i__ + c_dim1) as isize) as *mut ::core::ffi::c_double,
                lc,
                w.offset((iw + i__) as isize) as *mut ::core::ffi::c_double,
                e.offset(e_offset as isize) as *mut ::core::ffi::c_double,
                le,
                &raw const c__1,
                me,
            );
            i__2 = i__ + 1 as ::core::ffi::c_int;
            h12_(
                &raw const c__2,
                &raw mut i__,
                &raw mut i__2,
                n,
                c__.offset((i__ + c_dim1) as isize) as *mut ::core::ffi::c_double,
                lc,
                w.offset((iw + i__) as isize) as *mut ::core::ffi::c_double,
                g.offset(g_offset as isize) as *mut ::core::ffi::c_double,
                lg,
                &raw const c__1,
                mg,
            );
            i__ += 1;
        }
        *mode = 6 as ::core::ffi::c_int;
        i__2 = *mc;
        i__ = 1 as ::core::ffi::c_int;
        loop {
            if !(i__ <= i__2) {
                current_block = 18377268871191777778;
                break;
            }
            d__1 = *c__.offset((i__ + i__ * c_dim1) as isize);
            if (d__1).abs() < epmach {
                current_block = 6183821358311025092;
                break;
            }
            i__1 = i__ - 1 as ::core::ffi::c_int;
            *x.offset(i__ as isize) = (*d__.offset(i__ as isize)
                - ddot_sl__(
                    &raw mut i__1,
                    c__.offset((i__ + c_dim1) as isize) as *mut ::core::ffi::c_double,
                    *lc,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                ))
                / *c__.offset((i__ + i__ * c_dim1) as isize);
            i__ += 1;
        }
        match current_block {
            6183821358311025092 => {}
            _ => {
                *mode = 1 as ::core::ffi::c_int;
                *w.offset(mc1 as isize) = 0.0f64;
                i__2 = *mg;
                dcopy___(
                    &raw mut i__2,
                    w.offset(mc1 as isize) as *mut ::core::ffi::c_double,
                    0 as ::core::ffi::c_int,
                    w.offset(mc1 as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                if *mc == *n {
                    current_block = 7514739150765462945;
                } else {
                    i__2 = *me;
                    i__ = 1 as ::core::ffi::c_int;
                    while i__ <= i__2 {
                        *w.offset((if__ - 1 as ::core::ffi::c_int + i__) as isize) = *f
                            .offset(i__ as isize)
                            - ddot_sl__(
                                mc,
                                e.offset((i__ + e_dim1) as isize) as *mut ::core::ffi::c_double,
                                *le,
                                x.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            );
                        i__ += 1;
                    }
                    i__2 = *me;
                    i__ = 1 as ::core::ffi::c_int;
                    while i__ <= i__2 {
                        dcopy___(
                            &raw mut l,
                            e.offset((i__ + mc1 * e_dim1) as isize) as *mut ::core::ffi::c_double,
                            *le,
                            w.offset((ie - 1 as ::core::ffi::c_int + i__) as isize)
                                as *mut ::core::ffi::c_double,
                            *me,
                        );
                        i__ += 1;
                    }
                    i__2 = *mg;
                    i__ = 1 as ::core::ffi::c_int;
                    while i__ <= i__2 {
                        dcopy___(
                            &raw mut l,
                            g.offset((i__ + mc1 * g_dim1) as isize) as *mut ::core::ffi::c_double,
                            *lg,
                            w.offset((ig - 1 as ::core::ffi::c_int + i__) as isize)
                                as *mut ::core::ffi::c_double,
                            *mg,
                        );
                        i__ += 1;
                    }
                    if *mg > 0 as ::core::ffi::c_int {
                        i__2 = *mg;
                        i__ = 1 as ::core::ffi::c_int;
                        while i__ <= i__2 {
                            *h__.offset(i__ as isize) -= ddot_sl__(
                                mc,
                                g.offset((i__ + g_dim1) as isize) as *mut ::core::ffi::c_double,
                                *lg,
                                x.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            );
                            i__ += 1;
                        }
                        lsi_(
                            w.offset(ie as isize) as *mut ::core::ffi::c_double,
                            w.offset(if__ as isize) as *mut ::core::ffi::c_double,
                            w.offset(ig as isize) as *mut ::core::ffi::c_double,
                            h__.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            me,
                            me,
                            mg,
                            mg,
                            &raw mut l,
                            x.offset(mc1 as isize) as *mut ::core::ffi::c_double,
                            xnrm,
                            w.offset(mc1 as isize) as *mut ::core::ffi::c_double,
                            jw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
                            mode,
                        );
                        if *mc == 0 as ::core::ffi::c_int {
                            current_block = 6183821358311025092;
                        } else {
                            t = dnrm2___(
                                mc,
                                x.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            );
                            *xnrm = (*xnrm * *xnrm + t * t).sqrt();
                            if *mode != 1 as ::core::ffi::c_int {
                                current_block = 6183821358311025092;
                            } else {
                                current_block = 7514739150765462945;
                            }
                        }
                    } else {
                        *mode = 7 as ::core::ffi::c_int;
                        k = if *le >= *n { *le } else { *n };
                        t = (epmach).sqrt();
                        hfti_(
                            w.offset(ie as isize) as *mut ::core::ffi::c_double,
                            me,
                            me,
                            &raw mut l,
                            w.offset(if__ as isize) as *mut ::core::ffi::c_double,
                            &raw mut k,
                            &raw const c__1,
                            &raw mut t,
                            &raw mut krank,
                            xnrm,
                            w.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            w.offset((l + 1 as ::core::ffi::c_int) as isize)
                                as *mut ::core::ffi::c_double,
                            jw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
                        );
                        dcopy___(
                            &raw mut l,
                            w.offset(if__ as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            x.offset(mc1 as isize) as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                        if krank != l {
                            current_block = 6183821358311025092;
                        } else {
                            *mode = 1 as ::core::ffi::c_int;
                            current_block = 7514739150765462945;
                        }
                    }
                }
                match current_block {
                    6183821358311025092 => {}
                    _ => {
                        i__2 = *me;
                        i__ = 1 as ::core::ffi::c_int;
                        while i__ <= i__2 {
                            *f.offset(i__ as isize) = ddot_sl__(
                                n,
                                e.offset((i__ + e_dim1) as isize) as *mut ::core::ffi::c_double,
                                *le,
                                x.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            ) - *f.offset(i__ as isize);
                            i__ += 1;
                        }
                        i__2 = *mc;
                        i__ = 1 as ::core::ffi::c_int;
                        while i__ <= i__2 {
                            *d__.offset(i__ as isize) = ddot_sl__(
                                me,
                                e.offset((i__ * e_dim1 + 1 as ::core::ffi::c_int) as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                                f.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            ) - ddot_sl__(
                                mg,
                                g.offset((i__ * g_dim1 + 1 as ::core::ffi::c_int) as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                                w.offset(mc1 as isize) as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            );
                            i__ += 1;
                        }
                        i__ = *mc;
                        while i__ >= 1 as ::core::ffi::c_int {
                            i__2 = i__ + 1 as ::core::ffi::c_int;
                            h12_(
                                &raw const c__2,
                                &raw mut i__,
                                &raw mut i__2,
                                n,
                                c__.offset((i__ + c_dim1) as isize) as *mut ::core::ffi::c_double,
                                lc,
                                w.offset((iw + i__) as isize) as *mut ::core::ffi::c_double,
                                x.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                &raw const c__1,
                                &raw const c__1,
                                &raw const c__1,
                            );
                            i__ -= 1;
                        }
                        i__ = *mc;
                        while i__ >= 1 as ::core::ffi::c_int {
                            i__2 = i__ + 1 as ::core::ffi::c_int;
                            j = if i__2 <= *lc { i__2 } else { *lc };
                            i__2 = *mc - i__;
                            *w.offset(i__ as isize) = (*d__.offset(i__ as isize)
                                - ddot_sl__(
                                    &raw mut i__2,
                                    c__.offset((j + i__ * c_dim1) as isize)
                                        as *mut ::core::ffi::c_double,
                                    1 as ::core::ffi::c_int,
                                    w.offset(j as isize) as *mut ::core::ffi::c_double,
                                    1 as ::core::ffi::c_int,
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
    mut m: *mut ::core::ffi::c_int,
    mut meq: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut nl: *mut ::core::ffi::c_int,
    mut la: *mut ::core::ffi::c_int,
    mut l: *mut ::core::ffi::c_double,
    mut g: *mut ::core::ffi::c_double,
    mut a: *mut ::core::ffi::c_double,
    mut b: *mut ::core::ffi::c_double,
    mut xl: *const ::core::ffi::c_double,
    mut xu: *const ::core::ffi::c_double,
    mut x: *mut ::core::ffi::c_double,
    mut y: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
    mut jw: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
) {
    let one: ::core::ffi::c_double = 1.0f64;
    let mut a_dim1: ::core::ffi::c_int = 0;
    let mut a_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut i1: ::core::ffi::c_int = 0;
    let mut i2: ::core::ffi::c_int = 0;
    let mut i3: ::core::ffi::c_int = 0;
    let mut i4: ::core::ffi::c_int = 0;
    let mut m1: ::core::ffi::c_int = 0;
    let mut n1: ::core::ffi::c_int = 0;
    let mut n2: ::core::ffi::c_int = 0;
    let mut n3: ::core::ffi::c_int = 0;
    let mut ic: ::core::ffi::c_int = 0;
    let mut id: ::core::ffi::c_int = 0;
    let mut ie: ::core::ffi::c_int = 0;
    let mut if__: ::core::ffi::c_int = 0;
    let mut ig: ::core::ffi::c_int = 0;
    let mut ih: ::core::ffi::c_int = 0;
    let mut il: ::core::ffi::c_int = 0;
    let mut im: ::core::ffi::c_int = 0;
    let mut ip: ::core::ffi::c_int = 0;
    let mut iu: ::core::ffi::c_int = 0;
    let mut iw: ::core::ffi::c_int = 0;
    let mut diag: ::core::ffi::c_double = 0.;
    let mut mineq: ::core::ffi::c_int = 0;
    let mut xnorm: ::core::ffi::c_double = 0.0f64;
    y = y.offset(-1);
    x = x.offset(-1);
    xu = xu.offset(-1);
    xl = xl.offset(-1);
    g = g.offset(-1);
    l = l.offset(-1);
    b = b.offset(-1);
    a_dim1 = *la;
    a_offset = 1 as ::core::ffi::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    w = w.offset(-1);
    jw = jw.offset(-1);
    n1 = *n + 1 as ::core::ffi::c_int;
    mineq = *m - *meq;
    m1 = mineq + *n + *n;
    n2 = n1 * *n / 2 as ::core::ffi::c_int + 1 as ::core::ffi::c_int;
    if n2 == *nl {
        n2 = 0 as ::core::ffi::c_int;
    } else {
        n2 = 1 as ::core::ffi::c_int;
    }
    n3 = *n - n2;
    i2 = 1 as ::core::ffi::c_int;
    i3 = 1 as ::core::ffi::c_int;
    i4 = 1 as ::core::ffi::c_int;
    ie = 1 as ::core::ffi::c_int;
    if__ = *n * *n + 1 as ::core::ffi::c_int;
    i__1 = n3;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= i__1 {
        i1 = n1 - i__;
        diag = (*l.offset(i2 as isize)).sqrt();
        *w.offset(i3 as isize) = 0.0f64;
        dcopy___(
            &raw mut i1,
            w.offset(i3 as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            w.offset(i3 as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        i__2 = i1 - n2;
        dcopy___(
            &raw mut i__2,
            l.offset(i2 as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            w.offset(i3 as isize) as *mut ::core::ffi::c_double,
            *n,
        );
        i__2 = i1 - n2;
        dscal_sl__(
            &raw mut i__2,
            &raw mut diag,
            w.offset(i3 as isize) as *mut ::core::ffi::c_double,
            *n,
        );
        *w.offset(i3 as isize) = diag;
        i__2 = i__ - 1 as ::core::ffi::c_int;
        *w.offset((if__ - 1 as ::core::ffi::c_int + i__) as isize) = (*g.offset(i__ as isize)
            - ddot_sl__(
                &raw mut i__2,
                w.offset(i4 as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
                w.offset(if__ as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
            ))
            / diag;
        i2 = i2 + i1 - n2;
        i3 += n1;
        i4 += *n;
        i__ += 1;
    }
    if n2 == 1 as ::core::ffi::c_int {
        *w.offset(i3 as isize) = *l.offset(*nl as isize);
        *w.offset(i4 as isize) = 0.0f64;
        dcopy___(
            &raw mut n3,
            w.offset(i4 as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            w.offset(i4 as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        *w.offset((if__ - 1 as ::core::ffi::c_int + *n) as isize) = 0.0f64;
    }
    d__1 = -one;
    dscal_sl__(
        n,
        &raw mut d__1,
        w.offset(if__ as isize) as *mut ::core::ffi::c_double,
        1 as ::core::ffi::c_int,
    );
    ic = if__ + *n;
    id = ic + *meq * *n;
    if *meq > 0 as ::core::ffi::c_int {
        i__1 = *meq;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            dcopy___(
                n,
                a.offset((i__ + a_dim1) as isize) as *mut ::core::ffi::c_double,
                *la,
                w.offset((ic - 1 as ::core::ffi::c_int + i__) as isize)
                    as *mut ::core::ffi::c_double,
                *meq,
            );
            i__ += 1;
        }
        dcopy___(
            meq,
            b.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            w.offset(id as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        d__1 = -one;
        dscal_sl__(
            meq,
            &raw mut d__1,
            w.offset(id as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
    }
    ig = id + *meq;
    if mineq > 0 as ::core::ffi::c_int {
        i__1 = mineq;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            dcopy___(
                n,
                a.offset((*meq + i__ + a_dim1) as isize) as *mut ::core::ffi::c_double,
                *la,
                w.offset((ig - 1 as ::core::ffi::c_int + i__) as isize)
                    as *mut ::core::ffi::c_double,
                m1,
            );
            i__ += 1;
        }
    }
    ip = ig + mineq;
    i__1 = *n;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= i__1 {
        *w.offset((ip - 1 as ::core::ffi::c_int + i__) as isize) = 0.0f64;
        dcopy___(
            n,
            w.offset((ip - 1 as ::core::ffi::c_int + i__) as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            w.offset((ip - 1 as ::core::ffi::c_int + i__) as isize) as *mut ::core::ffi::c_double,
            m1,
        );
        i__ += 1;
    }
    i__1 = m1 + 1 as ::core::ffi::c_int;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= *n {
        if nlopt_isinf(*xl.offset(i__ as isize)) == 0 {
            *w.offset((ip - i__1 + i__ * i__1) as isize) = 1.0f64;
        }
        i__ += 1;
    }
    im = ip + *n;
    i__1 = *n;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= i__1 {
        *w.offset((im - 1 as ::core::ffi::c_int + i__) as isize) = 0.0f64;
        dcopy___(
            n,
            w.offset((im - 1 as ::core::ffi::c_int + i__) as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            w.offset((im - 1 as ::core::ffi::c_int + i__) as isize) as *mut ::core::ffi::c_double,
            m1,
        );
        i__ += 1;
    }
    i__1 = m1 + 1 as ::core::ffi::c_int;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= *n {
        if nlopt_isinf(*xu.offset(i__ as isize)) == 0 {
            *w.offset((im - i__1 + i__ * i__1) as isize) = -1.0f64;
        }
        i__ += 1;
    }
    ih = ig + m1 * *n;
    if mineq > 0 as ::core::ffi::c_int {
        dcopy___(
            &raw mut mineq,
            b.offset((*meq + 1 as ::core::ffi::c_int) as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            w.offset(ih as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        d__1 = -one;
        dscal_sl__(
            &raw mut mineq,
            &raw mut d__1,
            w.offset(ih as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
    }
    il = ih + mineq;
    iu = il + *n;
    i__ = 1 as ::core::ffi::c_int;
    while i__ <= *n {
        *w.offset((il - 1 as ::core::ffi::c_int + i__) as isize) =
            if nlopt_isinf(*xl.offset(i__ as isize)) != 0 {
                0 as ::core::ffi::c_int as ::core::ffi::c_double
            } else {
                *xl.offset(i__ as isize)
            };
        *w.offset((iu - 1 as ::core::ffi::c_int + i__) as isize) =
            if nlopt_isinf(*xu.offset(i__ as isize)) != 0 {
                0 as ::core::ffi::c_int as ::core::ffi::c_double
            } else {
                -*xu.offset(i__ as isize)
            };
        i__ += 1;
    }
    iw = iu + *n;
    i__1 = if 1 as ::core::ffi::c_int >= *meq {
        1 as ::core::ffi::c_int
    } else {
        *meq
    };
    lsei_(
        w.offset(ic as isize) as *mut ::core::ffi::c_double,
        w.offset(id as isize) as *mut ::core::ffi::c_double,
        w.offset(ie as isize) as *mut ::core::ffi::c_double,
        w.offset(if__ as isize) as *mut ::core::ffi::c_double,
        w.offset(ig as isize) as *mut ::core::ffi::c_double,
        w.offset(ih as isize) as *mut ::core::ffi::c_double,
        &raw mut i__1,
        meq,
        n,
        n,
        &raw mut m1,
        &raw mut m1,
        n,
        x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
        &raw mut xnorm,
        w.offset(iw as isize) as *mut ::core::ffi::c_double,
        jw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
        mode,
    );
    if *mode == 1 as ::core::ffi::c_int {
        dcopy___(
            m,
            w.offset(iw as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            y.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        dcopy___(
            &raw mut n3,
            w.offset((iw + *m) as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            y.offset((*m + 1 as ::core::ffi::c_int) as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        dcopy___(
            &raw mut n3,
            w.offset((iw + *m + *n) as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            y.offset((*m + n3 + 1 as ::core::ffi::c_int) as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        i__1 = *n;
        i__ = 1 as ::core::ffi::c_int;
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
    mut n: *mut ::core::ffi::c_int,
    mut a: *mut ::core::ffi::c_double,
    mut z__: *mut ::core::ffi::c_double,
    mut sigma: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
) {
    let one: ::core::ffi::c_double = 1.0f64;
    let four: ::core::ffi::c_double = 4.0f64;
    let epmach: ::core::ffi::c_double = 2.22e-16f64;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut t: ::core::ffi::c_double = 0.;
    let mut u: ::core::ffi::c_double = 0.;
    let mut v: ::core::ffi::c_double = 0.;
    let mut ij: ::core::ffi::c_int = 0;
    let mut tp: ::core::ffi::c_double = 0.;
    let mut beta: ::core::ffi::c_double = 0.;
    let mut gamma_: ::core::ffi::c_double = 0.;
    let mut alpha: ::core::ffi::c_double = 0.;
    let mut delta: ::core::ffi::c_double = 0.;
    w = w.offset(-1);
    z__ = z__.offset(-1);
    a = a.offset(-1);
    if !(*sigma == 0.0f64) {
        ij = 1 as ::core::ffi::c_int;
        t = one / *sigma;
        if !(*sigma > 0.0f64) {
            i__1 = *n;
            i__ = 1 as ::core::ffi::c_int;
            while i__ <= i__1 {
                *w.offset(i__ as isize) = *z__.offset(i__ as isize);
                i__ += 1;
            }
            i__1 = *n;
            i__ = 1 as ::core::ffi::c_int;
            while i__ <= i__1 {
                v = *w.offset(i__ as isize);
                t += v * v / *a.offset(ij as isize);
                i__2 = *n;
                j = i__ + 1 as ::core::ffi::c_int;
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
            i__ = 1 as ::core::ffi::c_int;
            while i__ <= i__1 {
                j = *n + 1 as ::core::ffi::c_int - i__;
                ij -= i__;
                u = *w.offset(j as isize);
                *w.offset(j as isize) = t;
                t -= u * u / *a.offset(ij as isize);
                i__ += 1;
            }
        }
        i__1 = *n;
        i__ = 1 as ::core::ffi::c_int;
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
                j = i__ + 1 as ::core::ffi::c_int;
                while j <= i__2 {
                    ij += 1;
                    u = *a.offset(ij as isize);
                    *a.offset(ij as isize) = gamma_ * u + beta * *z__.offset(j as isize);
                    *z__.offset(j as isize) -= v * u;
                    j += 1;
                }
            } else {
                i__2 = *n;
                j = i__ + 1 as ::core::ffi::c_int;
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
    mut m: *mut ::core::ffi::c_int,
    mut meq: *mut ::core::ffi::c_int,
    mut la: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut x: *mut ::core::ffi::c_double,
    mut xl: *const ::core::ffi::c_double,
    mut xu: *const ::core::ffi::c_double,
    mut f: *mut ::core::ffi::c_double,
    mut c__: *mut ::core::ffi::c_double,
    mut g: *mut ::core::ffi::c_double,
    mut a: *mut ::core::ffi::c_double,
    mut acc: *mut ::core::ffi::c_double,
    mut iter: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
    mut r__: *mut ::core::ffi::c_double,
    mut l: *mut ::core::ffi::c_double,
    mut x0: *mut ::core::ffi::c_double,
    mut mu: *mut ::core::ffi::c_double,
    mut s: *mut ::core::ffi::c_double,
    mut u: *mut ::core::ffi::c_double,
    mut v: *mut ::core::ffi::c_double,
    mut w: *mut ::core::ffi::c_double,
    mut iw: *mut ::core::ffi::c_int,
    mut state: *mut slsqpb_state,
) {
    let mut current_block: u64;
    let one: ::core::ffi::c_double = 1.0f64;
    let alfmin: ::core::ffi::c_double = 0.1f64;
    let hun: ::core::ffi::c_double = 100.0f64;
    let ten: ::core::ffi::c_double = 10.0f64;
    let two: ::core::ffi::c_double = 2.0f64;
    let mut a_dim1: ::core::ffi::c_int = 0;
    let mut a_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut d__1: ::core::ffi::c_double = 0.;
    let mut d__2: ::core::ffi::c_double = 0.;
    let mut i__: ::core::ffi::c_int = 0;
    let mut j: ::core::ffi::c_int = 0;
    let mut k: ::core::ffi::c_int = 0;
    let mut t: ::core::ffi::c_double = 0.;
    let mut f0: ::core::ffi::c_double = 0.;
    let mut h1: ::core::ffi::c_double = 0.;
    let mut h2: ::core::ffi::c_double = 0.;
    let mut h3: ::core::ffi::c_double = 0.;
    let mut h4: ::core::ffi::c_double = 0.;
    let mut n1: ::core::ffi::c_int = 0;
    let mut n2: ::core::ffi::c_int = 0;
    let mut n3: ::core::ffi::c_int = 0;
    let mut t0: ::core::ffi::c_double = 0.;
    let mut gs: ::core::ffi::c_double = 0.;
    let mut tol: ::core::ffi::c_double = 0.;
    let mut line: ::core::ffi::c_int = 0;
    let mut alpha: ::core::ffi::c_double = 0.;
    let mut iexact: ::core::ffi::c_int = 0;
    let mut incons: ::core::ffi::c_int = 0;
    let mut ireset: ::core::ffi::c_int = 0;
    let mut itermx: ::core::ffi::c_int = 0;
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
    a_offset = 1 as ::core::ffi::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    g = g.offset(-1);
    xu = xu.offset(-1);
    xl = xl.offset(-1);
    x = x.offset(-1);
    w = w.offset(-1);
    iw = iw.offset(-1);
    if *mode == -(1 as ::core::ffi::c_int) {
        i__1 = *n;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            *u.offset(i__ as isize) = *g.offset(i__ as isize)
                - ddot_sl__(
                    m,
                    a.offset((i__ * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                        as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    r__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                )
                - *v.offset(i__ as isize);
            i__ += 1;
        }
        k = 0 as ::core::ffi::c_int;
        i__1 = *n;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            h1 = 0.0f64;
            k += 1;
            i__2 = *n;
            j = i__ + 1 as ::core::ffi::c_int;
            while j <= i__2 {
                k += 1;
                h1 += *l.offset(k as isize) * *s.offset(j as isize);
                j += 1;
            }
            *v.offset(i__ as isize) = *s.offset(i__ as isize) + h1;
            i__ += 1;
        }
        k = 1 as ::core::ffi::c_int;
        i__1 = *n;
        i__ = 1 as ::core::ffi::c_int;
        while i__ <= i__1 {
            *v.offset(i__ as isize) = *l.offset(k as isize) * *v.offset(i__ as isize);
            k = k + n1 - i__;
            i__ += 1;
        }
        i__ = *n;
        while i__ >= 1 as ::core::ffi::c_int {
            h1 = 0.0f64;
            k = i__;
            i__1 = i__ - 1 as ::core::ffi::c_int;
            j = 1 as ::core::ffi::c_int;
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
            s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        h2 = ddot_sl__(
            n,
            s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
            v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        h3 = h2 * 0.2f64;
        if h1 < h3 {
            h4 = (h2 - h3) / (h2 - h1);
            h1 = h3;
            dscal_sl__(
                n,
                &raw mut h4,
                u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
            );
            d__1 = one - h4;
            daxpy_sl__(
                n,
                &raw mut d__1,
                v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
                u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                1 as ::core::ffi::c_int,
            );
        }
        d__1 = one / h1;
        ldl_(
            n,
            l.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            &raw mut d__1,
            v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
        );
        d__1 = -one / h2;
        ldl_(
            n,
            l.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            &raw mut d__1,
            u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
        );
        current_block = 17565015850456805979;
    } else if *mode == 0 as ::core::ffi::c_int {
        itermx = *iter;
        if *acc >= 0.0f64 {
            iexact = 0 as ::core::ffi::c_int;
        } else {
            iexact = 1 as ::core::ffi::c_int;
        }
        *acc = (*acc).abs();
        tol = ten * *acc;
        *iter = 0 as ::core::ffi::c_int;
        ireset = 0 as ::core::ffi::c_int;
        n1 = *n + 1 as ::core::ffi::c_int;
        n2 = n1 * *n / 2 as ::core::ffi::c_int;
        n3 = n2 + 1 as ::core::ffi::c_int;
        *s.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
        *mu.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
        dcopy___(
            n,
            s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        dcopy___(
            m,
            mu.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            0 as ::core::ffi::c_int,
            mu.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
            1 as ::core::ffi::c_int,
        );
        current_block = 14140784183277947939;
    } else {
        t = *f;
        i__1 = *m;
        j = 1 as ::core::ffi::c_int;
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
        match iexact + 1 as ::core::ffi::c_int {
            1 => {
                current_block = 15522146559271336860;
                match current_block {
                    15522146559271336860 => {
                        if nlopt_isfinite(h1) != 0 {
                            if h1 <= h3 / ten || line > 10 as ::core::ffi::c_int {
                                current_block = 2101818066333513453;
                            } else {
                                d__1 = h3 / (two * (h3 - h1));
                                alpha = if d__1 >= alfmin { d__1 } else { alfmin };
                                current_block = 11187960484546229445;
                            }
                        } else {
                            alpha = if alpha * 0.5f64 >= alfmin {
                                alpha * 0.5f64
                            } else {
                                alfmin
                            };
                            current_block = 11187960484546229445;
                        }
                    }
                    _ => {}
                }
                match current_block {
                    11187960484546229445 => {}
                    _ => {
                        h3 = 0.0f64;
                        i__1 = *m;
                        j = 1 as ::core::ffi::c_int;
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
                                s.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            ) < *acc)
                            && h3 < *acc
                        {
                            *mode = 0 as ::core::ffi::c_int;
                        } else {
                            *mode = -(1 as ::core::ffi::c_int);
                        }
                        current_block = 7663700016002789541;
                    }
                }
            }
            2 => {
                current_block = 15447661902298937816;
            }
            _ => {
                current_block = 2101818066333513453;
                match current_block {
                    15522146559271336860 => {
                        if nlopt_isfinite(h1) != 0 {
                            if h1 <= h3 / ten || line > 10 as ::core::ffi::c_int {
                                current_block = 2101818066333513453;
                            } else {
                                d__1 = h3 / (two * (h3 - h1));
                                alpha = if d__1 >= alfmin { d__1 } else { alfmin };
                                current_block = 11187960484546229445;
                            }
                        } else {
                            alpha = if alpha * 0.5f64 >= alfmin {
                                alpha * 0.5f64
                            } else {
                                alfmin
                            };
                            current_block = 11187960484546229445;
                        }
                    }
                    _ => {}
                }
                match current_block {
                    11187960484546229445 => {}
                    _ => {
                        h3 = 0.0f64;
                        i__1 = *m;
                        j = 1 as ::core::ffi::c_int;
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
                                s.offset(1 as ::core::ffi::c_int as isize)
                                    as *mut ::core::ffi::c_double,
                                1 as ::core::ffi::c_int,
                            ) < *acc)
                            && h3 < *acc
                        {
                            *mode = 0 as ::core::ffi::c_int;
                        } else {
                            *mode = -(1 as ::core::ffi::c_int);
                        }
                        current_block = 7663700016002789541;
                    }
                }
            }
        }
    }
    '_L330: loop {
        match current_block {
            15447661902298937816 => {
                *mode = 9 as ::core::ffi::c_int;
                return;
            }
            11187960484546229445 => {
                line += 1;
                h3 = alpha * h3;
                dscal_sl__(
                    n,
                    &raw mut alpha,
                    s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                dcopy___(
                    n,
                    x0.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                daxpy_sl__(
                    n,
                    &raw const one,
                    s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                i__1 = *n;
                i__ = 1 as ::core::ffi::c_int;
                while i__ <= i__1 {
                    if *x.offset(i__ as isize) < *xl.offset(i__ as isize) {
                        *x.offset(i__ as isize) = *xl.offset(i__ as isize);
                    } else if *x.offset(i__ as isize) > *xu.offset(i__ as isize) {
                        *x.offset(i__ as isize) = *xu.offset(i__ as isize);
                    }
                    i__ += 1;
                }
                *mode = if line == 1 as ::core::ffi::c_int {
                    -(2 as ::core::ffi::c_int)
                } else {
                    1 as ::core::ffi::c_int
                };
                current_block = 7663700016002789541;
            }
            17565015850456805979 => {
                *iter += 1;
                *mode = 9 as ::core::ffi::c_int;
                if *iter > itermx && itermx > 0 as ::core::ffi::c_int {
                    current_block = 7663700016002789541;
                    continue;
                }
                dcopy___(
                    n,
                    xl.offset(1 as ::core::ffi::c_int as isize) as *const ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                dcopy___(
                    n,
                    xu.offset(1 as ::core::ffi::c_int as isize) as *const ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                d__1 = -one;
                daxpy_sl__(
                    n,
                    &raw mut d__1,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                d__1 = -one;
                daxpy_sl__(
                    n,
                    &raw mut d__1,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                h4 = one;
                lsq_(
                    m,
                    meq,
                    n,
                    &raw mut n3,
                    la,
                    l.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    g.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    a.offset(a_offset as isize) as *mut ::core::ffi::c_double,
                    c__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    u.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    v.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    r__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    w.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    iw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
                    mode,
                );
                if *mode == 6 as ::core::ffi::c_int {
                    if *n == *meq {
                        *mode = 4 as ::core::ffi::c_int;
                    }
                }
                if *mode == 4 as ::core::ffi::c_int {
                    i__1 = *m;
                    j = 1 as ::core::ffi::c_int;
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
                    *s.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
                    dcopy___(
                        n,
                        s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        0 as ::core::ffi::c_int,
                        s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        1 as ::core::ffi::c_int,
                    );
                    h3 = 0.0f64;
                    *g.offset(n1 as isize) = 0.0f64;
                    *l.offset(n3 as isize) = hun;
                    *s.offset(n1 as isize) = one;
                    *u.offset(n1 as isize) = 0.0f64;
                    *v.offset(n1 as isize) = one;
                    incons = 0 as ::core::ffi::c_int;
                    loop {
                        lsq_(
                            m,
                            meq,
                            &raw mut n1,
                            &raw mut n3,
                            la,
                            l.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            g.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            a.offset(a_offset as isize) as *mut ::core::ffi::c_double,
                            c__.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            u.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            v.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            s.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            r__.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            w.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            iw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
                            mode,
                        );
                        h4 = one - *s.offset(n1 as isize);
                        if !(*mode == 4 as ::core::ffi::c_int) {
                            break;
                        }
                        *l.offset(n3 as isize) = ten * *l.offset(n3 as isize);
                        incons += 1;
                        if incons > 5 as ::core::ffi::c_int {
                            current_block = 7663700016002789541;
                            continue '_L330;
                        }
                    }
                    if *mode != 1 as ::core::ffi::c_int {
                        current_block = 7663700016002789541;
                        continue;
                    }
                } else if *mode != 1 as ::core::ffi::c_int {
                    current_block = 7663700016002789541;
                    continue;
                }
                i__1 = *n;
                i__ = 1 as ::core::ffi::c_int;
                while i__ <= i__1 {
                    *v.offset(i__ as isize) = *g.offset(i__ as isize)
                        - ddot_sl__(
                            m,
                            a.offset((i__ * a_dim1 + 1 as ::core::ffi::c_int) as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                            r__.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        );
                    i__ += 1;
                }
                f0 = *f;
                dcopy___(
                    n,
                    x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    x0.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                gs = ddot_sl__(
                    n,
                    g.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                    s.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                    1 as ::core::ffi::c_int,
                );
                h1 = (gs).abs();
                h2 = 0.0f64;
                i__1 = *m;
                j = 1 as ::core::ffi::c_int;
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
                *mode = 0 as ::core::ffi::c_int;
                if h1 < *acc && h2 < *acc {
                    current_block = 7663700016002789541;
                    continue;
                }
                h1 = 0.0f64;
                i__1 = *m;
                j = 1 as ::core::ffi::c_int;
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
                *mode = 8 as ::core::ffi::c_int;
                if h3 >= 0.0f64 {
                    current_block = 14140784183277947939;
                    continue;
                }
                line = 0 as ::core::ffi::c_int;
                alpha = one;
                if iexact == 1 as ::core::ffi::c_int {
                    current_block = 15447661902298937816;
                } else {
                    current_block = 11187960484546229445;
                }
            }
            14140784183277947939 => {
                ireset += 1;
                if ireset > 5 as ::core::ffi::c_int {
                    d__1 = *f - f0;
                    if ((d__1).abs() < tol
                        || dnrm2___(
                            n,
                            s.offset(1 as ::core::ffi::c_int as isize)
                                as *mut ::core::ffi::c_double,
                            1 as ::core::ffi::c_int,
                        ) < tol)
                        && h3 < tol
                    {
                        *mode = 0 as ::core::ffi::c_int;
                    } else {
                        *mode = 8 as ::core::ffi::c_int;
                    }
                    current_block = 7663700016002789541;
                } else {
                    *l.offset(1 as ::core::ffi::c_int as isize) = 0.0f64;
                    dcopy___(
                        &raw mut n2,
                        l.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        0 as ::core::ffi::c_int,
                        l.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
                        1 as ::core::ffi::c_int,
                    );
                    j = 1 as ::core::ffi::c_int;
                    i__1 = *n;
                    i__ = 1 as ::core::ffi::c_int;
                    while i__ <= i__1 {
                        *l.offset(j as isize) = one;
                        j = j + n1 - i__;
                        i__ += 1;
                    }
                    current_block = 17565015850456805979;
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
    mut m: *mut ::core::ffi::c_int,
    mut meq: *mut ::core::ffi::c_int,
    mut la: *mut ::core::ffi::c_int,
    mut n: *mut ::core::ffi::c_int,
    mut x: *mut ::core::ffi::c_double,
    mut xl: *const ::core::ffi::c_double,
    mut xu: *const ::core::ffi::c_double,
    mut f: *mut ::core::ffi::c_double,
    mut c__: *mut ::core::ffi::c_double,
    mut g: *mut ::core::ffi::c_double,
    mut a: *mut ::core::ffi::c_double,
    mut acc: *mut ::core::ffi::c_double,
    mut iter: *mut ::core::ffi::c_int,
    mut mode: *mut ::core::ffi::c_int,
    mut w: *mut ::core::ffi::c_double,
    mut l_w__: *mut ::core::ffi::c_int,
    mut jw: *mut ::core::ffi::c_int,
    mut l_jw__: *mut ::core::ffi::c_int,
    mut state: *mut slsqpb_state,
) {
    let mut a_dim1: ::core::ffi::c_int = 0;
    let mut a_offset: ::core::ffi::c_int = 0;
    let mut i__1: ::core::ffi::c_int = 0;
    let mut i__2: ::core::ffi::c_int = 0;
    let mut n1: ::core::ffi::c_int = 0;
    let mut il: ::core::ffi::c_int = 0;
    let mut im: ::core::ffi::c_int = 0;
    let mut ir: ::core::ffi::c_int = 0;
    let mut is: ::core::ffi::c_int = 0;
    let mut iu: ::core::ffi::c_int = 0;
    let mut iv: ::core::ffi::c_int = 0;
    let mut iw: ::core::ffi::c_int = 0;
    let mut ix: ::core::ffi::c_int = 0;
    let mut mineq: ::core::ffi::c_int = 0;
    c__ = c__.offset(-1);
    a_dim1 = *la;
    a_offset = 1 as ::core::ffi::c_int + a_dim1;
    a = a.offset(-(a_offset as isize));
    g = g.offset(-1);
    xu = xu.offset(-1);
    xl = xl.offset(-1);
    x = x.offset(-1);
    w = w.offset(-1);
    jw = jw.offset(-1);
    n1 = *n + 1 as ::core::ffi::c_int;
    mineq = *m - *meq + n1 + n1;
    il = (n1 * 3 as ::core::ffi::c_int + *m) * (n1 + 1 as ::core::ffi::c_int)
        + (n1 - *meq + 1 as ::core::ffi::c_int) * (mineq + 2 as ::core::ffi::c_int)
        + (mineq << 1 as ::core::ffi::c_int)
        + (n1 + mineq) * (n1 - *meq)
        + (*meq << 1 as ::core::ffi::c_int)
        + n1 * *n / 2 as ::core::ffi::c_int
        + (*m << 1 as ::core::ffi::c_int)
        + *n * 3 as ::core::ffi::c_int
        + (n1 << 2 as ::core::ffi::c_int)
        + 1 as ::core::ffi::c_int;
    i__1 = mineq;
    i__2 = n1 - *meq;
    im = if i__1 >= i__2 { i__1 } else { i__2 };
    if *l_w__ < il || *l_jw__ < im {
        *mode = (if 10 as ::core::ffi::c_int >= il {
            10 as ::core::ffi::c_int
        } else {
            il
        }) * 1000 as ::core::ffi::c_int;
        *mode += if 10 as ::core::ffi::c_int >= im {
            10 as ::core::ffi::c_int
        } else {
            im
        };
        return;
    }
    im = 1 as ::core::ffi::c_int;
    il = im
        + (if 1 as ::core::ffi::c_int >= *m {
            1 as ::core::ffi::c_int
        } else {
            *m
        });
    il = im + *la;
    ix = il + n1 * *n / 2 as ::core::ffi::c_int + 1 as ::core::ffi::c_int;
    ir = ix + *n;
    is = ir
        + *n
        + *n
        + (if 1 as ::core::ffi::c_int >= *m {
            1 as ::core::ffi::c_int
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
        x.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
        xl.offset(1 as ::core::ffi::c_int as isize) as *const ::core::ffi::c_double,
        xu.offset(1 as ::core::ffi::c_int as isize) as *const ::core::ffi::c_double,
        f,
        c__.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
        g.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_double,
        a.offset(a_offset as isize) as *mut ::core::ffi::c_double,
        acc,
        iter,
        mode,
        w.offset(ir as isize) as *mut ::core::ffi::c_double,
        w.offset(il as isize) as *mut ::core::ffi::c_double,
        w.offset(ix as isize) as *mut ::core::ffi::c_double,
        w.offset(im as isize) as *mut ::core::ffi::c_double,
        w.offset(is as isize) as *mut ::core::ffi::c_double,
        w.offset(iu as isize) as *mut ::core::ffi::c_double,
        w.offset(iv as isize) as *mut ::core::ffi::c_double,
        w.offset(iw as isize) as *mut ::core::ffi::c_double,
        jw.offset(1 as ::core::ffi::c_int as isize) as *mut ::core::ffi::c_int,
        state,
    );
    (*state).x0 = w.offset(ix as isize) as *mut ::core::ffi::c_double;
}
unsafe fn length_work(
    mut LEN_W: *mut ::core::ffi::c_int,
    mut LEN_JW: *mut ::core::ffi::c_int,
    mut M: ::core::ffi::c_int,
    mut MEQ: ::core::ffi::c_int,
    mut N: ::core::ffi::c_int,
) {
    let mut N1: ::core::ffi::c_int = N + 1 as ::core::ffi::c_int;
    let mut MINEQ: ::core::ffi::c_int = M - MEQ + N1 + N1;
    *LEN_W = (3 as ::core::ffi::c_int * N1 + M) * (N1 + 1 as ::core::ffi::c_int)
        + (N1 - MEQ + 1 as ::core::ffi::c_int) * (MINEQ + 2 as ::core::ffi::c_int)
        + 2 as ::core::ffi::c_int * MINEQ
        + (N1 + MINEQ) * (N1 - MEQ)
        + 2 as ::core::ffi::c_int * MEQ
        + N1
        + (N + 1 as ::core::ffi::c_int) * N / 2 as ::core::ffi::c_int
        + 2 as ::core::ffi::c_int * M
        + 3 as ::core::ffi::c_int * N
        + 3 as ::core::ffi::c_int * N1
        + 1 as ::core::ffi::c_int;
    *LEN_JW = MINEQ;
}

pub(crate) unsafe fn nlopt_slsqp<U>(
    mut n: ::core::ffi::c_uint,
    mut f: nlopt_func,
    mut f_data: *mut ::core::ffi::c_void,
    mut m: ::core::ffi::c_uint,
    mut fc: *mut nlopt_constraint,
    mut p: ::core::ffi::c_uint,
    mut h: *mut nlopt_constraint,
    mut lb: *const ::core::ffi::c_double,
    mut ub: *const ::core::ffi::c_double,
    mut x: *mut ::core::ffi::c_double,
    mut minf: *mut ::core::ffi::c_double,
    mut stop: *mut nlopt_stopping,
) -> nlopt_result {
    let mut current_block: u64;
    let mut state: slsqpb_state = slsqpb_state {
        t: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        f0: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        h1: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        h2: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        h3: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        h4: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        n1: 0 as ::core::ffi::c_int,
        n2: 0 as ::core::ffi::c_int,
        n3: 0 as ::core::ffi::c_int,
        t0: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        gs: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        tol: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        line: 0 as ::core::ffi::c_int,
        alpha: 0 as ::core::ffi::c_int as ::core::ffi::c_double,
        iexact: 0 as ::core::ffi::c_int,
        incons: 0 as ::core::ffi::c_int,
        ireset: 0 as ::core::ffi::c_int,
        itermx: 0 as ::core::ffi::c_int,
        x0: ::core::ptr::null_mut::<::core::ffi::c_double>(),
    };
    let mut mtot: ::core::ffi::c_uint = nlopt_count_constraints(m, fc);
    let mut ptot: ::core::ffi::c_uint = nlopt_count_constraints(p, h);
    let mut work: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut cgrad: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut c: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut grad: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut w: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut fcur: ::core::ffi::c_double = 0.;
    let mut xcur: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut fprev: ::core::ffi::c_double = 0.;
    let mut xprev: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut cgradtmp: *mut ::core::ffi::c_double = ::core::ptr::null_mut::<::core::ffi::c_double>();
    let mut mpi: ::core::ffi::c_int = mtot.wrapping_add(ptot) as ::core::ffi::c_int;
    let mut pi: ::core::ffi::c_int = ptot as ::core::ffi::c_int;
    let mut ni: ::core::ffi::c_int = n as ::core::ffi::c_int;
    let mut mpi1: ::core::ffi::c_int = if mpi > 0 as ::core::ffi::c_int {
        mpi
    } else {
        1 as ::core::ffi::c_int
    };
    let mut len_w: ::core::ffi::c_int = 0;
    let mut len_jw: ::core::ffi::c_int = 0;
    let mut jw: *mut ::core::ffi::c_int = ::core::ptr::null_mut::<::core::ffi::c_int>();
    let mut mode: ::core::ffi::c_int = 0 as ::core::ffi::c_int;
    let mut prev_mode: ::core::ffi::c_int = 0 as ::core::ffi::c_int;
    let mut acc: ::core::ffi::c_double = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
    let mut iter: ::core::ffi::c_int = 0 as ::core::ffi::c_int;
    let mut i: ::core::ffi::c_uint = 0;
    let mut ii: ::core::ffi::c_uint = 0;
    let mut ret: nlopt_result = NLOPT_SUCCESS;
    let mut feasible: ::core::ffi::c_int = 0;
    let mut feasible_cur: ::core::ffi::c_int = 0;
    let mut infeasibility: ::core::ffi::c_double = ::core::f64::INFINITY;
    let mut infeasibility_cur: ::core::ffi::c_double = ::core::f64::INFINITY;
    let mut max_cdim: ::core::ffi::c_uint = 0;
    let mut want_grad: ::core::ffi::c_int = 1 as ::core::ffi::c_int;
    if ptot > n {
        nlopt_stop_msg(
            stop,
            "slsqp: more equality constraints than variables"
        );
        return NLOPT_INVALID_ARGS;
    }
    max_cdim = if nlopt_max_constraint_dim(m, fc) >= nlopt_max_constraint_dim(p, h) {
        nlopt_max_constraint_dim(m, fc)
    } else {
        nlopt_max_constraint_dim(p, h)
    };
    length_work(&raw mut len_w, &raw mut len_jw, mpi, pi, ni);
    // work = malloc(
    //     (::std::mem::size_of::<::core::ffi::c_double>() as ::core::ffi::c_ulong)
    //         .wrapping_mul(
    //             (mpi1 as ::core::ffi::c_uint)
    //                 .wrapping_mul(n.wrapping_add(1 as ::core::ffi::c_int as ::core::ffi::c_uint))
    //                 .wrapping_add(mpi as ::core::ffi::c_uint)
    //                 .wrapping_add(n)
    //                 .wrapping_add(1 as ::core::ffi::c_int as ::core::ffi::c_uint)
    //                 .wrapping_add(n)
    //                 .wrapping_add(n)
    //                 .wrapping_add(max_cdim.wrapping_mul(n))
    //                 .wrapping_add(len_w as ::core::ffi::c_uint) as ::core::ffi::c_ulong,
    //         )
    //         .wrapping_add(
    //             (::std::mem::size_of::<::core::ffi::c_int>() as ::core::ffi::c_ulong)
    //                 .wrapping_mul(len_jw as ::core::ffi::c_uint as ::core::ffi::c_ulong),
    //         ),
    // ) as *mut ::core::ffi::c_double;

    let space_size = (mpi1 as ::core::ffi::c_uint)
        .wrapping_mul(n.wrapping_add(1 as ::core::ffi::c_int as ::core::ffi::c_uint))
        .wrapping_add(mpi as ::core::ffi::c_uint)
        .wrapping_add(n)
        .wrapping_add(1 as ::core::ffi::c_int as ::core::ffi::c_uint)
        .wrapping_add(n)
        .wrapping_add(n)
        .wrapping_add(max_cdim.wrapping_mul(n))
        .wrapping_add(len_w as ::core::ffi::c_uint)
        .wrapping_add(len_jw as ::core::ffi::c_uint); // size_of(::core::ffi::c_double) > size_of(::core::ffi::c_int)
    let mut space: Vec<::core::ffi::c_double> = vec![0.; usize::try_from(space_size).unwrap()];
    work = space.as_mut_ptr();

    if work.is_null() {
        return NLOPT_OUT_OF_MEMORY;
    }
    cgrad = work;
    c = cgrad.offset(
        (mpi1 as ::core::ffi::c_uint).wrapping_mul(n.wrapping_add(1 as ::core::ffi::c_uint))
            as isize,
    );
    grad = c.offset(mpi as isize);
    xcur = grad
        .offset(n as isize)
        .offset(1 as ::core::ffi::c_int as isize);
    xprev = xcur.offset(n as isize);
    cgradtmp = xprev.offset(n as isize);
    w = cgradtmp.offset(max_cdim.wrapping_mul(n) as isize);
    jw = w.offset(len_w as isize) as *mut ::core::ffi::c_int;
    memcpy(
        xcur as *mut ::core::ffi::c_void,
        x as *const ::core::ffi::c_void,
        (::core::mem::size_of::<::core::ffi::c_double>() as size_t).wrapping_mul(n as size_t),
    );
    memcpy(
        xprev as *mut ::core::ffi::c_void,
        x as *const ::core::ffi::c_void,
        (::core::mem::size_of::<::core::ffi::c_double>() as size_t).wrapping_mul(n as size_t),
    );
    *minf = ::core::f64::INFINITY;
    fcur = *minf;
    fprev = fcur;
    feasible_cur = 0 as ::core::ffi::c_int;
    feasible = feasible_cur;
    '_eval_f_and_grad: loop {
        want_grad = 1 as ::core::ffi::c_int;
        's_122: loop {
            let mut newgrad: *mut ::core::ffi::c_double =
                ::core::ptr::null_mut::<::core::ffi::c_double>();
            let mut newcgrad: *mut ::core::ffi::c_double =
                ::core::ptr::null_mut::<::core::ffi::c_double>();
            if want_grad != 0 {
                newgrad = grad;
                newcgrad = cgradtmp;
            }
            feasible_cur = 1 as ::core::ffi::c_int;
            infeasibility_cur = 0 as ::core::ffi::c_int as ::core::ffi::c_double;
            fcur = f.expect("non-null function pointer")(n, xcur, newgrad, f_data);
            *(*stop).nevals_p += 1;
            if nlopt_stop_forced(stop) != 0 {
                fcur = ::core::f64::INFINITY;
                ret = NLOPT_FORCED_STOP;
                break '_eval_f_and_grad;
            } else {
                if nlopt_isfinite(fcur) != 0 {
                    want_grad = 0 as ::core::ffi::c_int;
                    ii = 0 as ::core::ffi::c_uint;
                    i = 0 as ::core::ffi::c_uint;
                    while i < p {
                        let mut j: ::core::ffi::c_uint = 0;
                        let mut k: ::core::ffi::c_uint = 0;
                        nlopt_eval_constraint::<U>(
                            c.offset(ii as isize),
                            newcgrad,
                            h.offset(i as isize),
                            n,
                            xcur,
                        );
                        if nlopt_stop_forced(stop) != 0 {
                            ret = NLOPT_FORCED_STOP;
                            break '_eval_f_and_grad;
                        } else {
                            k = 0 as ::core::ffi::c_uint;
                            while k < (*h.offset(i as isize)).m {
                                infeasibility_cur =
                                    if infeasibility_cur >= (*c.offset(ii as isize)).abs() {
                                        infeasibility_cur
                                    } else {
                                        (*c.offset(ii as isize)).abs()
                                    };
                                feasible_cur = (feasible_cur != 0
                                    && (*c.offset(ii as isize)).abs()
                                        <= *(*h.offset(i as isize)).tol.offset(k as isize))
                                    as ::core::ffi::c_int;
                                if !newcgrad.is_null() {
                                    j = 0 as ::core::ffi::c_uint;
                                    while j < n {
                                        *cgrad.offset(
                                            j.wrapping_mul(mpi1 as ::core::ffi::c_uint)
                                                .wrapping_add(ii)
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
                    i = 0 as ::core::ffi::c_uint;
                    while i < m {
                        let mut j_0: ::core::ffi::c_uint = 0;
                        let mut k_0: ::core::ffi::c_uint = 0;
                        nlopt_eval_constraint::<U>(
                            c.offset(ii as isize),
                            newcgrad,
                            fc.offset(i as isize),
                            n,
                            xcur,
                        );
                        if nlopt_stop_forced(stop) != 0 {
                            ret = NLOPT_FORCED_STOP;
                            break '_eval_f_and_grad;
                        } else {
                            k_0 = 0 as ::core::ffi::c_uint;
                            while k_0 < (*fc.offset(i as isize)).m {
                                infeasibility_cur = if infeasibility_cur >= *c.offset(ii as isize) {
                                    infeasibility_cur
                                } else {
                                    *c.offset(ii as isize)
                                };
                                feasible_cur = (feasible_cur != 0
                                    && *c.offset(ii as isize)
                                        <= *(*fc.offset(i as isize)).tol.offset(k_0 as isize))
                                    as ::core::ffi::c_int;
                                if !newcgrad.is_null() {
                                    j_0 = 0 as ::core::ffi::c_uint;
                                    while j_0 < n {
                                        *cgrad.offset(
                                            j_0.wrapping_mul(mpi1 as ::core::ffi::c_uint)
                                                .wrapping_add(ii)
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
                            x as *mut ::core::ffi::c_void,
                            xcur as *const ::core::ffi::c_void,
                            (::core::mem::size_of::<::core::ffi::c_double>() as size_t)
                                .wrapping_mul(n as size_t),
                        );
                    }
                    if mode == -(1 as ::core::ffi::c_int) {
                        if nlopt_isinf(fprev) == 0 && feasible_cur != 0 {
                            if nlopt_stop_ftol(stop, fcur, fprev) != 0 {
                                ret = NLOPT_FTOL_REACHED;
                            } else if nlopt_stop_x(stop, xcur, xprev) != 0 {
                                ret = NLOPT_XTOL_REACHED;
                            }
                        }
                        fprev = fcur;
                        memcpy(
                            xprev as *mut ::core::ffi::c_void,
                            xcur as *const ::core::ffi::c_void,
                            (::core::mem::size_of::<::core::ffi::c_double>() as size_t)
                                .wrapping_mul(n as size_t),
                        );
                    }
                    if nlopt_stop_evals(stop) != 0 {
                        ret = NLOPT_MAXEVAL_REACHED;
                    } else if nlopt_stop_time(stop) != 0 {
                        ret = NLOPT_MAXTIME_REACHED;
                    } else if feasible != 0 && *minf < (*stop).minf_max {
                        ret = NLOPT_STOPVAL_REACHED;
                    }
                    if !(ret as ::core::ffi::c_int == NLOPT_SUCCESS as ::core::ffi::c_int) {
                        break '_eval_f_and_grad;
                    }
                    slsqp(
                        &raw mut mpi,
                        &raw mut pi,
                        &raw mut mpi1,
                        &raw mut ni,
                        xcur,
                        lb,
                        ub,
                        &raw mut fcur,
                        c,
                        grad,
                        cgrad,
                        &raw mut acc,
                        &raw mut iter,
                        &raw mut mode,
                        w,
                        &raw mut len_w,
                        jw,
                        &raw mut len_jw,
                        &raw mut state,
                    );
                    match mode {
                        -1 => {
                            if !(prev_mode == -(2 as ::core::ffi::c_int) && want_grad == 0) {
                                continue '_eval_f_and_grad;
                            }
                        }
                        -2 => {
                            continue '_eval_f_and_grad;
                        }
                        1 => {
                            break;
                        }
                        0 => {
                            break '_eval_f_and_grad;
                        }
                        8 => {
                            ret = NLOPT_ROUNDOFF_LIMITED;
                            if feasible_cur != 0 {
                                let mut save_ftol_rel: ::core::ffi::c_double = (*stop).ftol_rel;
                                let mut save_xtol_rel: ::core::ffi::c_double = (*stop).xtol_rel;
                                let mut save_ftol_abs: ::core::ffi::c_double = (*stop).ftol_abs;
                                (*stop).ftol_rel *=
                                    10 as ::core::ffi::c_int as ::core::ffi::c_double;
                                (*stop).ftol_abs *=
                                    10 as ::core::ffi::c_int as ::core::ffi::c_double;
                                (*stop).xtol_rel *=
                                    10 as ::core::ffi::c_int as ::core::ffi::c_double;
                                if nlopt_stop_ftol(stop, fcur, state.f0) != 0 {
                                    ret = NLOPT_FTOL_REACHED;
                                } else if nlopt_stop_x(stop, xcur, state.x0) != 0 {
                                    ret = NLOPT_XTOL_REACHED;
                                }
                                (*stop).ftol_rel = save_ftol_rel;
                                (*stop).ftol_abs = save_ftol_abs;
                                (*stop).xtol_rel = save_xtol_rel;
                            }
                            break '_eval_f_and_grad;
                        }
                        5 => {
                            current_block = 10253525168500281338;
                            break 's_122;
                        }
                        6 | 7 => {
                            current_block = 10253525168500281338;
                            break 's_122;
                        }
                        4 => {
                            current_block = 17546449896862090116;
                            break 's_122;
                        }
                        3 | 9 => {
                            current_block = 17546449896862090116;
                            break 's_122;
                        }
                        2 | _ => {
                            nlopt_stop_msg(
                                stop,
                                "bug: workspace is too small"
                            );
                            ret = NLOPT_INVALID_ARGS;
                            break '_eval_f_and_grad;
                        }
                    }
                }
            }
        }
        match current_block {
            10253525168500281338 => {
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
                x as *mut ::core::ffi::c_void,
                xprev as *const ::core::ffi::c_void,
                (::core::mem::size_of::<::core::ffi::c_double>() as size_t)
                    .wrapping_mul(n as size_t),
            );
        } else {
            *minf = fcur;
            memcpy(
                x as *mut ::core::ffi::c_void,
                xcur as *const ::core::ffi::c_void,
                (::core::mem::size_of::<::core::ffi::c_double>() as size_t)
                    .wrapping_mul(n as size_t),
            );
        }
    }
    // `space` (the buffer behind `work`) is dropped here.
    return ret;
}
