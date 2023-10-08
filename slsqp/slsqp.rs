#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![register_tool(c2rust)]
#![feature(register_tool)]
extern "C" {
    fn sqrt(_: libc::c_double) -> libc::c_double;
    fn fabs(_: libc::c_double) -> libc::c_double;
    fn memcpy(
        _: *mut libc::c_void,
        _: *const libc::c_void,
        _: libc::c_ulong,
    ) -> *mut libc::c_void;
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
pub unsafe extern "C" fn nlopt_isnan(mut x: libc::c_double) -> libc::c_int {
    return x.is_nan() as i32;
}
#[no_mangle]
pub unsafe extern "C" fn nlopt_isinf(mut x: libc::c_double) -> libc::c_int {
    return (fabs(x) >= ::std::f64::INFINITY * 0.99f64
        || if x.is_infinite() { if x.is_sign_positive() { 1 } else { -1 } } else { 0 }
            != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn nlopt_isfinite(mut x: libc::c_double) -> libc::c_int {
    return (fabs(x) <= 1.7976931348623157e+308f64 || x.is_finite() as i32 != 0)
        as libc::c_int;
}
unsafe extern "C" fn dcopy___(
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
unsafe extern "C" fn daxpy_sl__(
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
unsafe extern "C" fn ddot_sl__(
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
unsafe extern "C" fn dnrm2___(
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
        let mut xabs: libc::c_double = fabs(*dx.offset((incx * i) as isize));
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
    return xmax * sqrt(sum);
}
unsafe extern "C" fn dsrot_(
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
unsafe extern "C" fn dsrotg_(
    mut da: *mut libc::c_double,
    mut db: *mut libc::c_double,
    mut c: *mut libc::c_double,
    mut s: *mut libc::c_double,
) {
    let mut absa: libc::c_double = 0.;
    let mut absb: libc::c_double = 0.;
    let mut roe: libc::c_double = 0.;
    let mut scale: libc::c_double = 0.;
    absa = fabs(*da);
    absb = fabs(*db);
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
        r = (if roe < 0 as libc::c_int as libc::c_double { -scale } else { scale })
            * sqrt(tmpa * tmpa + tmpb * tmpb);
        *c = *da / r;
        *s = *db / r;
        *da = r;
        if *c != 0 as libc::c_int as libc::c_double && fabs(*c) <= *s {
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
unsafe extern "C" fn dscal_sl__(
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
unsafe extern "C" fn h12_(
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
        cl = fabs(d__1);
        if *mode == 2 as libc::c_int {
            if cl <= 0.0f64 {
                current_block = 5221028069996397600;
            } else {
                current_block = 15706997569754958587;
            }
        } else {
            i__1 = *m;
            j = *l1;
            while j <= i__1 {
                d__1 = *u.offset((j * u_dim1 + 1 as libc::c_int) as isize);
                sm = fabs(d__1);
                cl = if sm >= cl { sm } else { cl };
                j += 1;
            }
            if cl <= 0.0f64 {
                current_block = 5221028069996397600;
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
                cl *= sqrt(sm);
                if *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) > 0.0f64 {
                    cl = -cl;
                }
                *up = *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) - cl;
                *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize) = cl;
                current_block = 15706997569754958587;
            }
        }
        match current_block {
            5221028069996397600 => {}
            _ => {
                if !(*ncv <= 0 as libc::c_int) {
                    b = *up * *u.offset((*lpivot * u_dim1 + 1 as libc::c_int) as isize);
                    if !(b >= 0.0f64) {
                        b = one / b;
                        i2 = 1 as libc::c_int - *icv
                            + *ice * (*lpivot - 1 as libc::c_int);
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
                                sm
                                    += *c__.offset(i3 as isize)
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
                                    *c__.offset(i4 as isize)
                                        += sm
                                            * *u.offset((i__ * u_dim1 + 1 as libc::c_int) as isize);
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
unsafe extern "C" fn nnls_(
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
            &mut *x.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *x.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        'c_7799: loop {
            if iz1 > iz2 || nsetp >= *m {
                current_block = 12767425942910239388;
                break;
            }
            i__1 = iz2;
            iz = iz1;
            while iz <= i__1 {
                j = *indx.offset(iz as isize);
                i__2 = *m - nsetp;
                *w
                    .offset(
                        j as isize,
                    ) = ddot_sl__(
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
                    current_block = 12767425942910239388;
                    break 'c_7799;
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
                t = factor * fabs(d__1);
                d__1 = unorm + t;
                if !(d__1 - unorm <= 0.0f64) {
                    dcopy___(
                        m,
                        &mut *b.offset(1 as libc::c_int as isize),
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
                    if *z__.offset(npp1 as isize)
                        / *a.offset((npp1 + j * a_dim1) as isize) > 0.0f64
                    {
                        break;
                    }
                }
                *a.offset((npp1 + j * a_dim1) as isize) = asave;
                *w.offset(j as isize) = 0.0f64;
            }
            dcopy___(
                m,
                &mut *z__.offset(1 as libc::c_int as isize),
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
                &mut *w.offset(j as isize),
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
                            &mut d__1,
                            &mut *a.offset((jj * a_dim1 + 1 as libc::c_int) as isize),
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
                    current_block = 14151449179633537743;
                    break 'c_7799;
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
                    *x
                        .offset(
                            l as isize,
                        ) = (one - alpha) * *x.offset(l as isize)
                        + alpha * *z__.offset(ip as isize);
                    ip += 1;
                }
                if jj == 0 as libc::c_int {
                    break;
                }
                i__ = *indx.offset(jj as isize);
                'c_7847: loop {
                    *x.offset(i__ as isize) = 0.0f64;
                    jj += 1;
                    i__2 = nsetp;
                    j = jj;
                    while j <= i__2 {
                        ii = *indx.offset(j as isize);
                        *indx.offset((j - 1 as libc::c_int) as isize) = ii;
                        dsrotg_(
                            &mut *a
                                .offset((j - 1 as libc::c_int + ii * a_dim1) as isize),
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
                        current_block = 14151449179633537743;
                        break 'c_7799;
                    }
                    i__2 = nsetp;
                    jj = 1 as libc::c_int;
                    loop {
                        if !(jj <= i__2) {
                            break 'c_7847;
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
                    &mut *b.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *z__.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
            }
        }
        match current_block {
            14151449179633537743 => {
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
                &mut *w.offset(1 as libc::c_int as isize),
                0 as libc::c_int,
                &mut *w.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
            );
        }
    }
}
unsafe extern "C" fn ldp_(
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
            &mut *x.offset(1 as libc::c_int as isize),
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
                            *x
                                .offset(
                                    j as isize,
                                ) = fac
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
                            &mut *w.offset(1 as libc::c_int as isize),
                            0 as libc::c_int,
                            &mut *w.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        );
                        daxpy_sl__(
                            m,
                            &mut fac,
                            &mut *w.offset(iy as isize),
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
unsafe extern "C" fn lsi_(
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
            &mut i__3,
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
            if fabs(d__1) < epmach {
                current_block = 16236312898130099128;
                break 's_121;
            }
            i__3 = j - 1 as libc::c_int;
            *g
                .offset(
                    (i__ + j * g_dim1) as isize,
                ) = (*g.offset((i__ + j * g_dim1) as isize)
                - ddot_sl__(
                    &mut i__3,
                    &mut *g.offset((i__ + g_dim1) as isize),
                    *lg,
                    &mut *e.offset((j * e_dim1 + 1 as libc::c_int) as isize),
                    1 as libc::c_int,
                )) / *e.offset((j + j * e_dim1) as isize);
            j += 1;
        }
        *h__.offset(i__ as isize)
            -= ddot_sl__(
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
                    &mut *f.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                i__ = *n;
                while i__ >= 1 as libc::c_int {
                    i__2 = i__ + 1 as libc::c_int;
                    j = if i__2 <= *n { i__2 } else { *n };
                    i__2 = *n - i__;
                    *x
                        .offset(
                            i__ as isize,
                        ) = (*x.offset(i__ as isize)
                        - ddot_sl__(
                            &mut i__2,
                            &mut *e.offset((i__ + j * e_dim1) as isize),
                            *le,
                            &mut *x.offset(j as isize),
                            1 as libc::c_int,
                        )) / *e.offset((i__ + i__ * e_dim1) as isize);
                    i__ -= 1;
                }
                i__2 = *n + 1 as libc::c_int;
                j = if i__2 <= *me { i__2 } else { *me };
                i__2 = *me - *n;
                t = dnrm2___(&mut i__2, &mut *f.offset(j as isize), 1 as libc::c_int);
                *xnorm = sqrt(*xnorm * *xnorm + t * t);
            }
        }
        _ => {}
    };
}
unsafe extern "C" fn hfti_(
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
                current_block_56 = 14274665232496457393;
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
                    current_block_56 = 8942514791948061562;
                } else {
                    current_block_56 = 14274665232496457393;
                }
            }
            match current_block_56 {
                14274665232496457393 => {
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
                    *a
                        .offset(
                            (i__ + j * a_dim1) as isize,
                        ) = *a.offset((i__ + lmax * a_dim1) as isize);
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
                &mut i__3,
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
            if fabs(d__1) <= *tau {
                current_block = 8916208649952454714;
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
            *rnorm
                .offset(
                    jb as isize,
                ) = dnrm2___(
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
                        &mut i__2,
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
                    *b
                        .offset(
                            (i__ + jb * b_dim1) as isize,
                        ) = (*b.offset((i__ + jb * b_dim1) as isize)
                        - ddot_sl__(
                            &mut i__1,
                            &mut *a.offset((i__ + j * a_dim1) as isize),
                            *mda,
                            &mut *b.offset((j + jb * b_dim1) as isize),
                            1 as libc::c_int,
                        )) / *a.offset((i__ + i__ * a_dim1) as isize);
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
                        *b
                            .offset(
                                (l + jb * b_dim1) as isize,
                            ) = *b.offset((j + jb * b_dim1) as isize);
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
unsafe extern "C" fn lsei_(
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
        iw = (l + 1 as libc::c_int) * (*mg + 2 as libc::c_int)
            + (*mg << 1 as libc::c_int) + *mc;
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
                &mut i__3,
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
            if fabs(d__1) < epmach {
                current_block = 7757720978182437014;
                break;
            }
            i__1 = i__ - 1 as libc::c_int;
            *x
                .offset(
                    i__ as isize,
                ) = (*d__.offset(i__ as isize)
                - ddot_sl__(
                    &mut i__1,
                    &mut *c__.offset((i__ + c_dim1) as isize),
                    *lc,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                )) / *c__.offset((i__ + i__ * c_dim1) as isize);
            i__ += 1;
        }
        match current_block {
            7757720978182437014 => {}
            _ => {
                *mode = 1 as libc::c_int;
                *w.offset(mc1 as isize) = 0.0f64;
                i__2 = *mg;
                dcopy___(
                    &mut i__2,
                    &mut *w.offset(mc1 as isize),
                    0 as libc::c_int,
                    &mut *w.offset(mc1 as isize),
                    1 as libc::c_int,
                );
                if *mc == *n {
                    current_block = 13396851843167794670;
                } else {
                    i__2 = *me;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__2 {
                        *w
                            .offset(
                                (if__ - 1 as libc::c_int + i__) as isize,
                            ) = *f.offset(i__ as isize)
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
                            &mut *e.offset((i__ + mc1 * e_dim1) as isize),
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
                            &mut *g.offset((i__ + mc1 * g_dim1) as isize),
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
                            *h__.offset(i__ as isize)
                                -= ddot_sl__(
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
                            current_block = 7757720978182437014;
                        } else {
                            t = dnrm2___(
                                mc,
                                &mut *x.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            );
                            *xnrm = sqrt(*xnrm * *xnrm + t * t);
                            if *mode != 1 as libc::c_int {
                                current_block = 7757720978182437014;
                            } else {
                                current_block = 13396851843167794670;
                            }
                        }
                    } else {
                        *mode = 7 as libc::c_int;
                        k = if *le >= *n { *le } else { *n };
                        t = sqrt(epmach);
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
                            &mut *w.offset(if__ as isize),
                            1 as libc::c_int,
                            &mut *x.offset(mc1 as isize),
                            1 as libc::c_int,
                        );
                        if krank != l {
                            current_block = 7757720978182437014;
                        } else {
                            *mode = 1 as libc::c_int;
                            current_block = 13396851843167794670;
                        }
                    }
                }
                match current_block {
                    7757720978182437014 => {}
                    _ => {
                        i__2 = *me;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__2 {
                            *f
                                .offset(
                                    i__ as isize,
                                ) = ddot_sl__(
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
                            *d__
                                .offset(
                                    i__ as isize,
                                ) = ddot_sl__(
                                me,
                                &mut *e.offset((i__ * e_dim1 + 1 as libc::c_int) as isize),
                                1 as libc::c_int,
                                &mut *f.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            )
                                - ddot_sl__(
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
                            *w
                                .offset(
                                    i__ as isize,
                                ) = (*d__.offset(i__ as isize)
                                - ddot_sl__(
                                    &mut i__2,
                                    &mut *c__.offset((j + i__ * c_dim1) as isize),
                                    1 as libc::c_int,
                                    &mut *w.offset(j as isize),
                                    1 as libc::c_int,
                                )) / *c__.offset((i__ + i__ * c_dim1) as isize);
                            i__ -= 1;
                        }
                    }
                }
            }
        }
    }
}
unsafe extern "C" fn lsq_(
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
        diag = sqrt(*l.offset(i2 as isize));
        *w.offset(i3 as isize) = 0.0f64;
        dcopy___(
            &mut i1,
            &mut *w.offset(i3 as isize),
            0 as libc::c_int,
            &mut *w.offset(i3 as isize),
            1 as libc::c_int,
        );
        i__2 = i1 - n2;
        dcopy___(
            &mut i__2,
            &mut *l.offset(i2 as isize),
            1 as libc::c_int,
            &mut *w.offset(i3 as isize),
            *n,
        );
        i__2 = i1 - n2;
        dscal_sl__(&mut i__2, &mut diag, &mut *w.offset(i3 as isize), *n);
        *w.offset(i3 as isize) = diag;
        i__2 = i__ - 1 as libc::c_int;
        *w
            .offset(
                (if__ - 1 as libc::c_int + i__) as isize,
            ) = (*g.offset(i__ as isize)
            - ddot_sl__(
                &mut i__2,
                &mut *w.offset(i4 as isize),
                1 as libc::c_int,
                &mut *w.offset(if__ as isize),
                1 as libc::c_int,
            )) / diag;
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
            &mut *w.offset(i4 as isize),
            0 as libc::c_int,
            &mut *w.offset(i4 as isize),
            1 as libc::c_int,
        );
        *w.offset((if__ - 1 as libc::c_int + *n) as isize) = 0.0f64;
    }
    d__1 = -one;
    dscal_sl__(n, &mut d__1, &mut *w.offset(if__ as isize), 1 as libc::c_int);
    ic = if__ + *n;
    id = ic + *meq * *n;
    if *meq > 0 as libc::c_int {
        i__1 = *meq;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            dcopy___(
                n,
                &mut *a.offset((i__ + a_dim1) as isize),
                *la,
                &mut *w.offset((ic - 1 as libc::c_int + i__) as isize),
                *meq,
            );
            i__ += 1;
        }
        dcopy___(
            meq,
            &mut *b.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
            &mut *w.offset(id as isize),
            1 as libc::c_int,
        );
        d__1 = -one;
        dscal_sl__(meq, &mut d__1, &mut *w.offset(id as isize), 1 as libc::c_int);
    }
    ig = id + *meq;
    if mineq > 0 as libc::c_int {
        i__1 = mineq;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            dcopy___(
                n,
                &mut *a.offset((*meq + i__ + a_dim1) as isize),
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
            &mut *w.offset((ip - 1 as libc::c_int + i__) as isize),
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
            &mut *w.offset((im - 1 as libc::c_int + i__) as isize),
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
            &mut *b.offset((*meq + 1 as libc::c_int) as isize),
            1 as libc::c_int,
            &mut *w.offset(ih as isize),
            1 as libc::c_int,
        );
        d__1 = -one;
        dscal_sl__(&mut mineq, &mut d__1, &mut *w.offset(ih as isize), 1 as libc::c_int);
    }
    il = ih + mineq;
    iu = il + *n;
    i__ = 1 as libc::c_int;
    while i__ <= *n {
        *w
            .offset(
                (il - 1 as libc::c_int + i__) as isize,
            ) = if nlopt_isinf(*xl.offset(i__ as isize)) != 0 {
            0 as libc::c_int as libc::c_double
        } else {
            *xl.offset(i__ as isize)
        };
        *w
            .offset(
                (iu - 1 as libc::c_int + i__) as isize,
            ) = if nlopt_isinf(*xu.offset(i__ as isize)) != 0 {
            0 as libc::c_int as libc::c_double
        } else {
            -*xu.offset(i__ as isize)
        };
        i__ += 1;
    }
    iw = iu + *n;
    i__1 = if 1 as libc::c_int >= *meq { 1 as libc::c_int } else { *meq };
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
            &mut *w.offset(iw as isize),
            1 as libc::c_int,
            &mut *y.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        dcopy___(
            &mut n3,
            &mut *w.offset((iw + *m) as isize),
            1 as libc::c_int,
            &mut *y.offset((*m + 1 as libc::c_int) as isize),
            1 as libc::c_int,
        );
        dcopy___(
            &mut n3,
            &mut *w.offset((iw + *m + *n) as isize),
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
unsafe extern "C" fn ldl_(
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
unsafe extern "C" fn slsqpb_(
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
            *u
                .offset(
                    i__ as isize,
                ) = *g.offset(i__ as isize)
                - ddot_sl__(
                    m,
                    &mut *a.offset((i__ * a_dim1 + 1 as libc::c_int) as isize),
                    1 as libc::c_int,
                    &mut *r__.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                ) - *v.offset(i__ as isize);
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
                &mut h4,
                &mut *u.offset(1 as libc::c_int as isize),
                1 as libc::c_int,
            );
            d__1 = one - h4;
            daxpy_sl__(
                n,
                &mut d__1,
                &mut *v.offset(1 as libc::c_int as isize),
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
        current_block = 17386100730026577537;
    } else if *mode == 0 as libc::c_int {
        itermx = *iter;
        if *acc >= 0.0f64 {
            iexact = 0 as libc::c_int;
        } else {
            iexact = 1 as libc::c_int;
        }
        *acc = fabs(*acc);
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
            &mut *s.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *s.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        dcopy___(
            m,
            &mut *mu.offset(1 as libc::c_int as isize),
            0 as libc::c_int,
            &mut *mu.offset(1 as libc::c_int as isize),
            1 as libc::c_int,
        );
        current_block = 3938304504335184834;
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
                current_block = 13298064082901038618;
                match current_block {
                    13298064082901038618 => {
                        if nlopt_isfinite(h1) != 0 {
                            if h1 <= h3 / ten || line > 10 as libc::c_int {
                                current_block = 778573858929849943;
                            } else {
                                d__1 = h3 / (two * (h3 - h1));
                                alpha = if d__1 >= alfmin { d__1 } else { alfmin };
                                current_block = 2990001543308207002;
                            }
                        } else {
                            alpha = if alpha * 0.5f64 >= alfmin {
                                alpha * 0.5f64
                            } else {
                                alfmin
                            };
                            current_block = 2990001543308207002;
                        }
                    }
                    _ => {}
                }
                match current_block {
                    2990001543308207002 => {}
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
                        if (fabs(d__1) < *acc
                            || dnrm2___(
                                n,
                                &mut *s.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            ) < *acc) && h3 < *acc
                        {
                            *mode = 0 as libc::c_int;
                        } else {
                            *mode = -(1 as libc::c_int);
                        }
                        current_block = 5266944821211112966;
                    }
                }
            }
            2 => {
                current_block = 17680132173918469312;
            }
            _ => {
                current_block = 778573858929849943;
                match current_block {
                    13298064082901038618 => {
                        if nlopt_isfinite(h1) != 0 {
                            if h1 <= h3 / ten || line > 10 as libc::c_int {
                                current_block = 778573858929849943;
                            } else {
                                d__1 = h3 / (two * (h3 - h1));
                                alpha = if d__1 >= alfmin { d__1 } else { alfmin };
                                current_block = 2990001543308207002;
                            }
                        } else {
                            alpha = if alpha * 0.5f64 >= alfmin {
                                alpha * 0.5f64
                            } else {
                                alfmin
                            };
                            current_block = 2990001543308207002;
                        }
                    }
                    _ => {}
                }
                match current_block {
                    2990001543308207002 => {}
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
                        if (fabs(d__1) < *acc
                            || dnrm2___(
                                n,
                                &mut *s.offset(1 as libc::c_int as isize),
                                1 as libc::c_int,
                            ) < *acc) && h3 < *acc
                        {
                            *mode = 0 as libc::c_int;
                        } else {
                            *mode = -(1 as libc::c_int);
                        }
                        current_block = 5266944821211112966;
                    }
                }
            }
        }
    }
    'c_2994: loop {
        match current_block {
            17680132173918469312 => {
                *mode = 9 as libc::c_int;
                return;
            }
            2990001543308207002 => {
                line += 1;
                h3 = alpha * h3;
                dscal_sl__(
                    n,
                    &mut alpha,
                    &mut *s.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                dcopy___(
                    n,
                    &mut *x0.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                    &mut *x.offset(1 as libc::c_int as isize),
                    1 as libc::c_int,
                );
                daxpy_sl__(
                    n,
                    &one,
                    &mut *s.offset(1 as libc::c_int as isize),
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
                current_block = 5266944821211112966;
            }
            17386100730026577537 => {
                *iter += 1;
                *mode = 9 as libc::c_int;
                if *iter > itermx && itermx > 0 as libc::c_int {
                    current_block = 5266944821211112966;
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
                        &mut d__1,
                        &mut *x.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                        &mut *u.offset(1 as libc::c_int as isize),
                        1 as libc::c_int,
                    );
                    d__1 = -one;
                    daxpy_sl__(
                        n,
                        &mut d__1,
                        &mut *x.offset(1 as libc::c_int as isize),
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
                        &mut *u.offset(1 as libc::c_int as isize),
                        &mut *v.offset(1 as libc::c_int as isize),
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
                                *a
                                    .offset(
                                        (j + n1 * a_dim1) as isize,
                                    ) = -*c__.offset(j as isize);
                            } else {
                                d__1 = -*c__.offset(j as isize);
                                *a
                                    .offset(
                                        (j + n1 * a_dim1) as isize,
                                    ) = if d__1 >= 0.0f64 { d__1 } else { 0.0f64 };
                            }
                            j += 1;
                        }
                        *s.offset(1 as libc::c_int as isize) = 0.0f64;
                        dcopy___(
                            n,
                            &mut *s.offset(1 as libc::c_int as isize),
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
                                &mut *u.offset(1 as libc::c_int as isize),
                                &mut *v.offset(1 as libc::c_int as isize),
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
                                current_block = 5266944821211112966;
                                continue 'c_2994;
                            }
                        }
                        if *mode != 1 as libc::c_int {
                            current_block = 5266944821211112966;
                            continue;
                        }
                    } else if *mode != 1 as libc::c_int {
                        current_block = 5266944821211112966;
                        continue;
                    }
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        *v
                            .offset(
                                i__ as isize,
                            ) = *g.offset(i__ as isize)
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
                        &mut *x.offset(1 as libc::c_int as isize),
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
                    h1 = fabs(gs);
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
                        h3 = fabs(d__1);
                        d__1 = h3;
                        d__2 = (*mu.offset(j as isize) + h3) / two;
                        *mu.offset(j as isize) = if d__1 >= d__2 { d__1 } else { d__2 };
                        d__1 = *c__.offset(j as isize);
                        h1 += h3 * fabs(d__1);
                        j += 1;
                    }
                    *mode = 0 as libc::c_int;
                    if h1 < *acc && h2 < *acc {
                        current_block = 5266944821211112966;
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
                        h1
                            += *mu.offset(j as isize)
                                * (if d__1 >= h3 { d__1 } else { h3 });
                        j += 1;
                    }
                    t0 = *f + h1;
                    h3 = gs - h1 * h4;
                    *mode = 8 as libc::c_int;
                    if h3 >= 0.0f64 {
                        current_block = 3938304504335184834;
                        continue;
                    }
                    line = 0 as libc::c_int;
                    alpha = one;
                    if iexact == 1 as libc::c_int {
                        current_block = 17680132173918469312;
                    } else {
                        current_block = 2990001543308207002;
                    }
                }
            }
            3938304504335184834 => {
                ireset += 1;
                if ireset > 5 as libc::c_int {
                    d__1 = *f - f0;
                    if (fabs(d__1) < tol
                        || dnrm2___(
                            n,
                            &mut *s.offset(1 as libc::c_int as isize),
                            1 as libc::c_int,
                        ) < tol) && h3 < tol
                    {
                        *mode = 0 as libc::c_int;
                    } else {
                        *mode = 8 as libc::c_int;
                    }
                    current_block = 5266944821211112966;
                } else {
                    *l.offset(1 as libc::c_int as isize) = 0.0f64;
                    dcopy___(
                        &mut n2,
                        &mut *l.offset(1 as libc::c_int as isize),
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
                    current_block = 17386100730026577537;
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
    };
}
#[no_mangle]
pub unsafe extern "C" fn slsqp(
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
        + (mineq << 1 as libc::c_int) + (n1 + mineq) * (n1 - *meq)
        + (*meq << 1 as libc::c_int) + n1 * *n / 2 as libc::c_int
        + (*m << 1 as libc::c_int) + *n * 3 as libc::c_int + (n1 << 2 as libc::c_int)
        + 1 as libc::c_int;
    i__1 = mineq;
    i__2 = n1 - *meq;
    im = if i__1 >= i__2 { i__1 } else { i__2 };
    if *l_w__ < il || *l_jw__ < im {
        *mode = (if 10 as libc::c_int >= il { 10 as libc::c_int } else { il })
            * 1000 as libc::c_int;
        *mode += if 10 as libc::c_int >= im { 10 as libc::c_int } else { im };
        return;
    }
    im = 1 as libc::c_int;
    il = im + (if 1 as libc::c_int >= *m { 1 as libc::c_int } else { *m });
    il = im + *la;
    ix = il + n1 * *n / 2 as libc::c_int + 1 as libc::c_int;
    ir = ix + *n;
    is = ir + *n + *n + (if 1 as libc::c_int >= *m { 1 as libc::c_int } else { *m });
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
    let ref mut fresh0 = (*state).x0;
    *fresh0 = &mut *w.offset(ix as isize) as *mut libc::c_double;
}
