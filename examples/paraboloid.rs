use slsqp::{minimize, Func};

/// Problem cost function
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

fn main() {
    let x = vec![1., 1.];

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
        &x,
        &cons,
        (),
        200,
        &[(-10., 10.)],
        1e-4,
        0.0,
        0.0,
        &[0.0, 0.0],
    ) {
        Ok((status, x_opt, y_opt)) => {
            println!("status = {:?}", status);
            println!("x_opt = {:?}", x_opt);
            println!("y_opt = {}", y_opt);
        }
        Err((e, _, _)) => println!("Optim error: {:?}", e),
    }
}
