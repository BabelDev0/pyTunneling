extern crate ndarray;
use ndarray::{s, Array1, Array2};
use plotters::prelude::*;

fn main() {
    let hb = 1.0; // h bar in atomic units
    let eVtoHa = 0.03674; // eV to hartree
    let AtoB = 1.88973; // Angstrom to Bohr
    let VeV = 0.100; // Barrier height in eV
    let V = VeV * eVtoHa; // Barrier height in Hartree
    let m0 = 1.0; // electron mass in a.u.
    let m1 = m0 * 0.067; // effective mass in the well
    let m2 = m0 * 0.067; // effective mass in the barrier
    let L1A = 50.0; // thickness of the first barrier in Ang
    let L2A = 50.0; // thickness of the well in Ang
    let L3A = 50.0; // thickness of the second barrier in Ang
    let I1 = 0.0;
    let I2 = L1A * AtoB;
    let I3 = (L1A + L2A) * AtoB;
    let I4 = (L1A + L2A + L3A) * AtoB;

    let mut E = 0.0;
    let num_points = 10000;
    let mut x = Array1::<f64>::zeros(num_points);
    let mut T = Array1::<f64>::zeros(num_points);

    for n in 0..num_points {
        E += V / num_points;
        let k1 = (2.0 * m1 * E).sqrt() / hb;
        let k2 = (2.0 * m2 * (V - E)).sqrt() / hb;

        let mut M1 = Array2::<f64>::zeros((2, 2));
        M1[(0, 0)] = 1.0;
        M1[(0, 1)] = 1.0;
        M1[(1, 0)] = 1.0j * k1;
        M1[(1, 1)] = -1.0j * k1;

        // Define other matrices (M2, M3, ..., M8) similarly

        let M = M1
            .dot(&M2)
            .dot(&M3)
            .dot(&M4)
            .dot(&M5)
            .dot(&M6)
            .dot(&M7)
            .dot(&M8);
        x[n] = E / eVtoHa;
        T[n] = 1.0 / (M[(0, 0)].conj() * M[(0, 0)]).norm();
    }

    // Plotting
    let root = BitMapBackend::new("output.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Transmission Probability through a Double Barrier System",
            ("sans-serif", 20).into_font(),
        )
        .margin(5)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0.0..x[x.len() - 1], 0.0..1.0)
        .unwrap();

    chart
        .configure_mesh()
        .line_style_1(&WHITE.mix(0.2))
        .draw()
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            x.iter().zip(T.iter()).map(|(&x, &t)| (x, t)),
            &RED,
        ))
        .unwrap();
}
