use fast_float::parse;
use plotters::prelude::*;
use rayon::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
use std::error::Error;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::thread;

use crate::AVG_LEN;
use crate::STRINGS;
use crate::THREADS;
use crate::T_RES;

use crate::NotePeak;
use crate::NotesToIndex;
use crate::Peaks;

pub fn single_plot_data_norm(data: Vec<f32>) {
    draw_plot("plots/plot0.png", data.to_vec(), T_RES as f32, 1).unwrap();
}

pub fn plot_data_norm(data: &Vec<Vec<Vec<f32>>>, prefix: &str, sec: f32) {
    let mut filename = format!("plots/{prefix}");
    let mut handles = vec![];
    for i in 0..6 {
        let mut string_data = data[i].clone();
        let filename_base = filename.clone();

        handles.push(thread::spawn(move || {
            for n in 0..20 {
                let note_data = string_data.remove(0);
                let len = note_data.len() as f32;
                let filename = format!("{filename_base}{}{}.png", STRINGS[i], n);
                draw_plot(&filename, note_data, sec / len, 1).unwrap();
            }
        }));
    }

    let mut joined_data: Vec<f32> = Vec::new();
    for handle in handles {
        handle.join().unwrap();
    }
}

pub fn plot_data(data: Vec<Vec<Complex<f32>>>) {
    for note in 0..data.len() {
        let mut filename = String::from("plots/");
        filename.push_str(&note.to_string());
        filename.push_str("_conv.png");
        draw_plot(
            &filename,
            data[note].iter().map(|a| a.norm()).collect(),
            T_RES as f32,
            1,
        )
        .unwrap();
    }
}

pub fn draw_plot(
    plot_name: &str,
    data: Vec<f32>,
    time: f32,
    div: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(plot_name, (1024, 768)).into_drawing_area();
    let mut freq: Vec<f32> = Vec::new();
    let max = 0.005;
    let temp: Vec<f32> = vec![0.0, max];
    let max = temp.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
    let max = data.iter().max_by(|a, b| a.total_cmp(b)).unwrap();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(plot_name, ("sans-serif", 40))
        .build_cartesian_2d(0..(data.len() / div - 1), 0.0..*max)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_label_formatter(&|x| format!("{}", (*x as f32 * time) as f32))
        .draw()?;

    chart.draw_series(
        AreaSeries::new(
            (0..).zip(data.iter()).map(|(x, y)| (x, *y)),
            0.0,
            &RED.mix(0.2),
        )
        .border_style(&RED),
    )?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", plot_name);
    println!("max is {max}");
    Ok(())
}
