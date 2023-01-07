use std::io::Write;
use std::io::IoSlice;
use std::thread;
use std::path::Path;
use std::error::Error;
use rustfft::{FftPlanner, num_complex::Complex};
use fast_float::parse;
use rayon::prelude::*;
use plotters::prelude::*;

use crate::SAMPLE;
use crate::AVG_LEN;
use crate::THREADS;
use crate::T_RES;
use crate::F_RES;
use crate::STRINGS;

use crate::NotePeak;
use crate::NotesToIndex;
use crate::Peaks;

pub fn single_plot_data_norm(data: Vec<f32>){
    draw_plot("plots/plot.png",data.to_vec(), T_RES as f32, 1).unwrap();
}

pub fn plot_data_norm(data: Vec<Vec<Vec<f32>>>){
    for string in 0..data.len(){
        for note in 0..data[0].len(){
            let mut filename = String::from("plots/");
            filename.push_str(STRINGS[string]);
            filename.push_str(&note.to_string());
            filename.push_str(".png");
            draw_plot(&filename,data[string][note].to_vec(), T_RES * AVG_LEN as f32, 1).unwrap();
        }
    }
}

pub fn plot_data(data: Vec<Vec<Complex<f32>>>){
    for note in 0..data.len(){
        let mut filename = String::from("plots/");
        filename.push_str(&note.to_string());
        filename.push_str("_conv.png");
        draw_plot(&filename,data[note].iter()
                  .map(|a| a.norm()).collect(),
                  T_RES * AVG_LEN as f32, 1).unwrap();
    }
}

pub fn draw_plot(plot_name: &str, data: Vec<f32>,
             mul: f32, div: usize) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(plot_name, (1024, 768)).into_drawing_area();
    let mut freq: Vec<f32> = Vec::new();
    let max = 1.0;
    let temp: Vec<f32> = vec![0.0,max];
    let max = temp.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
    let max = data.iter().max_by(|a, b| a.total_cmp(b)).unwrap();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption("BRUH", ("sans-serif", 40))
        .build_cartesian_2d(0..(data.len()/div - 1), 0.0..*max)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_label_formatter(&|x| format!("{}", (*x as f32 * mul) as f32))
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
