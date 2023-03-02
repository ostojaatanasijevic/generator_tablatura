#![allow(warnings, unused)]
// KORISTI REALFFT
// ADAPTIRAJ INTENZITET SIGURNOSTI U NOTU U POREDJENU SA INTENZITETIOM PESME
// TEMPO SE RAČUNA AUTOKORELACIOJOM SIGNALA

mod cli;
mod fft;
mod fourier;
mod freq_generator;
mod harmonics;
mod misc;
mod plot;
mod post_processing;

use std::cell::RefCell;
use std::rc::Rc;

use cairo::glib::Receiver;
use gtk::prelude::*;
use plotters::prelude::*;
use plotters_cairo::CairoBackend;

use clap::Parser;
use fast_float::parse;
use harmonics::Note;
use misc::hertz_to_notes;
use num::complex::ComplexFloat;
use num::traits::Pow;
use num::FromPrimitive;
use rayon::vec;
use std::cmp;
use std::env;
use std::error::Error;
use std::f32::consts::PI;
use std::fs::File;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc;
use std::thread;
use std::time::Instant;

use realfft::RealFftPlanner;
use rustfft::{num_complex::Complex, FftPlanner};

use crate::harmonics::attenuate_neighbours;

pub const THREADS: usize = 8;
pub const AVG_LEN: usize = 1; // sample_len / 32; // mora da može da deli NFFT, da ne bi cureli podaci
pub const T_RES: f32 = 1.0 / (44100.0) * AVG_LEN as f32;
pub const STRINGS: [&str; 6] = ["e", "B", "G", "D", "A", "E"];

pub const HERZ: [&str; 44] = [
    "82.41", "87.31", "92.5", "98.0", "103.83", "110.0", "116.54", "123.47", "130.81", "138.59",
    "146.83", "155.56", "164.81", "174.61", "185.0", "196.0", "207.65", "220.0", "233.08",
    "246.94", "261.63", "277.18", "293.66", "311.13", "329.63", "349.23", "369.99", "392.0",
    "415.3", "440.0", "466.16", "493.88", "523.25", "554.37", "587.33", "622.25", "659.26",
    "698.46", "739.99", "783.99", "830.61", "880.0", "932.33", "987.77",
];

const OFFSET_TABLE: [usize; 6] = [0, 5, 10, 15, 19, 24];
const DIFF_TABLE: [usize; 6] = [20, 5, 4, 5, 5, 5];

//ADD THREAD DETECTION FOR INDIVIDUAL CPUs
fn main() {
    let args = cli::Args::parse();
    let (tx, rx): (mpsc::Sender<i32>, mpsc::Receiver<i32>) = mpsc::channel();

    let th = thread::spawn(move || {
        println!("gas");
        tx.send(1).unwrap();
    });

    rx.recv();
    th.join().unwrap();

    let note_intensity = process_song(&args);

    let application = gtk::Application::new(
        Some("io.github.plotters-rs.plotters-gtk-demo"),
        Default::default(),
    );

    application.connect_activate(move |app| {
        build_ui(app, &note_intensity, args.sec_to_run, &args);
    });

    application.run_with_args(&[""]);
}
const GLADE_UI_SOURCE: &'static str = include_str!("ui.glade");

#[derive(Clone)]
struct PlottingState {
    x_offset: f64,
    y_scale: f64,
    std_x: f64,
    std_y: f64,
    pitch: f64,
    time_frame: f64,
    time_res: f32,
    data_orig: Vec<Vec<Vec<f32>>>,
    data_att: Vec<Vec<Vec<f32>>>,
    data_out: Vec<Vec<Vec<f32>>>,
    low_threshold: f64,
    high_threshold: f64,
    args: cli::Args,
    time_processed: f64,
}

impl PlottingState {
    fn plot<'a, DB: DrawingBackend + 'a>(
        &self,
        backend: DB,
        data: &Vec<Vec<Vec<f32>>>,
        data_out: &Vec<Vec<Vec<f32>>>,
    ) -> Result<(), Box<dyn Error + 'a>> {
        let samples_per_sec = data[0][0].len() as f64 / self.time_processed;
        let x_offset = (self.x_offset * samples_per_sec as f64) as usize;

        let end_point = cmp::min(
            x_offset + (self.time_frame / self.time_processed * data[0][0].len() as f64) as usize,
            data[0][0].len(),
        );

        let mut max: f32 = 0.1;

        let root = backend.into_drawing_area();

        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .set_label_area_size(LabelAreaPosition::Left, 140)
            .set_label_area_size(LabelAreaPosition::Bottom, 60)
            .build_cartesian_2d(0..end_point - x_offset, 0.0..max)?;

        println!("Reploting");

        chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .x_label_formatter(&|x| format!("{}", ((*x + x_offset) as f64 / samples_per_sec)))
            .y_labels(0)
            .draw()?;

        for note in 0..data[0].len() {
            chart.draw_series(
                AreaSeries::new(
                    (0..)
                        .zip(data[self.pitch as usize][note][x_offset..end_point].iter())
                        .map(|(x, y)| (x, *y / 20.0 * self.y_scale as f32 + note as f32 / 200.0)),
                    note as f32 / 200.0,
                    &BLUE.mix(0.2),
                )
                .border_style(&BLUE),
            )?;

            chart.draw_series(
                AreaSeries::new(
                    (0..)
                        .zip(data_out[self.pitch as usize][note][x_offset..end_point].iter())
                        .map(|(x, y)| (x, *y / 20.0 * self.y_scale as f32 + note as f32 / 200.0)),
                    note as f32 / 200.0,
                    &RED.mix(0.2),
                )
                .border_style(&RED),
            )?;

            let mut x = 10;
            let y = ((20 - note) as i32 * root.get_pixel_range().1.end as i32 / 21) as i32;

            if y > 0 {
                root.draw(&Text::new(
                    misc::index_to_note(note + OFFSET_TABLE[5 - self.pitch as usize]),
                    (x, y),
                    ("sans-serif", 15.0).into_font(),
                ));
            }
        }

        root.present()?;
        Ok(())
    }
}

fn build_ui(app: &gtk::Application, data: &Vec<Vec<Vec<f32>>>, sec: f32, args: &cli::Args) {
    let window = fourier::calculate_window_function(args.nfft, &args.window_function); // blackman je bolji od hann
    let sample_notes = freq_generator::open_sample_notes(args.nfft);
    let mut all_notes = harmonics::generate_all_notes();
    harmonics::generate_note_network(&mut all_notes, &sample_notes, &window, args.nfft);
    harmonics::cross_polinate(&mut all_notes, &sample_notes, &window, args.nfft);

    print!("{:?}", &all_notes[0]);

    let builder = gtk::Builder::from_string(GLADE_UI_SOURCE);
    let window = builder.object::<gtk::Window>("MainWindow").unwrap();

    window.set_title("Generator tablatura");

    let drawing_area: gtk::DrawingArea = builder.object("MainDrawingArea").unwrap();
    let string_slider = builder.object::<gtk::Scale>("StringSlider").unwrap();
    let time_frame_slider = builder.object::<gtk::Scale>("TimeFrameSlider").unwrap();
    let time_slider = builder.object::<gtk::Scale>("TimeSlider").unwrap();
    let y_scale_slider = builder.object::<gtk::Scale>("MeanYScale").unwrap();
    let std_x_scale = builder.object::<gtk::Scale>("SDXScale").unwrap();
    let std_y_scale = builder.object::<gtk::Scale>("SDYScale").unwrap();

    let low_threshold_slider = builder.object::<gtk::Scale>("LowThresholdSlider").unwrap();
    let high_threshold_slider = builder.object::<gtk::Scale>("HighThresholdSlider").unwrap();

    let app_state = Rc::new(RefCell::new(PlottingState {
        x_offset: time_slider.value(),
        y_scale: y_scale_slider.value(),
        std_x: std_x_scale.value(),
        std_y: std_y_scale.value(),
        pitch: string_slider.value(),
        time_frame: time_frame_slider.value(),
        time_res: sec / data.len() as f32,
        data_orig: data.clone(),
        data_att: data.clone(),
        data_out: data.clone(),
        low_threshold: 0.0,
        high_threshold: 1.0,
        args: args.clone(),
        time_processed: args.sec_to_run as f64,
    }));

    window.set_application(Some(app));

    let state_cloned = app_state.clone();

    drawing_area.connect_draw(move |widget, cr| {
        let state = state_cloned.borrow().clone();
        let w = widget.allocated_width();
        let h = widget.allocated_height();
        let backend = CairoBackend::new(cr, (w as u32, h as u32)).unwrap();
        state
            .plot(backend, &state.data_att, &state.data_out)
            .unwrap();
        Inhibit(false)
    });

    let handle_change =
        |what: &gtk::Scale, how: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();
            what.connect_value_changed(move |target| {
                let mut state = app_state.borrow_mut();
                *how(&mut *state) = target.value();
                drawing_area.queue_draw();
            });
        };
    /*
    let process_more_data =
        |slider: &gtk::Scale,
         time: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>,
         data_orig: Box<dyn Fn(&mut PlottingState) -> &mut Vec<Vec<Vec<f32>>> + 'static>,
         time_processed: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();
            slider.connect_value_changed(move |target| {
                let mut state = app_state.borrow_mut();
                *time(&mut *state) = target.value();

                if (target.value() > state.time_processed - state.time_frame) {
                    let mut args = state.args.clone();
                    args.seek_offset = state.time_processed as f32;
                    args.sec_to_run = 5.0;
                    *time_processed(&mut *state) += args.sec_to_run as f64;
                    let temp = process_song(&args);
                    *data_orig(&mut *state) = weld_note_intensity(&state.data_orig, &temp);
                }

                drawing_area.queue_draw();
            });
        };
    process_more_data(
        &time_slider,
        Box::new(|s| &mut s.x_offset),
        Box::new(|s| &mut s.data_orig),
        Box::new(|s| &mut s.time_processed),
    );
    */

    let process_more_data_exp = |slider: &gtk::Scale,
                                 ploting_state: Box<
        dyn Fn(&mut PlottingState) -> &mut PlottingState + 'static,
    >| {
        let app_state = app_state.clone();
        let drawing_area = drawing_area.clone();
        slider.connect_value_changed(move |target| {
            let mut state = app_state.borrow_mut();
            ploting_state(&mut *state).x_offset = target.value();

            if (target.value() > state.time_processed - state.time_frame) {
                let mut args = state.args.clone();
                args.seek_offset = state.time_processed as f32;
                args.sec_to_run = 5.0;
                ploting_state(&mut *state).time_processed += args.sec_to_run as f64;
                let temp = process_song(&args);
                ploting_state(&mut *state).data_orig = weld_note_intensity(&state.data_orig, &temp);
            }

            drawing_area.queue_draw();
        });
    };

    let recalculate_attenuation =
        |slider: &gtk::Scale,
         ploting_state: Box<dyn Fn(&mut PlottingState) -> &mut PlottingState + 'static>,
         factor: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();
            let all_notes = all_notes.clone();
            slider.connect_value_changed(move |target| {
                let mut state = app_state.borrow_mut();
                *factor(&mut *state) = target.value();
                let temp = harmonics::attenuate_harmonics(
                    &state.data_orig,
                    &all_notes,
                    state.std_x as f32,
                    state.std_y as f32,
                );

                ploting_state(&mut *state).data_out = post_processing::schmitt(
                    &temp,
                    state.high_threshold as f32,
                    state.low_threshold as f32,
                );
                ploting_state(&mut *state).data_att = temp;

                drawing_area.queue_draw();
            });
        };

    let recalculate_thresholds =
        |what: &gtk::Scale,
         how: Box<dyn Fn(&mut PlottingState) -> &mut Vec<Vec<Vec<f32>>> + 'static>,
         factor: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();

            let all_notes = all_notes.clone();
            what.connect_value_changed(move |target| {
                let mut state = app_state.borrow_mut();
                *factor(&mut *state) = target.value();
                *how(&mut *state) = post_processing::schmitt(
                    &state.data_att,
                    state.high_threshold as f32,
                    state.low_threshold as f32,
                );

                drawing_area.queue_draw();
            });
        };

    handle_change(&string_slider, Box::new(|s| &mut s.pitch));
    handle_change(&time_frame_slider, Box::new(|s| &mut s.time_frame));

    process_more_data_exp(&time_slider, Box::new(|s| s));
    handle_change(&y_scale_slider, Box::new(|s| &mut s.y_scale));

    recalculate_attenuation(&std_x_scale, Box::new(|s| s), Box::new(|s| &mut s.std_x));
    recalculate_attenuation(&std_y_scale, Box::new(|s| s), Box::new(|s| &mut s.std_y));

    recalculate_thresholds(
        &low_threshold_slider,
        Box::new(|s| &mut s.data_out),
        Box::new(|s| &mut s.low_threshold),
    );

    recalculate_thresholds(
        &high_threshold_slider,
        Box::new(|s| &mut s.data_out),
        Box::new(|s| &mut s.high_threshold),
    );

    window.show_all();
}

fn process_song(args: &cli::Args) -> Vec<Vec<Vec<f32>>> {
    println!("Processing...");

    let h = post_processing::lp_filter(args.w, args.lenght_fir);
    let window = fourier::calculate_window_function(args.nfft, &args.window_function); // blackman je bolji od hann

    let sample_notes = freq_generator::open_sample_notes(args.nfft);
    let mut all_notes = harmonics::generate_all_notes();
    harmonics::generate_note_network(&mut all_notes, &sample_notes, &window, args.nfft);
    harmonics::cross_polinate(&mut all_notes, &sample_notes, &window, args.nfft);

    let song = fourier::open_song(&args.file_name, args.seek_offset, args.sec_to_run);
    let mut note_intensity =
        fft::threaded_interlaced_convolution_realfft(&song, &sample_notes, &window, &args);
    let mut note_intensity = post_processing::threaded_fft_fir_filtering(note_intensity, &h, &args);
    attenuate_neighbours(&note_intensity, &all_notes, 0.75, 1.0)
}

fn weld_note_intensity(
    first_part: &Vec<Vec<Vec<f32>>>,
    second_part: &Vec<Vec<Vec<f32>>>,
) -> Vec<Vec<Vec<f32>>> {
    let mut out = first_part.clone();

    for wire in 0..6 {
        for note in 0..20 {
            out[wire][note].extend(second_part[wire][note].to_vec());
        }
    }

    out
}
