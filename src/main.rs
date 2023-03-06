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
use gtk::builders::ScrollbarBuilder;
use gtk::ffi::gtk_widget_set_events;
use gtk::ffi::GTK_SCROLL_START;
use gtk::prelude::*;
use gtk::Scrollbar;
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
use std::cmp::max;
use std::cmp::min;
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

const GLADE_UI_SOURCE: &'static str = include_str!("ui.glade");
pub const THREADS: usize = 8;
pub const AVG_LEN: usize = 1; // sample_len / 32; // mora da može da deli NFFT, da ne bi cureli podaci
pub const T_RES: f32 = 1.0 / (44100.0) * AVG_LEN as f32;
pub const STRINGS: [&str; 6] = ["e", "B", "G", "D", "A", "E"];
pub const OFFSET_TABLE: [usize; 6] = [0, 5, 10, 15, 19, 24];
pub const DIFF_TABLE: [usize; 6] = [20, 5, 4, 5, 5, 5];

pub const BROJ_ZICA: usize = 6;
pub const BROJ_PRAGOVA: usize = 20;

pub const HERZ: [&str; 44] = [
    "82.41", "87.31", "92.5", "98.0", "103.83", "110.0", "116.54", "123.47", "130.81", "138.59",
    "146.83", "155.56", "164.81", "174.61", "185.0", "196.0", "207.65", "220.0", "233.08",
    "246.94", "261.63", "277.18", "293.66", "311.13", "329.63", "349.23", "369.99", "392.0",
    "415.3", "440.0", "466.16", "493.88", "523.25", "554.37", "587.33", "622.25", "659.26",
    "698.46", "739.99", "783.99", "830.61", "880.0", "932.33", "987.77",
];

//ADD THREAD DETECTION FOR INDIVIDUAL CPUs
fn main() {
    let args = cli::Args::parse();

    let application = gtk::Application::new(
        Some("ostojin.generator.tablatura.realdank"),
        Default::default(),
    );

    application.connect_activate(move |app| {
        build_ui(app, &args);
    });

    application.run_with_args(&[""]);
}

#[derive(Clone)]
struct PlottingState {
    time_offset: f64,
    y_slider: f64,
    attenuation_factor: f64,
    attenuation_power_factor: f64,
    current_string: f64,
    song: Vec<i16>,
    time_frame: f64,
    all_notes: Vec<Note>,
    data_orig: Vec<Vec<Vec<f32>>>,
    data_att: Vec<Vec<Vec<f32>>>,
    data_out: Vec<Vec<Vec<f32>>>,
    h: Vec<f32>,
    window_fn: Vec<f32>,
    sample_notes: Vec<Vec<Vec<i16>>>,
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
        let time_offset = (self.time_offset * samples_per_sec as f64) as usize;

        let end_point = min(
            time_offset + (self.time_frame * samples_per_sec) as usize,
            data[0][0].len(),
        );

        let mut max: f32 = 0.1;

        let root = backend.into_drawing_area();

        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .set_label_area_size(LabelAreaPosition::Left, 140)
            .set_label_area_size(LabelAreaPosition::Bottom, 60)
            .build_cartesian_2d(0..end_point - time_offset, 0.0..max)?;

        println!("Reploting");

        chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .x_label_formatter(&|x| format!("{:.2}", ((*x + time_offset) as f64 / samples_per_sec)))
            .y_labels(0)
            .draw()?;

        for note in 0..data[0].len() {
            chart.draw_series(
                AreaSeries::new(
                    (0..)
                        .zip(
                            data[self.current_string as usize][note][time_offset..end_point].iter(),
                        )
                        .map(|(x, y)| (x, *y / 20.0 * self.y_slider as f32 + note as f32 / 200.0)),
                    note as f32 / 200.0,
                    &BLUE.mix(0.2),
                )
                .border_style(&BLUE),
            )?;

            chart.draw_series(
                AreaSeries::new(
                    (0..)
                        .zip(
                            data_out[self.current_string as usize][note][time_offset..end_point]
                                .iter(),
                        )
                        .map(|(x, y)| (x, *y / 20.0 * self.y_slider as f32 + note as f32 / 200.0)),
                    note as f32 / 200.0,
                    &RED.mix(0.2),
                )
                .border_style(&RED),
            )?;

            let y = ((20 - note) as i32 * root.get_pixel_range().1.end as i32 / 21) as i32;
            root.draw(&Text::new(
                misc::index_to_note(note + OFFSET_TABLE[5 - self.current_string as usize]),
                (10, y),
                ("sans-serif", 15.0).into_font(),
            ));
        }
        root.present()?;
        Ok(())
    }
}

fn build_ui(app: &gtk::Application, args: &cli::Args) {
    let args = args.clone();
    let song = fourier::open_song(&args.file_name);
    let h = post_processing::lp_filter(args.w, args.lenght_fir);
    let window_fn = fourier::calculate_window_function(args.nfft, &args.window_function); // blackman je bolji od hann
    let sample_notes = freq_generator::open_sample_notes(args.nfft);
    let mut all_notes = harmonics::generate_all_notes();
    harmonics::generate_note_network(&mut all_notes, &sample_notes, &window_fn, args.nfft);
    harmonics::cross_polinate(&mut all_notes, &sample_notes, &window_fn, args.nfft);

    let mut data = cached_process_song(&song, &args, &all_notes, &h, &window_fn, &sample_notes);

    let builder = gtk::Builder::from_string(GLADE_UI_SOURCE);
    let window = builder.object::<gtk::Window>("MainWindow").unwrap();
    let main_box = builder.object::<gtk::Box>("MainBox").unwrap();

    window.set_title("Generator tablatura");

    let drawing_area: gtk::DrawingArea = builder.object("MainDrawingArea").unwrap();
    // Kreiranje instanci slidera iz ui.glade fajla
    let low_threshold_slider = builder.object::<gtk::Scale>("LowThresholdSlider").unwrap();
    let high_threshold_slider = builder.object::<gtk::Scale>("HighThresholdSlider").unwrap();
    let string_slider = builder.object::<gtk::Scale>("StringSlider").unwrap();
    let time_frame_slider = builder.object::<gtk::Scale>("TimeFrameSlider").unwrap();
    let time_slider = builder.object::<gtk::Scale>("TimeSlider").unwrap();
    let y_scale_slider = builder.object::<gtk::Scale>("YScaleSlider").unwrap();
    let attenuation_factor_slider = builder
        .object::<gtk::Scale>("AttenuationFactorSlider")
        .unwrap();
    let attenuation_power_factor_slider = builder
        .object::<gtk::Scale>("AttenuationPowerFactorSlider")
        .unwrap();

    let app_state = Rc::new(RefCell::new(PlottingState {
        attenuation_power_factor: attenuation_power_factor_slider.value(),
        attenuation_factor: attenuation_factor_slider.value(),
        time_processed: args.sec_to_run as f64,
        current_string: string_slider.value(),
        time_frame: time_frame_slider.value(),
        time_offset: time_slider.value(),
        y_slider: y_scale_slider.value(),
        data_att: data.clone(),
        data_out: data.clone(),
        high_threshold: 1.0,
        low_threshold: 0.0,
        data_orig: data,
        sample_notes,
        window_fn,
        all_notes,
        song,
        args,
        h,
    }));

    window.set_application(Some(app));

    let state_cloned = app_state.clone();

    drawing_area.connect_draw(move |widget, cr| {
        let w = widget.allocated_width() as u32;
        let h = widget.allocated_height() as u32;
        let state = state_cloned.borrow().clone();
        let backend = CairoBackend::new(cr, (w, h)).unwrap();
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

    let att =
        |state: &mut PlottingState,
         ploting_state: Box<dyn Fn(&mut PlottingState) -> &mut PlottingState + 'static>| {
            let data_att = harmonics::attenuate_harmonics(
                &state.data_orig,
                &state.all_notes,
                state.attenuation_factor as f32,
                state.attenuation_power_factor as f32,
            );

            ploting_state(&mut *state).data_att = data_att;
            ploting_state(&mut *state).data_out = post_processing::schmitt(
                &state.data_att,
                state.high_threshold as f32,
                state.low_threshold as f32,
            );
        };

    let handle_key =
        |window: &gtk::Window,
         ploting_state: Box<dyn Fn(&mut PlottingState) -> &mut PlottingState + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();

            window.connect_key_press_event(move |a, b| {
                let keyval = b.keycode().unwrap();
                let mut state = app_state.borrow_mut();

                ploting_state(&mut *state).current_string = match keyval {
                    111 => max(state.current_string as i32 - 1, 0) as f64, // up arrow key pressed
                    116 => min(state.current_string as i32 + 1, 5) as f64, // down arrow key
                    _ => state.current_string,
                };

                if ![113, 114].contains(&keyval) {
                    if [111, 116].contains(&keyval) {
                        drawing_area.queue_draw();
                    }
                    return Inhibit(true);
                }

                ploting_state(&mut *state).time_offset = match keyval {
                    113 => max(state.time_offset as i32 - 1, 0) as f64, // left arrow key pressed
                    114 => min(
                        state.time_offset as i32 + 1,
                        state.song.len() as i32 / 44100,
                    ) as f64, // right arrow key pressed
                    _ => state.time_offset,
                };

                if (state.time_offset > state.time_processed - state.time_frame) {
                    let now = Instant::now();

                    ploting_state(&mut *state).args.seek_offset = state.time_processed as f32;
                    ploting_state(&mut *state).time_processed += state.args.sec_to_run as f64;

                    let new_chunk = cached_process_song(
                        &state.song,
                        &state.args,
                        &state.all_notes,
                        &state.h,
                        &state.window_fn,
                        &state.sample_notes,
                    );

                    weld_note_intensity(&mut ploting_state(&mut *state).data_orig, new_chunk);
                    att(&mut state, Box::new(|s| s));

                    println!("Processing took {} seconds", now.elapsed().as_secs_f32());
                }

                drawing_area.queue_draw();
                Inhibit(true)
            });
        };

    let recalculate_attenuation =
        |slider: &gtk::Scale,
         ploting_state: Box<dyn Fn(&mut PlottingState) -> &mut PlottingState + 'static>,
         factor: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();
            slider.connect_value_changed(move |target| {
                let mut state = app_state.borrow_mut();
                *factor(&mut *state) = target.value();
                att(&mut state, Box::new(|s| s));
                drawing_area.queue_draw();
            });
        };

    let recalculate_thresholds =
        |what: &gtk::Scale,
         how: Box<dyn Fn(&mut PlottingState) -> &mut Vec<Vec<Vec<f32>>> + 'static>,
         factor: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();

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

    handle_key(&window, Box::new(|s| s));
    handle_change(&time_frame_slider, Box::new(|s| &mut s.time_frame));
    handle_change(&y_scale_slider, Box::new(|s| &mut s.y_slider));

    recalculate_attenuation(
        &attenuation_factor_slider,
        Box::new(|s| s),
        Box::new(|s| &mut s.attenuation_factor),
    );
    recalculate_attenuation(
        &attenuation_power_factor_slider,
        Box::new(|s| s),
        Box::new(|s| &mut s.attenuation_power_factor),
    );

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

fn cached_process_song(
    song: &Vec<i16>,
    args: &cli::Args,
    all_notes: &Vec<Note>,
    h: &Vec<f32>,
    window: &Vec<f32>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
) -> Vec<Vec<Vec<f32>>> {
    let song = fetch_song(song, args.seek_offset, args.sec_to_run);
    //let mut note_intensity =
    fft::threaded_interlaced_convolution_realfft(&song, sample_notes, &window, &args)
    //   let mut note_intensity = post_processing::threaded_fft_fir_filtering(note_intensity, &h, &args);
    //attenuate_neighbours(&note_intensity, &all_notes, 0.75, 1.0)
}

fn weld_note_intensity(first_part: &mut Vec<Vec<Vec<f32>>>, second_part: Vec<Vec<Vec<f32>>>) {
    for wire in 0..first_part.len() {
        for note in 0..first_part[0].len() {
            first_part[wire][note].extend(second_part[wire][note].to_vec());
        }
    }
}

fn fetch_song(song: &Vec<i16>, seek: f32, seconds: f32) -> Vec<i16> {
    let song_samples = (seconds * 44100.0) as usize;
    let seek_samples = (seek * 44100.0) as usize;

    let end_sample = std::cmp::min(seek_samples + song_samples, song.len());
    let seek_samples = std::cmp::min(seek_samples, song.len());
    song[seek_samples..end_sample].to_vec()
}
