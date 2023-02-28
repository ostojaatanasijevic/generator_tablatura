#![allow(warnings, unused)]
// KORISTI REALFFT
// ADAPTIRAJ INTENZITET SIGURNOSTI U NOTU U POREDJENU SA INTENZITETIOM PESME

mod cli;
mod fft;
mod fourier;
mod freq_generator;
mod misc;
mod plot;
mod post_processing;

use std::cell::RefCell;
use std::rc::Rc;

use gtk::prelude::*;
use plotters::prelude::*;
use plotters_cairo::CairoBackend;

use clap::Parser;
use fast_float::parse;
use misc::hertz_to_notes;
use num::complex::ComplexFloat;
use num::traits::Pow;
use num::FromPrimitive;
use rayon::vec;
use rustfft::{num_complex::Complex, FftPlanner};
use std::cmp;
use std::env;
use std::error::Error;
use std::f32::consts::PI;
use std::fs::File;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::thread;
use std::time::Instant;

pub struct NotesToIndex {
    note: i32,
    index: i32,
}

pub struct Peaks {
    freq: f32,
    ampl: f32,
}

#[derive(Debug)]
pub struct NotePeak {
    pub time: f32,
    pub ampl: f32,
    pub index: usize,
}

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

// IZVORI ZABUNE
// 1. e0 = B5 = G9 itd
// rešenje: sistemom eliminacije odrediti koja nota zvoni na kojoj žici
// TEMPO SE RAČUNA AUTOKORELACIOJOM SIGNALA
//

#[derive(Debug, Clone)]
pub struct Note {
    name: String,
    freq: f32,
    harmonics: Vec<(usize, f32)>,
}

const offset_table: [usize; 6] = [0, 5, 10, 15, 19, 24];

fn open_sample_notes(sample_len: usize) -> Vec<Vec<Vec<i16>>> {
    let fs = 44100.0;
    let secs = sample_len as f32 / fs;
    let mut out = vec![vec![Vec::new(); 20]; 6];

    for string in 0..6 {
        for note in 0..20 {
            let freq = HERZ[offset_table[5 - string] + note]
                .parse::<f32>()
                .unwrap();

            out[string][note] = freq_generator::sin(freq, fs, secs); // 2.0
        }
    }

    out
}

fn generate_volume_map(song: &Vec<i16>, decemation_len: usize) -> Vec<f32> {
    // first abs
    let mut before: Vec<f32> = Vec::new();

    for sample in song.iter() {
        before.push(sample.abs() as f32);
    }

    let out = post_processing::block_max_decemation(&before, decemation_len);
    out
}

fn vector_multiply(a: &Vec<f32>, b: &Vec<f32>) -> Vec<f32> {
    let mut out: Vec<f32> = Vec::new();
    let len = cmp::min(a.len(), b.len());
    println!("first len {}, second {}", a.len(), b.len());

    for t in 0..len {
        out.push(a[t] * b[t]);
    }

    out
}

//ADD THREAD DETECTION FOR INDIVIDUAL CPUs
fn main() {
    let args = cli::Args::parse();
    let sample_len = args.nfft;

    let h = post_processing::lp_filter(args.w, args.lenght_fir);
    let sec_to_run: f32 = args.sec_to_run;
    let song = fourier::open_song(&args.file_name, sec_to_run);
    let window = fourier::calculate_window_function(sample_len, &args.window_function); // blackman je bolji od hann

    //let volume_map = generate_volume_map(&song, args.decemation_len);

    let sample_notes = open_sample_notes(sample_len);
    let mut all_notes = generate_all_notes();
    generate_note_network(&mut all_notes, &sample_notes, &window, sample_len);

    for note in all_notes.iter() {
        println!("Za notu: {}", note.name);
        for harm in note.harmonics.iter() {
            println!(
                "Harmonik {} je intenziteta {}",
                all_notes[harm.0].name, harm.1
            );
        }
    }

    println!("starting convolution...");
    let mut note_intensity = Vec::new();
    let mut handles = vec![];
    for i in 0..6 {
        let song = song.clone();
        let sample_notes = sample_notes[i].clone();
        let window = window.clone();
        let conv_type = args.conv_type.clone();

        handles.push(thread::spawn(move || {
            let mut notes_on_string = vec![Vec::new(); 20];
            for n in 0..20 {
                notes_on_string[n] =
                    fft::interlaced_convolution(&song, &sample_notes[n], &window, None, 512);
            }

            notes_on_string
        }));
    }

    let mut joined_data: Vec<f32> = Vec::new();
    for handle in handles {
        let tmp = handle.join().unwrap();
        note_intensity.push(tmp);
    }

    println!("convolution complete!");

    println!("threaded filtering started...");

    let window = fourier::calculate_window_function(h.len(), &args.window_function); // blackman je bolji od hann
    let mut handles = vec![];
    for i in 0..6 {
        let h = h.clone();
        let mut string_data = note_intensity.remove(0);
        let window = window.clone();
        let conv_type = args.conv_type.clone();

        handles.push(thread::spawn(move || {
            let mut out_string_data = vec![Vec::new(); 20];

            for n in 0..20 {
                let mut note_data = string_data.remove(0);
                //FFT method
                //za nepravilan broj padduj do 2_n
                out_string_data[n] = fft::convolve(&note_data, &h, &window, None); // applying low pass filter
                out_string_data[n] =
                    post_processing::block_max_decemation(&out_string_data[n], args.decemation_len);
            }

            (i, out_string_data)
        }));
    }

    let mut note_intensity = vec![Vec::new(); 6];
    for handle in handles {
        let tmp = handle.join().unwrap();
        note_intensity[tmp.0] = tmp.1;
    }

    println!("fir filters applied");
    /*
    // find the max on E0 - E4, that's a sure bet
    let mut temp = post_processing::eliminate_by_string(&note_intensity[5][0..5].to_vec());
    for n in 0..5 {
        note_intensity[5][n] = temp.remove(0);
    }
    // e15 - 20 as well
    let mut temp = post_processing::eliminate_by_string(&note_intensity[0][15..20].to_vec());
    for n in 15..20 {
        note_intensity[0][n] = temp.remove(0);
    }
    */
    let application = gtk::Application::new(
        Some("io.github.plotters-rs.plotters-gtk-demo"),
        Default::default(),
    );

    application.connect_activate(move |app| {
        build_ui(app, &note_intensity, args.sec_to_run, &all_notes);
    });

    application.run_with_args(&[""]);
}

fn generate_all_notes() -> Vec<Note> {
    let mut all_notes: Vec<Note> = Vec::new();
    let mut counter = 0;
    for string in STRINGS.iter().rev() {
        for n in 0..20 {
            let mut name_new = String::from(*string);
            name_new.push_str(&n.to_string());
            //println!("{name_new}");
            all_notes.push(Note {
                name: name_new,
                freq: HERZ[offset_table[counter] + n].parse().unwrap(),
                harmonics: Vec::new(),
            });
        }
        counter += 1;
    }

    all_notes
}

fn generate_note_network(
    all_notes: &mut Vec<Note>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
    window: &Vec<f32>,
    sample_len: usize,
) {
    let end = all_notes.len();

    for note in 0..end {
        //calculate fft and detect peaks
        let midi_note = fourier::open_sample_note(&all_notes[note]);

        let mut auto_conv = fft::convolve(
            &midi_note,
            &sample_notes[5 - note / 20][note % 20],
            &window,
            None,
        );

        // zero index compare
        //let mut baseline = auto_conv[0];
        // average compare
        let mut baseline = auto_conv.iter().sum::<f32>();
        // max based compare
        //let baseline = auto_conv.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
        //for each higher harmonic convolve and save ratio
        for h in 0..end {
            let ratio = all_notes[h].freq / all_notes[note].freq;
            let rounded_ratio = ratio.round();

            //ako nije celobrojni umnozak, šibaj dalje
            if (ratio - rounded_ratio).abs() > 0.03 || ratio < 1.5 {
                continue;
            }

            //get corresponding wav sine of freq harmonic

            let mut conv =
                fft::convolve(&midi_note, &sample_notes[5 - h / 20][h % 20], &window, None);
            // zero index compare
            //let mut intensity = conv[0];
            // average compare
            let mut intensity = conv.iter().sum::<f32>();
            // max based compare
            //let intensity = conv.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
            /*
            println!(
                "ratio of relation: {}_{} is {}",
                &all_notes[note].name,
                &all_notes[h].name,
                intensity / baseline
            );
            */

            all_notes[note].harmonics.push((h, intensity / baseline));
        }
    }
}

fn attenuate_harmonics(
    note_intensity: &Vec<Vec<Vec<f32>>>,
    all_notes: &Vec<Note>,
    factor: f32,
    power_of_harmonics: f32,
) -> Vec<Vec<Vec<f32>>> {
    let mut out = note_intensity.clone();

    // E A D G B e
    for note in 0..120 {
        for hb in all_notes[note].harmonics.iter() {
            let wire = 5 - hb.0 / 20;
            let tab = hb.0 % 20;

            for t in 0..out[wire][tab].len() {
                out[wire][tab][t] -=
                    factor * out[5 - note / 20][note % 20][t] * (hb.1).pow(power_of_harmonics);

                if out[wire][tab][t] < 0.0 {
                    out[wire][tab][t] = 0.0;
                }
            }

            println!(
                "Attenuating {} with {} and a ratio of {}, bfr",
                &all_notes[hb.0].name,
                &all_notes[note].name,
                hb.1 * factor
            );
        }
    }

    out
}

fn add_harmonics(
    note_intensity: &Vec<Vec<Vec<f32>>>,
    all_notes: &Vec<Note>,
    factor: f32,
    power_of_harmonics: f32,
) -> Vec<Vec<Vec<f32>>> {
    let mut out = note_intensity.clone();

    // E A D G B e
    for note in (0..120).rev() {
        for hb in all_notes[note].harmonics.iter() {
            let wire = 5 - note / 20;
            let tab = note % 20;

            for t in 0..out[wire][tab].len() {
                out[wire][tab][t] +=
                    factor * out[5 - hb.0 / 20][hb.0 % 20][t] * (hb.1).pow(power_of_harmonics);
            }

            println!(
                "Attenuating {} with {} and a ratio of {}, bfr",
                &all_notes[hb.0].name,
                &all_notes[note].name,
                hb.1 * factor
            );
        }
    }

    out
}
const GLADE_UI_SOURCE: &'static str = include_str!("ui.glade");

#[derive(Clone)]
struct PlottingState {
    x_offset: f64,
    y_scale: f64,
    std_x: f64,
    std_y: f64,
    pitch: f64,
    roll: f64,
    note: usize,
    string: usize,
    time_res: f32,
    data_orig: Box<Vec<Vec<Vec<f32>>>>,
    data_curr: Box<Vec<Vec<Vec<f32>>>>,
}

impl PlottingState {
    fn plot_pdf<'a, DB: DrawingBackend + 'a>(
        &self,
        backend: DB,
        data: &Vec<Vec<Vec<f32>>>,
        time: f32,
    ) -> Result<(), Box<dyn Error + 'a>> {
        let x_offset = (self.x_offset * 10.0) as usize;

        let end_point = cmp::min(
            x_offset + (self.roll as f32 / 20.0 * data[0][0].len() as f32) as usize,
            data[0][0].len(),
        );

        let mut max: f32 = 0.1;

        let root = backend.into_drawing_area();

        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .set_label_area_size(LabelAreaPosition::Bottom, 60)
            .build_cartesian_2d(0..end_point - x_offset, 0.0..max)?;

        println!("{}:{}:{}", x_offset, end_point, end_point - x_offset);

        chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .x_label_formatter(&|x| format!("{}", ((*x + x_offset) as f32 * self.time_res) as f32))
            .y_labels(0)
            .draw()?;

        for note in 0..data[0].len() {
            chart.draw_series(
                AreaSeries::new(
                    (0..)
                        .zip(data[self.pitch as usize][note][x_offset..end_point].iter())
                        .map(|(x, y)| (x, *y / 20.0 * self.y_scale as f32 + note as f32 / 200.0)),
                    note as f32 / 200.0,
                    &RED.mix(0.2),
                )
                .border_style(&RED),
            )?;

            root.draw(&Text::new(
                format!("{}{}", STRINGS[self.pitch as usize], note),
                (
                    10,
                    ((20 - note) as i32 * root.get_pixel_range().1.end as i32 / 21) as i32,
                ),
                ("sans-serif", 15.0).into_font(),
            ));
        }

        root.present()?;
        Ok(())
    }
}

fn build_ui(app: &gtk::Application, data: &Vec<Vec<Vec<f32>>>, sec: f32, all_notes: &Vec<Note>) {
    let builder = gtk::Builder::from_string(GLADE_UI_SOURCE);
    let window = builder.object::<gtk::Window>("MainWindow").unwrap();

    window.set_title("Generator tablatura");

    let drawing_area: gtk::DrawingArea = builder.object("MainDrawingArea").unwrap();
    let string_slider = builder.object::<gtk::Scale>("StringSlider").unwrap();
    let yaw_scale = builder.object::<gtk::Scale>("YawScale").unwrap();
    let time_slider = builder.object::<gtk::Scale>("TimeSlider").unwrap();
    let y_scale_slider = builder.object::<gtk::Scale>("MeanYScale").unwrap();
    let std_x_scale = builder.object::<gtk::Scale>("SDXScale").unwrap();
    let std_y_scale = builder.object::<gtk::Scale>("SDYScale").unwrap();

    let note_default = 2;
    let string_default = 0;

    let app_state = Rc::new(RefCell::new(PlottingState {
        x_offset: time_slider.value(),
        y_scale: y_scale_slider.value(),
        std_x: std_x_scale.value(),
        std_y: std_y_scale.value(),
        pitch: string_slider.value(),
        roll: yaw_scale.value(),
        note: note_default,
        string: string_default,
        time_res: sec / data.len() as f32,
        data_orig: Box::new(data.clone()),
        data_curr: Box::new(data.clone()),
    }));

    window.set_application(Some(app));

    // detect std_x
    let mut out = std_x_scale.connect_value_changed(move |target| println!("change"));

    let state_cloned = app_state.clone();
    let mut data = data.clone();
    let mut data_att = data.clone();
    drawing_area.connect_draw(move |widget, cr| {
        let state = state_cloned.borrow().clone();
        let w = widget.allocated_width();
        let h = widget.allocated_height();
        let backend = CairoBackend::new(cr, (w as u32, h as u32)).unwrap();
        state.plot_pdf(backend, &*state.data_curr, sec).unwrap();
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

    let recalculate_linear_attenuation = |what: &gtk::Scale,
                                          how: Box<
        dyn Fn(&mut PlottingState) -> &mut Box<Vec<Vec<Vec<f32>>>> + 'static,
    >,
                                          factor: Box<
        dyn Fn(&mut PlottingState) -> &mut f64 + 'static,
    >| {
        let app_state = app_state.clone();
        let drawing_area = drawing_area.clone();
        let all_notes = all_notes.clone();
        what.connect_value_changed(move |target| {
            let mut state = app_state.borrow_mut();
            *factor(&mut *state) = target.value();
            **how(&mut *state) =
                attenuate_harmonics(&data, &all_notes, target.value() as f32, state.std_y as f32);
            drawing_area.queue_draw();
        });
    };

    let recalculate_exponential_attenuation =
        |what: &gtk::Scale,
         how: Box<dyn Fn(&mut PlottingState) -> &mut Box<Vec<Vec<Vec<f32>>>> + 'static>,
         factor: Box<dyn Fn(&mut PlottingState) -> &mut f64 + 'static>| {
            let app_state = app_state.clone();
            let drawing_area = drawing_area.clone();
            let all_notes = all_notes.clone();
            what.connect_value_changed(move |target| {
                let mut state = app_state.borrow_mut();
                *factor(&mut *state) = target.value();
                **how(&mut *state) = attenuate_harmonics(
                    &state.data_orig,
                    &all_notes,
                    state.std_x as f32,
                    target.value() as f32,
                );
                //                data = data.clone();
                drawing_area.queue_draw();
            });
        };

    handle_change(&string_slider, Box::new(|s| &mut s.pitch));
    handle_change(&yaw_scale, Box::new(|s| &mut s.roll));
    handle_change(&time_slider, Box::new(|s| &mut s.x_offset));
    handle_change(&y_scale_slider, Box::new(|s| &mut s.y_scale));

    recalculate_linear_attenuation(
        &std_x_scale,
        Box::new(|s| &mut s.data_curr),
        Box::new(|s| &mut s.std_x),
    );
    recalculate_exponential_attenuation(
        &std_y_scale,
        Box::new(|s| &mut s.data_curr),
        Box::new(|s| &mut s.std_y),
    );

    window.show_all();
}
