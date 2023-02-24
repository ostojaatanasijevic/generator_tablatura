#![allow(warnings, unused)]

// dodaj prefilter
// radi dtft za svaku frekvencija, time se može postići velika preciznost
//
// skroz nova fora al ne baš
// radi konvoluciju sa sinusima
// u takvom slučaju nema dvoznačnosti o tome odakle se pojavljuje harmonik, pre ili posle
// može samo posle
// tkd kreneš od E na dole i ubijaš harmonike proporcijonalno njihovim intezitetima za zadatu notu
// odakle intenziteti?
// odradiš fft za svaku notu i odmeriš peakove, mnogo prostije nego dosadašnja tehnika

// dobra ideja zasad debilno implementirana

mod fft;
mod fourier;
mod legacy_fourier;
mod misc;
mod plot;
mod post_processing;

use clap::Parser;
use fast_float::parse;
use misc::hertz_to_notes;
use num::complex::ComplexFloat;
use num::FromPrimitive;
use plotters::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
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

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
///Generator tablatura
struct Args {
    #[arg(short, long, default_value_t = 3072)]
    ///Lenght of fft sample size
    nfft: usize,
    #[arg(short, long, default_value_t = 0.1)]
    ///Filter cutoff frequency
    w: f32,
    ///Circular convolution type: add, save
    #[arg(short, long, default_value = "add")]
    conv_type: String,

    #[arg(short, long)]
    ///Song file path
    file_name: String,

    #[arg(short, long, default_value_t = 10.0)]
    ///Number of seconds to analyze
    sec_to_run: f32,

    #[arg(short, long, default_value = "blackman")]
    ///Window function
    window_function: String,
}

#[derive(Debug)]
pub struct NotePeak {
    time: f32,
    ampl: f32,
    index: usize,
}

pub const THREADS: usize = 8;
pub const SAMPLE: usize = 1024 * 3;
pub const F_RES: f32 = 44100.0 / SAMPLE as f32;
pub const NFFT: usize = SAMPLE;
pub const AVG_LEN: usize = 1; // SAMPLE / 32; // mora da može da deli NFFT, da ne bi cureli podaci
pub const T_RES: f32 = 1.0 / (44100.0) * AVG_LEN as f32;
pub const STRINGS: [&str; 6] = ["e", "B", "G", "D", "A", "E"];

pub const HERZ: [&str; 44] = [
    "82.41", "87.31", "92.5", "98.0", "103.83", "110.0", "116.54", "123.47", "130.81", "138.59",
    "146.83", "155.56", "164.81", "174.61", "185.0", "196.0", "207.65", "220.0", "233.08",
    "246.94", "261.63", "277.18", "293.66", "311.13", "329.63", "349.23", "369.99", "392.0",
    "415.3", "440.0", "466.16", "493.88", "523.25", "554.37", "587.33", "622.25", "659.26",
    "698.46", "739.99", "783.99", "830.61", "880.0", "932.33", "987.77",
];

//
//
// IZVORI ZABUNE
// 1. e0 = B5 = G9 itd
// rešenje: sistemom eliminacije odrediti koja nota zvoni na kojoj žici
// 2. nota E0 ima harmonike na E14, B0, itd
// rešenje: možda oduzeti grafike: B0 -= E0 * faktor_harmonika
//
//

// TEMPO SE RAČUNA AUTOKORELACIOJOM SIGNALA
//
// RAZLOG ZAŠTO CUSTOM sample_len nije rešilo u potpunosti problem šuma je prost
// za više note, šum najviše proizvode bass note, one kojima takodje odgovara
// novoodabrani sample_len, mada, rešili su se neki šumovi

// NAUČI MATI SOLO od Djordjeta Balaševića

// KORISTI REALFFT
#[derive(Debug, Clone)]
pub struct Note {
    name: String,
    freq: f32,
    harmonics: Vec<(usize, f32)>,
}

const offset_table: [usize; 6] = [0, 5, 10, 15, 19, 24];

//ADD THREAD DETECTION FOR INDIVIDUAL CPUs
fn main() {
    let args = Args::parse();

    let h = post_processing::lp_filter(args.w, 100);
    let sec_to_run: f32 = args.sec_to_run;
    let song = fourier::open_song(&args.file_name, sec_to_run);
    let window = fourier::calculate_window_function(SAMPLE, &args.window_function); // blackman je bolji od hann

    let sample_ffts = legacy_fourier::calculate_sample_ffts(&window, SAMPLE, SAMPLE);
    let mut all_notes = generate_all_notes();
    generate_note_network(&mut all_notes, &sample_ffts, &window);

    println!("starting convolution...");
    let mut note_intensity = Vec::new();
    let mut handles = vec![];
    for i in 0..6 {
        let song = song.clone();
        let sample_ffts = sample_ffts[i].clone();
        let window = window.clone();
        let conv_type = args.conv_type.clone();

        handles.push(thread::spawn(move || {
            let mut notes_on_string = vec![Vec::new(); 20];
            for n in 0..20 {
                notes_on_string[n] = fft::convolve(&song, &sample_ffts[n], &window, &conv_type);
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

    //  let mut note_intensity = legacy_fourier::threaded_dtft_and_conv(&song, &sample_ffts, &window, "add");
    for s in 0..6 {
        for n in 0..20 {
            post_processing::fir_filter(&h, &mut note_intensity[s][n]);
        }
    }

    println!("fir filters applied");

    plot::plot_data_norm(&note_intensity, "before_");
    let mut note_intensity = attenuate_harmonics(&note_intensity, &all_notes, 0.5);

    /*

    let peaks = post_processing::find_peaks(&note_intensity[0][0], 20000.0);
    for peak in peaks.iter() {
        println!("{:?}", peak);
    }
    */

    plot::plot_data_norm(&note_intensity, "after_");
    plot::draw_plot("plots/fir.png", h, 1.0, 1);
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

//modify to convolve
fn generate_note_network(
    all_notes: &mut Vec<Note>,
    sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
    window: &Vec<f32>,
) {
    let end = all_notes.len();

    for note in 0..end {
        //calculate fft and detect peaks
        let midi_note = fourier::open_sample_note(&all_notes[note]);

        let mut auto_conv = fft::convolve(
            &midi_note,
            &sample_ffts[5 - note / 20][note % 20],
            &window,
            "add",
        );

        let mut baseline = auto_conv.iter().sum::<f32>();

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
                fft::convolve(&midi_note, &sample_ffts[5 - h / 20][h % 20], &window, "add");

            let mut intensity = conv.iter().sum::<f32>();

            /*
            plot::draw_plot(
                &format!("plots/{}_{}.png", &all_notes[note].name, &all_notes[h].name),
                conv,
                1.0,
                1,
            );
            */
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
) -> Vec<Vec<Vec<f32>>> {
    let mut out = note_intensity.clone();

    // E A D G B e
    for note in 0..120 {
        for hb in all_notes[note].harmonics.iter() {
            let wire = 5 - hb.0 / 20;
            let tab = hb.0 % 20;

            for t in 0..note_intensity[0][0].len() {
                //                if (note_intensity[note / 20][note % 20][t] - note_intensity[wire][tab][t]).abs() < 200.0 {
                out[wire][tab][t] -= factor * note_intensity[5 - note / 20][note % 20][t] * hb.1;
                //               }
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
