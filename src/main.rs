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

mod fourier;
mod legacy_fourier;
mod misc;
mod plot;
mod post_processing;

use fast_float::parse;
use misc::hertz_to_notes;
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

#[derive(Debug)]
pub struct NotePeak {
    time: f32,
    ampl: f32,
    index: usize,
}

pub const THREADS: usize = 8;
pub const SAMPLE: usize = 4096 * 2;
pub const F_RES: f32 = 44100.0 / SAMPLE as f32;
pub const T_RES: f32 = 1.0 / (44100.0);
pub const NFFT: usize = SAMPLE;
pub const AVG_LEN: usize = SAMPLE / 8; // mora da može da deli NFFT, da ne bi cureli podaci
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

fn main() {
    let args: Vec<String> = env::args().collect();
    let song = fourier::open_song(&args[1], 15.0);
    // blackman je bolji od hann
    let window = fourier::calculate_window_function(SAMPLE, "blackman");

    let mut all_notes = generate_all_notes();
    generate_note_network(&mut all_notes, &window);

    //ADD THREAD DETECTION FOR INDIVIDUAL CPUs
    let sample_ffts = legacy_fourier::calculate_sample_ffts(&window, SAMPLE, 0);
    let mut note_intensity =
        legacy_fourier::threaded_dtft_and_conv(&song, &sample_ffts, &window, "circular");

    plot::plot_data_norm(&note_intensity, "before_");
    let note_intensity = attenuate_harmonics(&note_intensity, &all_notes, 1.0);
    plot::plot_data_norm(&note_intensity, "after_");
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

fn generate_note_network(all_notes: &mut Vec<Note>, window: &Vec<f32>) {
    let end = all_notes.len();

    for note in 0..end {
        //calculate fft and detect peaks
        let data = fourier::open_sample_note(&all_notes[note]);
        let dft_data = fourier::dtft(&data, &window, SAMPLE)[0..SAMPLE / 16].to_vec();

        let peaks = post_processing::find_peaks(&dft_data, 0.041);

        //find start
        let mut base_ampl = 1.0;
        for peak in peaks.iter() {
            let peak_freq = peak.index as f32 * F_RES;
            if (peak_freq - all_notes[note].freq).abs() < F_RES * 1.5 {
                println!("found base harmonic of {}", all_notes[note].name);
                base_ampl = peak.ampl;
            }
        }

        if base_ampl == 1.0 {
            println!("failed to find {}", all_notes[note].name);
            println!("first peak: {}", peaks[0].index as f32 * F_RES);
        }

        //calculate ratios
        for peak in peaks.iter() {
            let peak_freq = peak.index as f32 * F_RES;
            for index in 0..120 {
                if (all_notes[index].freq - peak_freq).abs() > 5.0 {
                    continue;
                }

                let harmonic: i32 = (all_notes[index].freq / (all_notes[note].freq)).round() as i32;
                if harmonic == 1 {
                    continue;
                }

                let ratio = peak.ampl / base_ampl;
                all_notes[note].harmonics.push((index, ratio));
                println!(
                    "base: {}, harmonic: {}, {}th ratio: {}",
                    all_notes[note].name, all_notes[index].name, harmonic, ratio
                );
            }
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
