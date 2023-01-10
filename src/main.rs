#![allow(warnings, unused)]

mod misc;
mod fourier;
mod legacy_fourier;
mod post_processing;
mod plot;

use std::io::Write;
use std::io::IoSlice;
use std::time::Instant;
use std::fs::File;
use std::thread;
use std::path::Path;
use std::error::Error;
use std::f32::consts::PI;
use rustfft::{FftPlanner, num_complex::Complex};
use plotters::prelude::*;
use fast_float::parse;
use misc::hertz_to_notes;
use num::{FromPrimitive};

pub struct NotesToIndex{
    note: i32,
    index: i32,
}

pub struct Peaks{
    freq: f32,
    ampl: f32,
}

#[derive(Debug)]
pub struct NotePeak{
    time: f32,
    ampl: f32,
    index: usize,
}

pub const THREADS: usize = 8; 
pub const SAMPLE: usize = 4096;
pub const F_RES: f32 = 2.0 * 44100.0 / SAMPLE as f32;
pub const T_RES: f32 = 1.0 / (44100.0);
pub const NFFT: usize = SAMPLE;
pub const AVG_LEN: usize = SAMPLE / 4; // mora da može da deli NFFT, da ne bi cureli podaci
pub const STRINGS: [&str; 6] = ["e","B","G","D","A","E"];

pub const HERZ: [&str; 44] = [ "82.41", "87.31","92.5","98.0","103.83","110.0","116.54","123.47","130.81",
"138.59","146.83","155.56","164.81","174.61","185.0","196.0","207.65","220.0","233.08","246.94","261.63",
"277.18","293.66","311.13","329.63","349.23","369.99","392.0","415.3","440.0","466.16","493.88","523.25",
"554.37","587.33","622.25","659.26","698.46","739.99","783.99","830.61","880.0","932.33","987.77"];

//
//
// IZVORI ZABUNE
// 1. e0 = B5 = G9 itd
// rešenje: sistemom eliminacije odrediti koja nota zvoni na kojoj žici
// 2. nota E0 ima harmonike na E14, B0, itd
// rešenje: možda oduzeti grafike: B0 -= E0 * faktor_harmonika
//
//
//

// TEMPO SE RAČUNA AUTOKORELACIOJOM SIGNALA
//
// RAZLOG ZAŠTO CUSTOM sample_len nije rešilo u potpunosti problem šuma je prost
// za više note, šum najviše proizvode bass note, one kojima takodje odgovara
// novoodabrani sample_len, mada, rešili su se neki šumovi

// KORISTI REALFFT 
#[derive(Debug,Clone)]
pub struct Note{
    name: String,
    freq: f32,
    harmonics: Vec<usize>,
}

const offset_table: [usize; 6] = [0,5,10,15,19,24];

fn main(){
    //ADD THREAD DETECTION FOR INDIVIDUAL CPUs

    let mut all_notes: Vec<Note> = Vec::new();
    let mut counter = 0;
    for string in STRINGS.iter().rev(){
        for n in 0..20{
            let mut name = String::from(*string);
            name.push_str(&n.to_string());
            println!("{name}");
            all_notes.push(Note{ name: name, freq: HERZ[offset_table[counter] + n].parse().unwrap(), harmonics: Vec::new()});
        }
        counter += 1;
    }

    for note in 0..all_notes.len(){
        for rev in note..all_notes.len(){
            let f1 = all_notes[rev].freq;
            let f2 = all_notes[note].freq;
            let r = f1 / f2;
            let m = r - r.floor();

            if (m < 0.1 && r > 1.5) || (m > 0.9 && r > 1.5){
                all_notes[note].harmonics.push(rev);
                println!("{:?} is {r}th harmonic of", all_notes[rev]);
                println!("{:?}\n", all_notes[note]);
            }
        }
    }

    let song_data = fourier::open_song("songs/januar.wav", 10.0);
    println!("broj samplova pesme je: {}", song_data.len());
 
    let song = fourier::open_song("songs/januar.wav", 10.0);
    let window = fourier::calculate_window_function(SAMPLE, "hann");
    let sample_ffts = legacy_fourier::calculate_sample_ffts(&window, SAMPLE);

    let mut note_intensity = legacy_fourier::threaded_dtft_and_conv(&song, &sample_ffts, &window, "circular");
    /* 
    for i in 0..20{
        for t in 0..note_intensity[0][0].len(){
            note_intensity[0][i][t] -= note_intensity[5][i][t] / 5.0;
        }
    }
    */

    plot::plot_data_norm(note_intensity);
}
