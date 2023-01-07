#![allow(warnings, unused)]

mod misc;
mod fourier;
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
pub const INTERPOL: usize = 1; 
pub const SAMPLE: usize = 8192;
pub const F_RES: f32 = 2.0 * 44100.0 / SAMPLE as f32;
pub const T_RES: f32 = 1.0 / (44100.0);
pub const NFFT: usize = SAMPLE;
pub const AVG_LEN: usize = SAMPLE / 4; // mora da može da deli NFFT, da ne bi cureli podaci
pub const STRINGS: [&str; 6] = ["e","B","G","D","A","E"];

pub const HERZ: [&str; 43] = [ "82.41", "87.31","92.5","98.0","103.83","110.0","116.54","123.47","130.81",
"138.59","146.83","155.56","164.81","174.61","185.0","196.0","207.65","220.0","233.08","246.94","261.63",
"277.18","293.66","311.13","329.63","349.23","369.99","392.0","415.3","440.0","466.16","493.88","523.25",
"554.37","587.33","622.25","659.26","698.46","739.99","830.61","880.0","932.33","987.77"];

//
//
// IZVORI ZABUNE
// 1. e0 = B5 = G9 itd
// rešenje: sistemom eliminacije odrediti koja nota zvoni na kojoj žici
// 2. nota E0 ima harmonike na E14, B0, itd
// rešenje: možda oduzeti grafike: B0 -= E0 * faktor_harmonika
//
//
struct Note{
    name: String,
    freq: f32,
    harmonics: Vec<Note>,
}

fn main(){
    //ADD THREAD DETECTION FOR INDIVIDUAL CPUs
    //let lookup = hertz_to_notes();
    let window = fourier::calculate_window_function(SAMPLE, "hann");
    let sample_ffts = fourier::calculate_sample_ffts(&window);

    let song_data = fourier::open_song("songs/januar.wav", 10.0);
    println!("broj samplova pesme je: {}", song_data.len());
    //MOŽDA TI NI NE TREBA KONVOLICIJA; MOŽDA JE DTFT DOVOLJAN ZA POSTIZANJE PRECIZNE PROCENE
    //FREKVENCIJE
    /* 
    let note_intensity = fourier::dtft(&song_data, &window);
    plot::single_plot_data_norm(note_intensity); 
    */ 
    // UZMI U OBZIR KAKO PROZORSKA FUNKCIJA UTIČE NA INTENZITET SIGNALA
    let note_intensity = fourier::threaded_dtft_and_conv(&song_data, &sample_ffts, &window, "circular");
    /* 
    let schmitt_out = post_processing::schmitt(&note_intensity, 8.0, 4.0);
    let diracs = post_processing::to_dirac(&schmitt_out); // rising edge detection
    
    for t in 0..diracs[0].len(){
        for string in 0..diracs.len(){
        for note in 0..diracs[0].len(){
            if diracs[string][note][t] == 1.0{
                println!("{} {}", t as f32 * T_RES * AVG_LEN as f32, misc::index_to_note(note));
            }
        }
        }
    }
    */
    // AKO SU NA ISTOJ ŽICI 2 NOTE onda PROCES eliminacije
    plot::plot_data_norm(note_intensity); 
    
    //parse_peaks(&note_intensity);
}
