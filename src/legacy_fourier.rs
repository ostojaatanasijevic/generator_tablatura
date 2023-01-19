use fast_float::parse;
use plotters::prelude::*;
use rayon::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
use std::error::Error;
use std::f32::consts::PI;
use std::fs::File;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::thread;

use crate::offset_table;
use crate::post_processing::block_average_decemation;
use crate::post_processing::block_max_decemation;
use crate::Note;
use crate::AVG_LEN;
use crate::HERZ;
use crate::NFFT;
use crate::SAMPLE;
use crate::STRINGS;
use crate::THREADS;

const BROJ_ZICA: usize = 6;
const BROJ_PRAGOVA: usize = 20;

pub fn threaded_dtft_and_conv(
    song: &Vec<i16>,
    sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
    window: &Vec<f32>,
    convolution_type: &str,
) -> Vec<Vec<Vec<f32>>> {
    let mut chunks_of_the_song: Vec<Vec<i16>> = vec![];
    let mut chunk_lenght = (song.len() / THREADS) / SAMPLE;
    chunk_lenght = chunk_lenght * SAMPLE;

    for c in 0..THREADS {
        chunks_of_the_song.push(song[c * chunk_lenght..(c + 1) * chunk_lenght].to_vec());
    }

    let mut handles = vec![];
    for i in 0..THREADS {
        let gas = chunks_of_the_song[i].clone();
        let win = window.clone();
        let sam = sample_ffts.clone();
        let ct = String::from(convolution_type);
        handles.push(thread::spawn(move || {
            convolution_per_note(&gas, &sam, &win, &ct)
        }));
    }

    let mut joined_data: Vec<Vec<Vec<f32>>> = vec![vec![Vec::new(); BROJ_PRAGOVA]; BROJ_ZICA];
    for handle in handles {
        let tmp = handle.join().unwrap();

        for string in 0..tmp.len() {
            for note in 0..tmp[string].len() {
                joined_data[string][note].extend(&tmp[string][note]);
            }
        }
    }

    println!("threaded dtft and conv done!");
    joined_data
}

pub fn convolution_per_note(
    input_chunk: &Vec<i16>,
    sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
    window: &Vec<f32>,
    convolution_type: &str,
) -> Vec<Vec<Vec<f32>>> {
    let mut start_index: usize = 0;
    let mut stop_index: usize = SAMPLE;
    let mut nfft: usize = SAMPLE;

    if convolution_type == "circular" {
        start_index = 0;
        stop_index = SAMPLE;
        nfft = SAMPLE;
    } else if convolution_type == "save" {
        start_index = SAMPLE;
        stop_index = SAMPLE * 2;
        nfft = SAMPLE * 2;
    } else if convolution_type == "add" {
        start_index = 0;
        stop_index = SAMPLE * 2;
        nfft = SAMPLE * 2;
    }

    let mut final_buffer = vec![vec![Vec::<f32>::new(); BROJ_PRAGOVA]; BROJ_ZICA];
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let ifft = planner.plan_fft_inverse(nfft);

    let chunk_lenght = input_chunk.len();
    let num_of_chunks: usize = chunk_lenght / SAMPLE;

    for string in 0..sample_ffts.len() {
        for note in 0..sample_ffts[0].len() {
            let mut out = vec![0.0; input_chunk.len() + nfft];

            for c in 0..num_of_chunks {
                let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft];
                //NE IDE OVAKO ZA SAVE
                for i in 0..SAMPLE {
                    pesma_fft[i].re = (input_chunk[c * SAMPLE + i] as f32) * window[i] / 65536.0;
                }

                fft.process(&mut pesma_fft);

                let mut s_buffer: Vec<Complex<f32>> = pesma_fft
                    .iter()
                    .zip(sample_ffts[string][note].iter())
                    .map(|(x, y)| x * y.conj())
                    .collect();

                ifft.process(&mut s_buffer);

                let current: Vec<f32> = s_buffer[start_index..stop_index]
                    .iter()
                    .map(|a| a.norm())
                    .collect();

                for i in 0..current.len() {
                    out[c * SAMPLE + i] += current[i];
                }

                // decemate here; save RAM
                let decemeted =
                    block_max_decemation(&out[c * SAMPLE..(c + 1) * SAMPLE].to_vec(), AVG_LEN);
                final_buffer[string][note].extend(decemeted);
            }
        }
    }

    println!("dtft and conv done!");
    final_buffer
}

pub fn dtft_and_conv(
    input_chunk: &Vec<i16>,
    sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
    window: &Vec<f32>,
    convolution_type: &str,
) -> Vec<Vec<Vec<f32>>> {
    let mut start_index: usize = 0;
    let mut stop_index: usize = SAMPLE;
    let mut nfft: usize = SAMPLE;

    if convolution_type == "circular" {
        start_index = 0;
        stop_index = SAMPLE;
        nfft = SAMPLE;
    } else if convolution_type == "save" {
        start_index = SAMPLE;
        stop_index = SAMPLE * 2;
        nfft = SAMPLE * 2;
    } else if convolution_type == "add" {
        start_index = 0;
        stop_index = SAMPLE * 2;
        nfft = SAMPLE * 2;
    }

    let mut final_buffer = vec![vec![Vec::<f32>::new(); BROJ_PRAGOVA]; BROJ_ZICA];
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let ifft = planner.plan_fft_inverse(nfft);

    let mut sf = sample_ffts.clone();
    for string in 0..sample_ffts.len() {
        for note in 0..sample_ffts[0].len() {
            let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft - SAMPLE];
            sf[string][note].extend(pesma_fft);
        }
    }

    let chunk_lenght = input_chunk.len();
    let num_of_chunks: usize = chunk_lenght / SAMPLE;

    for c in 0..num_of_chunks {
        let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft];

        for i in 0..SAMPLE {
            pesma_fft[i].re = (input_chunk[c * SAMPLE + i] as f32) * window[i] / 65536.0;
        }

        fft.process(&mut pesma_fft);

        for string in 0..sample_ffts.len() {
            for note in 0..sample_ffts[0].len() {
                let mut s_buffer = vec![Complex { re: 0.0, im: 0.0 }; nfft];
                let mut out = vec![0.0; chunk_lenght + nfft * 2];

                s_buffer = pesma_fft
                    .iter()
                    .zip(sf[string][note].iter())
                    .map(|(x, y)| x * y.conj())
                    .collect();

                ifft.process(&mut s_buffer);

                let current: Vec<f32> = s_buffer[start_index..stop_index]
                    .iter()
                    .map(|a| a.norm())
                    .collect();

                for t in 0..current.len() {
                    out[c * SAMPLE + t] += current[t];
                }

                let decemeted =
                    block_max_decemation(&out[c * SAMPLE..(c + 1) * SAMPLE].to_vec(), SAMPLE / 4);
                // decemate here; save RAM
                final_buffer[string][note].extend(decemeted);
            }
        }
    }

    println!("dtft and conv done!");
    final_buffer
}

pub fn threaded_conv(
    fft_data: &Vec<Vec<Complex<f32>>>,
    sample_ffts: &Vec<Vec<Complex<f32>>>,
    percent: f32,
) -> Vec<Vec<Complex<f32>>> {
    let mut chunks_of_fft: Vec<Vec<Vec<Complex<f32>>>> = Vec::new();
    let mut chunk_lenght = ((fft_data.len() / THREADS) as f32 * percent) as usize;

    for c in 0..THREADS {
        chunks_of_fft.push(fft_data[c * chunk_lenght..(c + 1) * chunk_lenght].to_vec());
    }

    let mut handles = vec![];
    for i in 0..THREADS {
        let gas = chunks_of_fft[i].clone();
        let sam = sample_ffts.clone();
        handles.push(thread::spawn(move || {
            conv_with_samples(&gas, &sam, &(i as i32))
        }));
    }

    let mut joined_data: Vec<Vec<Complex<f32>>> = vec![Vec::new(); sample_ffts.len()];
    for handle in handles {
        let tmp = handle.join().unwrap();

        for note in 0..sample_ffts.len() {
            joined_data[note].extend(&tmp[note]);
        }
    }

    println!("threaded conv done!");
    joined_data
}

fn conv_with_samples(
    fft_data: &Vec<Vec<Complex<f32>>>,
    sample_ffts: &Vec<Vec<Complex<f32>>>,
    thread_index: &i32,
) -> Vec<Vec<Complex<f32>>> {
    //FFT DATA [time][0..SAMPLE]
    let mut final_buffer = vec![Vec::<Complex<f32>>::new(); sample_ffts.len()];
    let mut planner = FftPlanner::<f32>::new();
    let ifft = planner.plan_fft_inverse(SAMPLE);

    for t in 0..fft_data.len() {
        for note in 0..sample_ffts.len() {
            //let mut time = Instant::now();
            let mut s_buffer = vec![Complex { re: 0.0, im: 0.0 }; SAMPLE];

            s_buffer = fft_data[t]
                .iter()
                .zip(sample_ffts[note].iter())
                .map(|(x, y)| x * y.conj())
                .collect();

            ifft.process(&mut s_buffer);
            //MAYBE SLOW; PROLLY NOT
            final_buffer[note].extend(s_buffer);
            //println!("{:?}",time.elapsed());
        }
    }

    final_buffer
}

pub fn threaded_fourier(filename: &str, window: &Vec<f32>) -> Vec<Vec<Complex<f32>>> {
    let mut fft_mem: Vec<Vec<Complex<f32>>> = vec![Vec::new(); SAMPLE];
    let mut fajl_pesme = File::open(Path::new(filename)).unwrap();
    let (header, pesma) = wav::read(&mut fajl_pesme).unwrap();

    if !pesma.is_sixteen() {
        panic!("Wav file : {filename} isn't 16 bit! ");
    }
    let data = pesma.as_sixteen().unwrap();

    let mut chunks_of_the_song: Vec<Vec<i16>> = vec![];
    let mut chunk_lenght = data.len() / THREADS;

    for c in 0..THREADS {
        chunks_of_the_song.push(data[c * chunk_lenght..(c + 1) * chunk_lenght].to_vec());
    }

    let mut handles = vec![];
    for i in 0..THREADS {
        let gas = chunks_of_the_song[i].clone();
        let win = window.clone();
        handles.push(thread::spawn(move || fourier(gas, &win, &(i as i32))));
    }

    let mut joined_data: Vec<Vec<Complex<f32>>> = vec![];
    for handle in handles {
        let mut tmp = handle.join().unwrap();
        joined_data.extend(tmp);
    }

    println!("threaded fourier done!");
    joined_data
}

fn fourier(data: Vec<i16>, window: &Vec<f32>, thread_index: &i32) -> Vec<Vec<Complex<f32>>> {
    let mut fft_mem: Vec<Vec<Complex<f32>>> = vec![];
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(SAMPLE);

    for j in 0..data.len() / SAMPLE {
        let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; SAMPLE];

        for i in 0..SAMPLE {
            pesma_fft[i].re = (data[i + j * SAMPLE as usize] as f32) * window[i];
        }
        fft.process(&mut pesma_fft);
        fft_mem.push(pesma_fft);

        //let mut t = Instant::now();
        //println!("{:?}",t.elapsed());
    }

    fft_mem
}

pub fn calculate_sample_ffts(
    window: &Vec<f32>,
    sample_len: usize,
    padd: usize,
) -> Vec<Vec<Vec<Complex<f32>>>> {
    let start: usize = 0;
    let stop: usize = 20;
    let mut samples_fft = vec![vec![Vec::<Complex<f32>>::new(); stop - start]; 6];

    for string in 0..STRINGS.len() {
        for note in start..stop {
            //use format!
            let filename = format!("midi/{}/{}{}.wav", STRINGS[string], STRINGS[string], &note);
            let mut file = File::open(Path::new(&filename))
                .expect(&format!("Can't open file named {filename}"));
            let (_, raw_data) =
                wav::read(&mut file).expect(&format!("Can't read file, im retarded: ~{filename}~"));
            let data = raw_data
                .as_sixteen()
                .expect(&format!("Waw file: {} is not 16 bit", &filename));

            let mut planner = FftPlanner::<f32>::new();
            let fft = planner.plan_fft_forward(sample_len + padd);
            let mut fft_data = vec![Complex { re: 0.0, im: 0.0 }; sample_len + padd];

            for i in 0..sample_len {
                fft_data[i].re = (data[i] as f32) * window[i] / 65536.0;
            }

            fft.process(&mut fft_data);
            samples_fft[string][note] = fft_data; //PODELI SA SAMPLES
        }
    }

    println!("samples calculated!");
    samples_fft
}
