use crate::HERZ;
use crate::OFFSET_TABLE;
use std::f32::consts::PI;

pub fn sin(freq: f32, fs: f32, len_s: f32) -> Vec<i16> {
    let mut out: Vec<i16> = Vec::new();
    let n = (len_s * fs) as usize;
    for t in 0..n {
        let angle = 2.0 * PI * t as f32 * freq / fs;
        out.push((angle.sin() * 100000.0) as i16);
    }

    out
}

pub fn sin_coherant(freq: f32, fs: f32, len_s: f32) -> Vec<i16> {
    let mut out: Vec<i16> = Vec::new();
    let len_s = len_s * fs;
    let n_t = (len_s * freq / fs);
    let n_ratio: f32 = (len_s / n_t); //.floor().max(1.0);
    let n = (n_t * n_ratio) as usize;

    let mut n = (fs / freq);
    let ratio = (len_s / n) as usize;
    n *= ratio as f32;
    println!("{}:{}:{}", n_t, n_ratio, n);
    for t in 0..n as usize {
        let angle = 2.0 * PI * t as f32 * freq / fs;
        out.push((angle.sin() * 100000.0) as i16);
    }
    println!("len is : {}", out.len());
    out
}

pub fn open_sample_notes(sample_len: usize) -> Vec<Vec<Vec<i16>>> {
    let fs = 44100.0;
    let secs = sample_len as f32 / fs;
    let mut out = vec![vec![Vec::new(); 20]; 6];

    for string in 0..6 {
        for note in 0..20 {
            let freq = HERZ[OFFSET_TABLE[5 - string] + note]
                .parse::<f32>()
                .unwrap();

            out[string][note] = sin(freq, fs, secs); // 2.0
        }
    }

    out
}
