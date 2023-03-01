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
