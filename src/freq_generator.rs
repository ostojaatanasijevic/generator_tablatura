pub const PI: f32 = 3.1415;

pub fn sin(freq: f32, fs: f32, len_s: f32) -> Vec<i16> {
    let mut out: Vec<i16> = Vec::new();
    let n = (len_s * fs) as usize;
    for t in 0..n {
        let angle = 2.0 * PI * t as f32 * freq / fs;
        out.push((angle.sin() * 100000.0) as i16);
    }

    out
}
