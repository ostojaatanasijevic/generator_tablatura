use crate::HERZ;
use crate::OFFSET_TABLE;
use crate::STRINGS;

use crate::BROJ_PRAGOVA;
use crate::BROJ_ZICA;

#[derive(Debug, Clone)]
pub struct Note {
    pub name: String,
    pub freq: f32,
    pub harmonics: Vec<(usize, f32)>,
    pub komsije: Vec<(usize, f32)>,
}

pub fn generate_all_notes() -> Vec<Note> {
    let mut all_notes: Vec<Note> = Vec::new();
    let mut counter = 0;
    for string in STRINGS.iter().rev() {
        for n in 0..BROJ_PRAGOVA {
            let mut name_new = String::from(*string);
            name_new.push_str(&n.to_string());
            //println!("{name_new}");
            all_notes.push(Note {
                name: name_new,
                freq: HERZ[OFFSET_TABLE[counter] + n].parse().unwrap(),
                harmonics: Vec::new(),
                komsije: Vec::new(),
            });
        }
        counter += 1;
    }

    all_notes
}

pub fn generate_note_network(
    all_notes: &mut Vec<Note>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
    window: &Vec<f32>,
    sample_len: usize,
) {
    let end = all_notes.len();

    for note in 0..end {
        //calculate fft and detect peaks
        let midi_note = crate::fourier::open_midi_note(&all_notes[note]);

        let mut auto_conv = crate::fft::convolve(
            &midi_note,
            &sample_notes[5 - note / BROJ_PRAGOVA][note % BROJ_PRAGOVA],
            &window,
            None,
        );

        // zero index compare
        //let mut baseline = auto_conv[0];
        // average compare
        let mut baseline = auto_conv.iter().map(|x| x.abs()).sum::<f32>();
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

            let mut conv = crate::fft::convolve(
                &midi_note,
                &sample_notes[5 - h / BROJ_PRAGOVA][h % BROJ_PRAGOVA],
                &window,
                None,
            );
            // zero index compare
            //let mut intensity = conv[0];
            // average compare
            let mut intensity = conv.iter().map(|x| x.abs()).sum::<f32>();
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

pub fn cross_polinate(
    all_notes: &mut Vec<Note>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
    window: &Vec<f32>,
    sample_len: usize,
) {
    let end = all_notes.len();

    for note in 0..end {
        let mut auto_conv = crate::fft::convolve(
            &sample_notes[5 - note / BROJ_PRAGOVA][note % BROJ_PRAGOVA],
            &sample_notes[5 - note / BROJ_PRAGOVA][note % BROJ_PRAGOVA],
            &window,
            None,
        );

        // zero index compare
        //let mut baseline = auto_conv[0];
        // average compare
        let mut baseline = auto_conv.iter().map(|x| x.abs()).sum::<f32>();
        // max based compare
        //let baseline = auto_conv.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
        //for each higher harmonic convolve and save ratio
        for h in 0..end {
            if all_notes[h].freq == all_notes[note].freq {
                continue;
            }
            let mut diff = all_notes[note].freq - all_notes[h].freq;
            diff = diff.abs();
            if diff > BROJ_PRAGOVA as f32 {
                continue;
            }

            let mut conv = crate::fft::convolve(
                &sample_notes[5 - note / BROJ_PRAGOVA][note % BROJ_PRAGOVA],
                &sample_notes[5 - h / BROJ_PRAGOVA][h % BROJ_PRAGOVA],
                &window,
                None,
            );
            // zero index compare
            //let mut intensity = conv[0];
            // average compare
            let mut intensity = conv.iter().map(|x| x.abs()).sum::<f32>();
            // max based compare
            //let intensity = conv.iter().max_by(|a, b| a.total_cmp(b)).unwrap();

            println!(
                "ratio of relation: {}_{} is {}",
                &all_notes[note].name,
                &all_notes[h].name,
                intensity / baseline
            );

            all_notes[note].komsije.push((h, intensity / baseline));
        }
    }
}

pub fn attenuate_neighbours(
    note_intensity: &Vec<Vec<Vec<f32>>>,
    all_notes: &Vec<Note>,
    factor: f32,
    power_of_harmonics: f32,
) -> Vec<Vec<Vec<f32>>> {
    let mut out = note_intensity.clone();

    // E A D G B e
    for note in 0..6 * BROJ_PRAGOVA {
        for hb in all_notes[note].komsije.iter() {
            let wire = 5 - hb.0 / BROJ_PRAGOVA;
            let tab = hb.0 % BROJ_PRAGOVA;

            for t in 0..out[wire][tab].len() {
                out[wire][tab][t] -= factor
                    * out[5 - note / BROJ_PRAGOVA][note % BROJ_PRAGOVA][t]
                    * (hb.1).powf(power_of_harmonics);

                if out[wire][tab][t] < 0.0 {
                    out[wire][tab][t] = 0.0;
                }
            }
        }
    }

    out
}

pub fn attenuate_harmonics(
    note_intensity: &Vec<Vec<Vec<f32>>>,
    all_notes: &Vec<Note>,
    factor: f32,
    power_of_harmonics: f32,
) -> Vec<Vec<Vec<f32>>> {
    let mut out = note_intensity.clone();

    // E A D G B e
    for note in 0..6 * BROJ_PRAGOVA {
        for hb in all_notes[note].harmonics.iter() {
            let wire = 5 - hb.0 / BROJ_PRAGOVA;
            let tab = hb.0 % BROJ_PRAGOVA;

            for t in 0..out[wire][tab].len() {
                out[wire][tab][t] -= factor
                    * out[5 - note / BROJ_PRAGOVA][note % BROJ_PRAGOVA][t]
                    * (hb.1).powf(power_of_harmonics);

                if out[wire][tab][t] < 0.0 {
                    out[wire][tab][t] = 0.0;
                }
            }

            /*
            println!(
                "Attenuating {} with {} and a ratio of {}, bfr",
                &all_notes[hb.0].name,
                &all_notes[note].name,
                hb.1 * factor
            );
            */
        }
    }

    out
}

pub fn add_harmonics(
    note_intensity: &Vec<Vec<Vec<f32>>>,
    all_notes: &Vec<Note>,
    factor: f32,
    power_of_harmonics: f32,
) -> Vec<Vec<Vec<f32>>> {
    let mut out = note_intensity.clone();

    // E A D G B e
    for note in (0..6 * BROJ_PRAGOVA).rev() {
        for hb in all_notes[note].harmonics.iter() {
            let wire = 5 - note / BROJ_PRAGOVA;
            let tab = note % BROJ_PRAGOVA;

            for t in 0..out[wire][tab].len() {
                out[wire][tab][t] += factor
                    * out[5 - hb.0 / BROJ_PRAGOVA][hb.0 % BROJ_PRAGOVA][t]
                    * (hb.1).powf(power_of_harmonics);
            }
            /*
            println!(
                "Attenuating {} with {} and a ratio of {}, bfr",
                &all_notes[hb.0].name,
                &all_notes[note].name,
                hb.1 * factor
            );
            */
        }
    }

    out
}
