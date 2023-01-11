/**
*/


//
//
// IZVORI ZABUNE
// 1. e0 = B5 = G9 itd
// rešenje: sistemom eliminacije odrediti koja nota zvoni na kojoj žici
// 2. nota E0 ima harmonike na E14, B0, itd
// rešenje: možda oduzeti grafike: B0 -= E0 * faktor_harmonika
//
//


#[derive(Debug,Clone)]
pub struct Note{
    name: String,
    freq: f32,
    harmonics: Vec<usize>,
}

const offset_table: [usize; 6] = [0,5,10,15,19,24];



/*
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

*/
