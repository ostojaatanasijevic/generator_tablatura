pub fn hertz_to_notes() -> Vec<String> {
    let mut notes = vec![String::new(); 43];

    for i in 0..43 {
        if i > 23 {
            notes[i].push_str(":e");
            notes[i].push_str(&(i - 24).to_string());
        }
        if i > 18 && i < 39 {
            notes[i].push_str(":B");
            notes[i].push_str(&(i - 19).to_string());
        }
        if i > 14 && i < 35 {
            notes[i].push_str(":G");
            notes[i].push_str(&(i - 15).to_string());
        }
        if i > 9 && i < 30 {
            notes[i].push_str(":D");
            notes[i].push_str(&(i - 10).to_string());
        }
        if i > 4 && i < 25 {
            notes[i].push_str(":A");
            notes[i].push_str(&(i - 5).to_string());
        }
        if i < 20 {
            notes[i].push_str("E");
            notes[i].push_str(&i.to_string());
        }

        if notes[i].chars().nth(0).unwrap() == ':' {
            notes.remove(0);
        }
    }

    for i in 0..43 {
        println!("{}", notes[i]);
    }

    notes
}

pub fn index_to_note(i: usize) -> String {
    let mut notes = Vec::new();
    if i > 23 {
        notes.push(format!("e{}", (i - 24)));
    }
    if i > 18 && i < 39 {
        notes.push(format!("B{}", (i - 19)));
    }
    if i > 14 && i < 35 {
        notes.push(format!("G{}", (i - 15)));
    }
    if i > 9 && i < 30 {
        notes.push(format!("D{}", (i - 10)));
    }
    if i > 4 && i < 25 {
        notes.push(format!("A{}", (i - 5)));
    }
    if i < 20 {
        notes.push(format!("E{}", i));
    }

    let mut note = String::from(&notes[0]);
    for i in 1..notes.len() {
        note.push_str(":");
        note.push_str(&notes[i]);
    }

    note
}
