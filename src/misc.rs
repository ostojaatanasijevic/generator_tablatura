pub fn hertz_to_notes() -> Vec<String>{
    let mut notes = vec![String::new();43];
  
    for i in 0..43{
        if i < 20{
            notes[i].push_str("E");
            notes[i].push_str(&i.to_string());
        }
        if i > 4 && i < 25{
            notes[i].push_str(":A");
            notes[i].push_str(&(i - 5 ).to_string());
        }
        if i > 9 && i < 30{
            notes[i].push_str(":D");
            notes[i].push_str(&(i - 10 ).to_string());
        }
        if i > 14 && i < 35{
            notes[i].push_str(":G");
            notes[i].push_str(&(i - 15).to_string());
        }
        if i > 18 && i < 39{
            notes[i].push_str(":B");
            notes[i].push_str(&(i - 19).to_string());
        }
        if i > 23{
            notes[i].push_str(":e");
            notes[i].push_str(&(i - 24).to_string());
        }

        if notes[i].chars().nth(0).unwrap() == ':' { notes[i].remove(0); }
    }

    for i in 0..43{
        println!("{}", notes[i]);
    }

    notes
}

pub fn index_to_note(i: usize) -> String{
    let mut note = String::new(); 

        if i < 20{
            note.push_str("E");
            note.push_str(&i.to_string());
        }
        if i > 4 && i < 25{
            note.push_str(":A");
            note.push_str(&(i - 5 ).to_string());
        }
        if i > 9 && i < 30{
            note.push_str(":D");
            note.push_str(&(i - 10 ).to_string());
        }
        if i > 14 && i < 35{
            note.push_str(":G");
            note.push_str(&(i - 15).to_string());
        }
        if i > 18 && i < 39{
            note.push_str(":B");
            note.push_str(&(i - 19).to_string());
        }
        if i > 23{
            note.push_str(":e");
            note.push_str(&(i - 24).to_string());
        }

        note
}


