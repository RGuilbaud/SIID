#CODE GENETIQUE
code_genet<-function(full_data,i){
  if(i%%3 == 0){
    if(full_data$minor_seq[i-2] == "T"){
      if(full_data$minor_seq[i-1] == "T"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "F"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "F"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "L"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "L"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "C"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "S"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "S"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "S"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "S"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "A"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "Y"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "Y"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "*"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "*"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "G"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "C"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "C"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "*"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "W"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "del") {
        full_data$minor[i] = "del"
      } else {
        full_data$minor[i] = "too_low"
      }
    } else if(full_data$minor_seq[i-2] == "C"){
      if(full_data$minor_seq[i-1] == "T"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "L"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "L"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "L"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "L"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "C"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "P"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "P"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "P"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "P"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "A"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "H"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "H"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "Q"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "Q"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "G"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "R"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "R"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "R"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "R"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "del") {
        full_data$minor[i] = "del"
      } else {
        full_data$minor[i] = "too_low"
      }
    } else if(full_data$minor_seq[i-2] == "A"){
      if(full_data$minor_seq[i-1] == "T"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "I"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "I"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "I"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "M"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "C"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "T"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "T"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "T"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "T"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "A"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "N"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "N"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "K"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "K"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "G"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "S"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "S"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "R"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "R"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "del") {
        full_data$minor[i] = "del"
      } else {
        full_data$minor[i] = "too_low"
      }
    } else if(full_data$minor_seq[i-2] == "G"){
      if(full_data$minor_seq[i-1] == "T"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "V"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "V"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "V"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "V"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "C"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "A"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "A"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "A"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "A"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "A"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "D"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "D"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "E"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "E"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "G"){
        if(full_data$minor_seq[i] == "T"){
          full_data$minor[i] = "G"
        } else if(full_data$minor_seq[i] == "C"){
          full_data$minor[i] = "G"
        } else if(full_data$minor_seq[i] == "A"){
          full_data$minor[i] = "G"
        } else if(full_data$minor_seq[i] == "G"){
          full_data$minor[i] = "G"
        } else if(full_data$minor_seq[i] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-1] == "del") {
        full_data$minor[i] = "del"
      } else {
        full_data$minor[i] = "too_low"
      }
    } else if(full_data$minor_seq[i-2] == "del") {
      full_data$minor[i] = "del"
    } else {
      full_data$minor[i] = "too_low"
    }
  } else {
    full_data$minor[i] = NA
  }
  
  if(i%%3 == 0){
    if(full_data$major_seq[i-2] == "T"){
      if(full_data$major_seq[i-1] == "T"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "F"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "F"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "L"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "L"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "C"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "S"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "S"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "S"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "S"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "A"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "Y"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "Y"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "*"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "*"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "G"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "C"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "C"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "*"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "W"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "del") {
        full_data$major[i] = "del"
      } else {
        full_data$major[i] = "too_low"
      }
    } else if(full_data$major_seq[i-2] == "C"){
      if(full_data$major_seq[i-1] == "T"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "L"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "L"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "L"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "L"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "C"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "P"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "P"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "P"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "P"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "A"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "H"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "H"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "Q"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "Q"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "G"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "R"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "R"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "R"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "R"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "del") {
        full_data$major[i] = "del"
      } else {
        full_data$major[i] = "too_low"
      }
    } else if(full_data$major_seq[i-2] == "A"){
      if(full_data$major_seq[i-1] == "T"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "I"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "I"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "I"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "M"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "C"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "T"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "T"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "T"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "T"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "A"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "N"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "N"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "K"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "K"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "G"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "S"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "S"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "R"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "R"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "del") {
        full_data$major[i] = "del"
      } else {
        full_data$major[i] = "too_low"
      }
    } else if(full_data$major_seq[i-2] == "G"){
      if(full_data$major_seq[i-1] == "T"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "V"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "V"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "V"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "V"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "C"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "A"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "A"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "A"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "A"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "A"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "D"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "D"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "E"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "E"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "G"){
        if(full_data$major_seq[i] == "T"){
          full_data$major[i] = "G"
        } else if(full_data$major_seq[i] == "C"){
          full_data$major[i] = "G"
        } else if(full_data$major_seq[i] == "A"){
          full_data$major[i] = "G"
        } else if(full_data$major_seq[i] == "G"){
          full_data$major[i] = "G"
        } else if(full_data$major_seq[i] == "del") {
          full_data$major[i] = "del"
        } else {
          full_data$major[i] = "too_low"
        }
      } else if(full_data$major_seq[i-1] == "del") {
        full_data$major[i] = "del"
      } else {
        full_data$major[i] = "too_low"
      }
    } else if(full_data$major_seq[i-2] == "del") {
      full_data$major[i] = "del"
    } else {
      full_data$major[i] = "too_low"
    }
  } else {
    full_data$major[i] = NA
  }
  return(full_data)
}