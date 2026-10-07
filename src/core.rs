use regex::regex;

pub trait InstrumentPattern {
    fn name(&self) -> String;
    fn regex(&self) -> regex::Regex;
}
pub struct InstrumentPatterns {
    pub name: String,
    pub regex: regex::Regex,
}

impl InstrumentPatterns {
    pub fn is_match(&self, instrument_id: &str) -> bool {
        self.regex.is_match(instrument_id)
    }
}

pub fn instruments() -> Vec<InstrumentPatterns> {
    vec![
        InstrumentPatterns {
            name: "Illumina NovaSeq X Plus".to_string(),
            regex: regex!("^@LH[0-9]{5}").to_owned(),
        },
        InstrumentPatterns {
            name: "Illumina NovaSeq X".to_string(),
            regex: regex!("^@LL[0-9]{5}").to_owned(),
        },
        InstrumentPatterns {
            name: "Illumina NovaSeq 6000".to_string(),
            regex: regex!("^@A[0-9]{5}").to_owned(),
        },
        InstrumentPatterns {
            name: "Illumina NovaSeq 6000Dx".to_string(),
            regex: regex!("^@ADX[0-9]{5}").to_owned(),
        },
        InstrumentPatterns {
            name: "Illumina HiSeq X".to_string(),
            regex: regex!("^@ST-E[0-9]{5}").to_owned(),
        },
    ]
}

pub fn identify_instrument(fastq_header: &str) -> String {
    let instruments = instruments();
    identify_instrument_from_patterns(fastq_header, &instruments)
}

fn identify_instrument_from_patterns(
    fastq_header: &str,
    patterns: &[InstrumentPatterns],
) -> String {
    let hits: Vec<_> = patterns
        .iter()
        .filter(|pattern| pattern.is_match(fastq_header))
        .collect();

    match hits.len() {
        0 => "unknown".to_string(),
        1 => hits.first().expect("bug in extraction of instrument name. Please report this error message in a new github issue").name.clone(),
        _ => format!(
            "ambiguous:{}",
            hits.iter()
                .map(|instrument| instrument.name.as_str())
                .collect::<Vec<_>>()
                .join(";")
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::{
        InstrumentPatterns, identify_instrument, identify_instrument_from_patterns, instruments,
    };
    use std::io::ErrorKind;

    const CASES: [(&str, &str); 5] = [
        ("Illumina NovaSeq X Plus", "@LH12345"),
        ("Illumina NovaSeq X", "@LL12345"),
        ("Illumina NovaSeq 6000", "@A12345"),
        ("Illumina NovaSeq 6000Dx", "@ADX12345"),
        ("Illumina HiSeq X", "@ST-E12345"),
    ];

    #[test]
    fn matches_instrument_codes_at_start_of_fastq_headers() {
        let patterns = instruments();

        for (name, code) in CASES {
            let pattern = patterns
                .iter()
                .find(|pattern| pattern.name == name)
                .unwrap();
            assert!(pattern.is_match(code), "{name} should match {code}");

            let header = format!("{code}:1:FLOWCELL:1:1101:1000:1000");
            assert!(pattern.is_match(&header), "{name} should match {header}");
        }
    }

    #[test]
    fn rejects_other_instruments_and_leading_characters() {
        let patterns = instruments();

        for (name, code) in CASES {
            let header = format!("{code}:1:FLOWCELL:1:1101:1000:1000");
            for pattern in &patterns {
                assert_eq!(
                    pattern.is_match(&header),
                    pattern.name == name,
                    "unexpected match for {} against {header}",
                    pattern.name
                );
                assert!(
                    !pattern.is_match(&format!("prefix{header}")),
                    "{} should reject characters before @",
                    pattern.name
                );
            }
        }
    }

    #[test]
    fn rejects_malformed_instrument_codes() {
        let patterns = instruments();

        for (name, code) in CASES {
            let pattern = patterns
                .iter()
                .find(|pattern| pattern.name == name)
                .unwrap();
            let prefix = code.strip_suffix("12345").unwrap();

            for malformed in [format!("{prefix}1234"), format!("{prefix}12A45")] {
                assert!(
                    !pattern.is_match(&malformed),
                    "{name} should reject {malformed}"
                );
            }
        }
    }

    #[test]
    fn identifies_each_supported_instrument() {
        for (name, code) in CASES {
            let header = format!("{code}:1:FLOWCELL:1:1101:1000:1000");
            assert_eq!(identify_instrument(&header), name);
        }
    }

    #[test]
    fn returns_unknown_when_no_instrument_matches() {
        assert_eq!(identify_instrument("@ZZ12345:1:FLOWCELL"), "unknown");
    }

    #[test]
    fn returns_ambiguous_when_multiple_patterns_match() {
        let patterns = [
            InstrumentPatterns {
                name: "First".to_string(),
                regex: regex::Regex::new("^@LH[0-9]{5}").unwrap(),
            },
            InstrumentPatterns {
                name: "Second".to_string(),
                regex: regex::Regex::new("^@LH").unwrap(),
            },
        ];

        assert_eq!(
            identify_instrument_from_patterns("@LH12345:1:FLOWCELL", &patterns),
            "ambiguous:First;Second"
        );
    }

    #[test]
    fn lists_all_ambiguous_matches_in_pattern_order() {
        let patterns = [
            InstrumentPatterns {
                name: "First".to_string(),
                regex: regex::Regex::new("^@LH[0-9]{5}").unwrap(),
            },
            InstrumentPatterns {
                name: "Other".to_string(),
                regex: regex::Regex::new("^@LL").unwrap(),
            },
            InstrumentPatterns {
                name: "Second".to_string(),
                regex: regex::Regex::new("^@LH").unwrap(),
            },
            InstrumentPatterns {
                name: "First".to_string(),
                regex: regex::Regex::new("^@LH123").unwrap(),
            },
        ];

        assert_eq!(
            identify_instrument_from_patterns("@LH12345:1:FLOWCELL", &patterns),
            "ambiguous:First;Second;First"
        );
    }

    #[test]
    fn identify_instrument_works_for_common_instruments() {
        assert_eq!(identify_instrument("@A00119"), "Illumina NovaSeq 6000");
        assert_eq!(identify_instrument("@A00488"), "Illumina NovaSeq 6000");
        assert_eq!(identify_instrument("@ST-E00185"), "Illumina HiSeq X");
        assert_eq!(identify_instrument("@ST-E00180"), "Illumina HiSeq X");
        assert_eq!(identify_instrument("@A00119"), "Illumina NovaSeq 6000");
    }
}
