use std::{
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

use flate2::{Compression, write::GzEncoder};

static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

#[cfg(windows)]
const HOME_ENV: &str = "USERPROFILE";
#[cfg(not(windows))]
const HOME_ENV: &str = "HOME";

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new() -> Self {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "seqsource-manifest-{}-{timestamp}-{}",
            std::process::id(),
            NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self { path }
    }

    fn file(&self, name: &str, contents: &[u8]) -> PathBuf {
        let path = self.path.join(name);
        fs::write(&path, contents).unwrap();
        path
    }

    fn gzip(&self, name: &str, contents: &[u8]) -> PathBuf {
        let path = self.path.join(name);
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        encoder.write_all(contents).unwrap();
        encoder.finish().unwrap();
        path
    }

    fn manifest(&self, headers: &[&str], rows: &[&[&str]]) -> PathBuf {
        let path = self.path.join("manifest.tsv");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .flexible(true)
            .from_path(&path)
            .unwrap();
        writer.write_record(headers).unwrap();
        for row in rows {
            writer.write_record(*row).unwrap();
        }
        writer.flush().unwrap();
        path
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

fn run_manifest(path: &Path, options: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_seqsource"))
        .arg("manifest")
        .args(options)
        .arg(path)
        .output()
        .unwrap()
}

fn records(bytes: &[u8]) -> (csv::StringRecord, Vec<csv::StringRecord>) {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(bytes);
    let header = reader.headers().unwrap().clone();
    let rows = reader.records().map(Result::unwrap).collect();
    (header, rows)
}

#[test]
fn manifest_preserves_order_and_reports_each_fastq_outcome() {
    let directory = TempDir::new();
    let known = directory.file("known.fastq", b"@A00119:1:FLOWCELL\nACGT\n+\n!!!!\n");
    let unknown = directory.gzip("unknown.fastq.gz", b"@ZZ12345:1:FLOWCELL\nACGT\n+\n!!!!\n");
    directory.file("empty.fastq", b"");
    directory.gzip("empty.fastq.gz", b"");
    directory.file("corrupt.fastq.gz", b"not gzip data");
    directory.file("bad_header.fastq", b"A00119:1:FLOWCELL\n");
    fs::create_dir(directory.path.join("a_directory.fastq")).unwrap();

    let manifest = directory.manifest(
        &["extra", "fastq", "id"],
        &[
            &["ignored", "known.fastq", "known"],
            &["ignored", "missing.fastq", "missing"],
            &["ignored", "unknown.fastq.gz", "unknown"],
            &["ignored", "empty.fastq", "empty_plain"],
            &["ignored", "empty.fastq.gz", "empty_gzip"],
            &["ignored", "corrupt.fastq.gz", "corrupt"],
            &["ignored", "bad_header.fastq", "bad_header"],
            &["ignored", "a_directory.fastq", "directory"],
        ],
    );
    let output = run_manifest(&manifest, &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let (header, rows) = records(&output.stdout);
    assert_eq!(
        header,
        csv::StringRecord::from(vec![
            "id",
            "fastq",
            "instrument_guess",
            "code",
            "fastq_header",
            "fastq_check_success",
            "failure_reason"
        ])
    );
    assert_eq!(rows.len(), 8);
    assert_eq!(
        rows.iter()
            .map(|row| row.get(0).unwrap())
            .collect::<Vec<_>>(),
        [
            "known",
            "missing",
            "unknown",
            "empty_plain",
            "empty_gzip",
            "corrupt",
            "bad_header",
            "directory"
        ]
    );
    assert_eq!(rows[0].get(1), Some("known.fastq"));
    assert_eq!(rows[0].get(5), Some("true"));
    assert_eq!(rows[2].get(2), Some("unknown"));
    assert_eq!(rows[2].get(3), Some("unknown"));
    assert_eq!(rows[2].get(5), Some("true"));
    assert_eq!(rows[2].get(6), Some(""));

    for (index, path) in [(0, known), (2, unknown)] {
        let fastq_output = Command::new(env!("CARGO_BIN_EXE_seqsource"))
            .current_dir(&directory.path)
            .arg("fastq")
            .arg("--id")
            .arg(rows[index].get(0).unwrap())
            .arg(path.file_name().unwrap())
            .output()
            .unwrap();
        assert!(fastq_output.status.success());
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(fastq_output.stdout.as_slice());
        assert_eq!(reader.headers().unwrap(), &header);
        let fastq_row = reader.records().next().unwrap().unwrap();
        assert_eq!(fastq_row, rows[index]);
    }

    for (index, expected) in [
        (1, "could not access FASTQ file"),
        (3, "FASTQ file is empty"),
        (4, "FASTQ file is empty"),
        (5, "could not read first FASTQ header"),
        (6, "first FASTQ line must start"),
        (7, "not a regular file"),
    ] {
        assert_eq!(rows[index].get(2), Some(""));
        assert_eq!(rows[index].get(3), Some(""));
        assert_eq!(rows[index].get(4), Some(""));
        assert_eq!(rows[index].get(5), Some("false"));
        assert!(rows[index].get(6).unwrap().contains(expected));
    }
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("Samples: total=8 successful=2 failed=6"));
    assert!(stderr.contains(
        "missing_file=1 empty_fastq=2 invalid_fastq=2 non_file_path=1 other_read_error=0"
    ));
    assert!(stderr.contains("Elapsed: "));

    for options in [&["-t", "1"][..], &["--threads", "4"][..]] {
        let threaded = run_manifest(&manifest, options);
        assert!(threaded.status.success());
        assert_eq!(threaded.stdout, output.stdout);
        let threaded_stderr = String::from_utf8(threaded.stderr).unwrap();
        assert_eq!(
            threaded_stderr.lines().take(2).collect::<Vec<_>>(),
            stderr.lines().take(2).collect::<Vec<_>>()
        );
    }
}

#[test]
fn manifest_expands_home_paths_and_preserves_input_paths() {
    let directory = TempDir::new();
    let home = TempDir::new();
    let reads = home.path.join("reads");
    fs::create_dir(&reads).unwrap();
    fs::write(reads.join("home.fastq"), b"@A00119:1\n").unwrap();
    directory.file("relative.fastq", b"@A00119:1\n");
    let absolute = directory.file("absolute.fastq", b"@A00119:1\n");
    directory.file("~otheruser.fastq", b"@A00119:1\n");

    let absolute_text = absolute.to_str().unwrap();
    let entries = vec![
        vec!["home", "~/reads/home.fastq"],
        vec!["relative", "relative.fastq"],
        vec!["absolute", absolute_text],
        vec!["other_user_literal", "~otheruser.fastq"],
    ];
    #[cfg(windows)]
    let entries = {
        let mut entries = entries;
        entries.push(vec!["home_backslash", r"~\reads\home.fastq"]);
        entries
    };
    let refs: Vec<&[&str]> = entries.iter().map(Vec::as_slice).collect();
    let manifest = directory.manifest(&["id", "fastq"], &refs);

    let output = Command::new(env!("CARGO_BIN_EXE_seqsource"))
        .arg("manifest")
        .arg(&manifest)
        .env(HOME_ENV, &home.path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (_, rows) = records(&output.stdout);
    assert_eq!(rows.len(), entries.len());
    for (row, entry) in rows.iter().zip(&entries) {
        assert_eq!(row.get(0), Some(entry[0]));
        assert_eq!(row.get(1), Some(entry[1]));
        assert_eq!(row.get(5), Some("true"));
    }
}

#[test]
fn manifest_reports_unavailable_or_invalid_home_before_output() {
    let directory = TempDir::new();
    let manifest = directory.manifest(&["id", "fastq"], &[&["sample", "~/reads.fastq"]]);

    for (home_value, expected) in [
        (None, "unset or empty"),
        (Some("relative"), "absolute path"),
    ] {
        let mut command = Command::new(env!("CARGO_BIN_EXE_seqsource"));
        command.arg("manifest").arg(&manifest);
        match home_value {
            Some(value) => {
                command.env(HOME_ENV, value);
            }
            None => {
                command.env_remove(HOME_ENV);
            }
        }
        let output = command.output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("manifest record 2"), "{stderr}");
        assert!(stderr.contains("~/reads.fastq"), "{stderr}");
        assert!(stderr.contains(HOME_ENV), "{stderr}");
        assert!(stderr.contains(expected), "{stderr}");
    }
}

#[test]
fn manifest_quotes_input_and_output_fields() {
    let directory = TempDir::new();
    let filename = "name\twith\"quote.fastq";
    directory.file(filename, b"@A00119:1\tsample \"one\"\n");
    let manifest = directory.manifest(&["id", "fastq"], &[&["id\t\"one\"", filename]]);
    let output = run_manifest(&manifest, &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (_, rows) = records(&output.stdout);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].get(0), Some("id\t\"one\""));
    assert_eq!(rows[0].get(1), Some(filename));
    assert_eq!(rows[0].get(4), Some("@A00119:1\tsample \"one\""));
    assert_eq!(rows[0].get(5), Some("true"));
}

#[test]
fn malformed_manifest_fails_before_writing_stdout() {
    let directory = TempDir::new();
    for (headers, rows) in [
        (vec!["id"], vec![vec!["sample"]]),
        (vec!["id", "id", "fastq"], vec![vec!["a", "b", "x.fastq"]]),
        (vec!["id", "fastq"], vec![vec!["", "x.fastq"]]),
        (vec!["id", "fastq"], vec![vec!["sample", ""]]),
        (vec!["id", "fastq"], vec![vec!["sample"]]),
        (
            vec!["id", "fastq"],
            vec![vec!["first", "missing.fastq"], vec!["", "x.fastq"]],
        ),
    ] {
        let refs: Vec<&[&str]> = rows.iter().map(Vec::as_slice).collect();
        let manifest = directory.manifest(&headers, &refs);
        let output = run_manifest(&manifest, &[]);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!output.stderr.is_empty());
    }
}

#[test]
fn header_only_manifest_has_zero_counts() {
    let directory = TempDir::new();
    let manifest = directory.manifest(&["id", "fastq"], &[]);
    let output = run_manifest(&manifest, &[]);
    assert!(output.status.success());
    let (_, rows) = records(&output.stdout);
    assert!(rows.is_empty());
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("Samples: total=0 successful=0 failed=0")
    );
    let threaded = run_manifest(&manifest, &["-t", "4"]);
    assert!(threaded.status.success());
    assert_eq!(threaded.stdout, output.stdout);
    assert!(
        String::from_utf8(threaded.stderr)
            .unwrap()
            .contains("Samples: total=0 successful=0 failed=0")
    );
}

#[test]
fn invalid_thread_counts_are_rejected_before_output() {
    let directory = TempDir::new();
    let manifest = directory.manifest(&["id", "fastq"], &[]);
    for count in ["0", "not-a-number", "-1"] {
        let output = run_manifest(&manifest, &["-t", count]);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!output.stderr.is_empty());
    }
}

#[test]
fn manifest_help_describes_the_command_and_argument() {
    let top_level = Command::new(env!("CARGO_BIN_EXE_seqsource"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(top_level.status.success());
    assert!(
        String::from_utf8(top_level.stdout)
            .unwrap()
            .contains("Analyse FASTQ paths from a TSV manifest")
    );

    let manifest_help = Command::new(env!("CARGO_BIN_EXE_seqsource"))
        .args(["manifest", "--help"])
        .output()
        .unwrap();
    assert!(manifest_help.status.success());
    let help = String::from_utf8(manifest_help.stdout).unwrap();
    assert!(help.contains("TSV file with id and fastq columns"));
    assert!(help.contains("-t, --threads <N>"));
    assert!(help.contains("Number of FASTQ files to check at once"));
    assert!(help.contains("[default: 1]"));
}
