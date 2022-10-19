//! This is a simple program that converts a BAM/CRAM file into lossy CRAM file with quality scores of Illumina 8-bins

// Author: Shohei Kojima @ RIKEN

const VERSION: &str = "version 0.0.1";
extern crate my_rust;
extern crate clap;
extern crate num_cpus;
extern crate rust_htslib;

// Illumina 8 binning accourding to https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf
// 0 = !, 6 = ', 15 = 0, 22 = 7, 27 = <, 33 = B, 37 = F, 40 = I
const TO_8_BINS: [u8; 110] = [
    0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 22, 22, 22, 22, 22, 27, 27, 27, 27, 27, 33, 33, 33, 33, 33, 37,
    37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40,
];

// 4 binning used in 1000 Genome Project and others, e.g. TOPMed https://topmed.nhlbi.nih.gov/sites/default/files/CORE_YR_1_WGS_CRAM_data_file_standardsv4.pdf
// 0 = !, 2 = #, 3 = $, 4 = %, 5 = &, 6 = ', 10 = +, 20 = 5, 30 = ?
const TO_4_BINS: [u8; 110] = [
    0, 0, 2, 3, 4, 5, 6, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
    30, 30, 30, 30, 30, 30,
];

use std::path::PathBuf;
use std::process::exit;
use clap::{AppSettings, Parser};
use rust_htslib::{bam, bam::Read};

fn main() {
    // initial file checks
    let args = Args::parse();
    file_exist_check(&args.i);
    let is_cram = bam_format_check(&args);
    file_absence_check(&args.o, args.n);
    
    // thread number check
    if args.p == 0 {
        eprintln!("Number of threads must be 1 or more.");
        exit(1);
    }
    if args.p > num_cpus::get() {
        eprintln!(
            "Number of threads ({}) exceeds the number of threads ({}) in your machine.",
            args.p,
            num_cpus::get()
        );
        exit(1);
    } else {
        println!(
            "Number of threads ({}) was set; {} threads found in your machine.",
            args.p,
            num_cpus::get()
        );
    }
    
    // set thread numbers
    let mut nthread_in = 0;
    let mut nthread_out = 0;
    if args.n >= 2 {
        nthread_in = args.n / 2;
        nthread_out = if (args.n % 2) == 0 {
            args.n / 2
        } else {
            1 + (args.n / 2)
        };
    }
    
    // make output file
    let mut infile = bam::Reader::from_path(&args.i).unwrap();
    if is_cram {
        infile.set_reference(&args.r).unwrap();
    }
    let header = bam::Header::from_template(infile.header());
    let mut outfile = bam::Writer::from_path(&args.o, &header, bam::Format::Cram).unwrap();
    outfile.set_reference(&args.r).unwrap();
    if nthread_in > 0 {
        infile.set_threads(nthread_in as usize).unwrap();
        outfile.set_threads(nthread_out as usize).unwrap();
    }
    
    // apply binning
    let convarr = if args.e >= 1 {
        println!("Will apply Illumina 8-binning instead of 4-binning.");
        TO_8_BINS
    } else {
        TO_4_BINS
    };
    for record in infile.records().map(|x| x.expect("Failure parsing input file")) {
        let data =
            unsafe {
                std::slice::from_raw_parts_mut(record.inner.data, record.inner().l_data as usize)
            };
        let offset = qname_capacity(&record) + record.cigar_len() * 4 + (record.seq_len() + 1) / 2;
        for i in offset..(offset + record.seq_len()) {
            data[i] = convarr[data[i] as usize];
        }
        outfile.write(&record).unwrap();
    }
}

#[derive(Parser, Debug)]
#[clap(author = "Author: Shohei Kojima @ RIKEN", version = VERSION, about = "Converts a BAM/CRAM file to a lossy compressed CRAM file with 8- or 4-binning for quality scores.", setting = AppSettings::DeriveDisplayOrder)]
#[clap(propagate_version = true)]
struct Args {
    /// Specify an input file [Required]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub i: String,
    
    /// Specify an output file [Required]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub o: String,
    
    /// Specify a refernece genome file [Required]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub r: String,
    
    /// Apply Illumina 8-binning instead of 4-binning [Optional]
    #[clap(short, long, action = clap::ArgAction::Count)]
    pub e: u8,
    
    /// Number of threads [Optional]
    #[clap(short, long, value_parser, value_name = "NUMBER", default_value_t = 1)]
    pub p: usize,
    
    /// Specify when you don't want to overwrite [Optional]
    #[clap(short, long, action = clap::ArgAction::Count)]
    pub n: u8,
}

/// judge whether the input file is BAM or CRAM
fn bam_format_check(args: &Args) -> bool {
    let mut is_cram: bool = false;
    match PathBuf::from(&args.i)
        .extension()
        .unwrap()
        .to_str()
        .unwrap()
    {
        "cram" => {
            is_cram = true;
            println!("Input file ({}) is considered as CRAM format.", args.i);
        }
        "bam" => println!("Input file ({}) is considered as BAM format.", args.i),
        _ => {
            eprintln!("Input file ({}) does not have extension .bam or .cram", args.i);
            exit(1);
        }
    }
    is_cram
}

/// re-implementation of fn qname_capacity(&self) in rust-htslib-0.40.2/src/bam/record.rs
fn qname_capacity(record: &bam::record::Record) -> usize {
    record.inner().core.l_qname as usize
}

/// Check whether a file exists.
/// Exit with 1 if the file does not exist.
fn file_exist_check(file_path: &str) {
    if !my_rust::utils::path_exists(file_path) {
        eprintln!("Input file ({}) does not exist.", file_path);
        exit(1);
    } else {
        println!("Input file ({}) found.", file_path);
    }
}

/// Check whether a file does not exist.
fn file_absence_check(file_path: &str, overwrite: u8) {
    if my_rust::utils::path_exists(file_path) {
        if overwrite == 0 {
            println!("Warn: Output will be overwritten in {}.", file_path);
        } else {
            eprintln!("Output file ({}) already exists.", file_path);
            exit(1);
        }
    } else {
        println!("Output will be written in ({}).", file_path);
    }
}
