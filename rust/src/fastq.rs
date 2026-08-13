//! Fastq reading, matching `help_functions.readfq`.
//!
//! The reference's reader is a hand-rolled state machine with a few habits that
//! are load-bearing and easy to miss:
//!
//! * the record name is everything after the leading `@`/`>`, with **spaces
//!   replaced by underscores** --- so `@read1 len=42` becomes `read1_len=42`;
//! * sequence and quality may span multiple lines;
//! * quality lines are consumed until their total length reaches the sequence
//!   length, not until a fixed line count;
//! * `r_id` is the 1-based position in the file, and correction is
//!   order-dependent, so the order here is part of the contract.

use std::io::{self, BufRead};

/// One fastq record, as the reference models it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    /// 1-based input order. Reads are processed in `sorted(reads)` order.
    pub r_id: usize,
    /// Header without the leading sigil, spaces replaced by underscores.
    pub acc: String,
    pub seq: String,
    pub qual: String,
}

/// Read every record, preserving input order.
pub fn read_fastq<R: BufRead>(reader: R) -> io::Result<Vec<Read>> {
    let mut lines = Vec::new();
    for line in reader.lines() {
        lines.push(line?);
    }

    let mut reads = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        // Find the next header.
        if !lines[i].starts_with('@') && !lines[i].starts_with('>') {
            i += 1;
            continue;
        }
        let acc = lines[i][1..].replace(' ', "_");
        i += 1;

        // Sequence lines run until the next line opening a new section.
        let mut seq = String::new();
        while i < lines.len() {
            let l = &lines[i];
            if l.starts_with('@') || l.starts_with('+') || l.starts_with('>') {
                break;
            }
            seq.push_str(l);
            i += 1;
        }

        // A '+' line means fastq; anything else means this was a fasta record.
        if i < lines.len() && lines[i].starts_with('+') {
            i += 1;
            let mut qual = String::new();
            while i < lines.len() && qual.len() < seq.len() {
                qual.push_str(&lines[i]);
                i += 1;
            }
            reads.push(Read {
                r_id: reads.len() + 1,
                acc,
                seq,
                qual,
            });
        } else {
            reads.push(Read {
                r_id: reads.len() + 1,
                acc,
                seq,
                qual: String::new(),
            });
        }
    }
    Ok(reads)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Vec<Read> {
        read_fastq(io::Cursor::new(s)).unwrap()
    }

    #[test]
    fn reads_a_simple_record() {
        let r = parse("@r1\nACGT\n+\n!!!!\n");
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].r_id, 1);
        assert_eq!(r[0].acc, "r1");
        assert_eq!(r[0].seq, "ACGT");
        assert_eq!(r[0].qual, "!!!!");
    }

    #[test]
    fn r_id_is_one_based_input_order() {
        let r = parse("@a\nAC\n+\n!!\n@b\nGT\n+\n!!\n@c\nTT\n+\n!!\n");
        assert_eq!(
            r.iter()
                .map(|x| (x.r_id, x.acc.as_str()))
                .collect::<Vec<_>>(),
            vec![(1, "a"), (2, "b"), (3, "c")]
        );
    }

    #[test]
    fn spaces_in_header_become_underscores() {
        // readfq does `last[1:].replace(" ", "_")`, so this is observable in the
        // output fastq headers.
        let r = parse("@read1 length=42 foo\nACGT\n+\n!!!!\n");
        assert_eq!(r[0].acc, "read1_length=42_foo");
    }

    #[test]
    fn multiline_sequence_and_quality_are_joined() {
        let r = parse("@r1\nACGT\nACGT\n+\n!!!!\n????\n");
        assert_eq!(r[0].seq, "ACGTACGT");
        assert_eq!(r[0].qual, "!!!!????");
    }

    #[test]
    fn quality_stops_once_it_covers_the_sequence() {
        // The reference reads quality lines until the accumulated length
        // reaches len(seq), then hands control back.
        let r = parse("@r1\nACGT\n+\n!!!!\n@r2\nGGGG\n+\n####\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].qual, "!!!!");
        assert_eq!(r[1].seq, "GGGG");
        assert_eq!(r[1].qual, "####");
    }

    #[test]
    fn empty_input_yields_nothing() {
        assert!(parse("").is_empty());
    }
}
