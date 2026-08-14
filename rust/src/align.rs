//! Global alignment with an edlib-compatible CIGAR, matching
//! `edlib.align(query, target, task="path", mode="NW")`.
//!
//! `get_best_corrections` aligns every supporting segment back to the consensus
//! and consumes the **CIGAR**, not just the distance. That makes this the one
//! place where matching edlib's *choice* of alignment matters: many alignments
//! achieve the same optimal edit distance, and a different tie-break yields a
//! different-but-equally-optimal CIGAR — which changes the MSA and therefore the
//! corrected output.
//!
//! # The tie-break
//!
//! Determined empirically rather than guessed. Running the traceback with all
//! six preference orderings against 677 real `task="path"` calls recorded from
//! the reference:
//!
//! | preference (backwards from the end) | CIGARs matching |
//! | --- | --- |
//! | insert, delete, diagonal | **677 / 677** |
//! | delete, insert, diagonal | 637 / 677 |
//! | delete, diagonal, insert | 281 / 677 |
//! | insert, diagonal, delete | 131 / 677 |
//! | diagonal, insert, delete | 56 / 677 |
//! | diagonal, delete, insert | 55 / 677 |
//!
//! So edlib prefers **insertion, then deletion, then the diagonal**, walking
//! backwards from `(n, m)`. Note how badly the intuitive diagonal-first ordering
//! does — 8% — which is why this was measured.
//!
//! Every ordering reproduces the edit *distance*; only the path differs. That is
//! the whole point: the distance is unique, the alignment is not.
//!
//! CIGAR operations use edlib's extended alphabet: `=` match, `X` mismatch,
//! `I` insertion to target (consumes query), `D` deletion from target.

/// Alignment of `query` against `target`: edit distance and extended CIGAR.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub edit_distance: usize,
    pub cigar: String,
}

/// Globally align `query` to `target`, returning an edlib-compatible CIGAR.
///
/// Straightforward `O(n*m)` Needleman-Wunsch with unit costs. The segments this
/// runs on are short — anchor spans are bounded by `--xmax`, 80 by default — so
/// the quadratic table is small. edlib is bit-parallel and would be faster on
/// long inputs; if that ever matters, the oracle test below is what makes a
/// faster implementation safe to swap in.
pub fn global(query: &[u8], target: &[u8]) -> Alignment {
    let (n, m) = (query.len(), target.len());

    let mut dp = vec![0u32; (n + 1) * (m + 1)];
    let at = |i: usize, j: usize| i * (m + 1) + j;
    for i in 0..=n {
        dp[at(i, 0)] = i as u32;
    }
    for j in 0..=m {
        dp[at(0, j)] = j as u32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let sub = dp[at(i - 1, j - 1)] + u32::from(query[i - 1] != target[j - 1]);
            dp[at(i, j)] = sub.min(dp[at(i - 1, j)] + 1).min(dp[at(i, j - 1)] + 1);
        }
    }

    // Traceback, preferring insertion, then deletion, then the diagonal.
    let (mut i, mut j) = (n, m);
    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    while i > 0 || j > 0 {
        if i > 0 && dp[at(i - 1, j)] + 1 == dp[at(i, j)] {
            ops.push(b'I');
            i -= 1;
        } else if j > 0 && dp[at(i, j - 1)] + 1 == dp[at(i, j)] {
            ops.push(b'D');
            j -= 1;
        } else if i > 0 && j > 0 {
            let same = query[i - 1] == target[j - 1];
            ops.push(if same { b'=' } else { b'X' });
            i -= 1;
            j -= 1;
        } else {
            break;
        }
    }
    ops.reverse();

    Alignment {
        edit_distance: dp[at(n, m)] as usize,
        cigar: run_length_encode(&ops),
    }
}

/// One run of a CIGAR: `(length, operation)`.
pub type CigarOp = (usize, u8);

/// Parse an extended CIGAR into runs, matching `help_functions.cigar_to_seq`'s
/// `re.split(r'[=DXSMI]+', cigar)` walk.
///
/// The reference's regex tolerates `S` and `M` when splitting but its loop then
/// falls through to `print("error"); sys.exit()` for anything that is not
/// `=`, `X`, `I` or `D`. Nothing in the pipeline produces those operations —
/// edlib emits only the extended four — so this returns `None` instead of
/// aborting the process, and callers treat it as a bug rather than an input.
pub fn parse_cigar(cigar: &str) -> Option<Vec<CigarOp>> {
    let mut ops = Vec::new();
    let mut len = 0usize;
    let mut seen_digit = false;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            len = len * 10 + usize::from(c - b'0');
            seen_digit = true;
        } else {
            if !seen_digit || !matches!(c, b'=' | b'X' | b'I' | b'D') {
                return None;
            }
            ops.push((len, c));
            len = 0;
            seen_digit = false;
        }
    }
    if seen_digit {
        return None; // trailing length with no operation
    }
    Some(ops)
}

/// Expand a CIGAR into the two gapped strings, matching
/// `help_functions.cigar_to_seq(cigar, query, ref)`.
///
/// Returns `(query_aligned, ref_aligned)` — the reference's return order, which
/// `get_best_corrections` immediately unpacks as `read_alignment,
/// ref_alignment`. `I` is an insertion with respect to the reference (a gap in
/// `ref_aligned`); `D` is a deletion (a gap in `query_aligned`).
///
/// Returns `None` if the CIGAR does not parse or does not cover both sequences.
pub fn cigar_to_seq(cigar: &str, query: &[u8], reference: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    let ops = parse_cigar(cigar)?;
    let (mut q_index, mut r_index) = (0usize, 0usize);
    let (mut q_aln, mut r_aln) = (Vec::new(), Vec::new());

    for (len, op) in ops {
        match op {
            b'=' | b'X' => {
                q_aln.extend_from_slice(query.get(q_index..q_index + len)?);
                r_aln.extend_from_slice(reference.get(r_index..r_index + len)?);
                q_index += len;
                r_index += len;
            }
            b'I' => {
                r_aln.extend(std::iter::repeat_n(b'-', len));
                q_aln.extend_from_slice(query.get(q_index..q_index + len)?);
                q_index += len;
            }
            b'D' => {
                r_aln.extend_from_slice(reference.get(r_index..r_index + len)?);
                q_aln.extend(std::iter::repeat_n(b'-', len));
                r_index += len;
            }
            _ => return None,
        }
    }
    // Python would silently accept a CIGAR that stops short; nothing produces
    // one, and a partial alignment would corrupt the MSA rather than fail.
    if q_index != query.len() || r_index != reference.len() {
        return None;
    }
    Some((q_aln, r_aln))
}

/// `[b'=', b'=', b'X']` becomes `"2=1X"`.
fn run_length_encode(ops: &[u8]) -> String {
    let mut out = String::new();
    let mut k = 0;
    while k < ops.len() {
        let c = ops[k];
        let start = k;
        while k < ops.len() && ops[k] == c {
            k += 1;
        }
        out.push_str(&(k - start).to_string());
        out.push(c as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(q: &str, t: &str) -> String {
        global(q.as_bytes(), t.as_bytes()).cigar
    }

    #[test]
    fn identical_sequences_are_all_matches() {
        let a = global(b"ACGT", b"ACGT");
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.cigar, "4=");
    }

    #[test]
    fn a_substitution_shows_as_x() {
        let a = global(b"ACGT", b"AGGT");
        assert_eq!(a.edit_distance, 1);
        assert_eq!(a.cigar, "1=1X2=");
    }

    #[test]
    fn empty_query_is_all_deletions() {
        let a = global(b"", b"ACGT");
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.cigar, "4D");
    }

    #[test]
    fn empty_target_is_all_insertions() {
        let a = global(b"ACGT", b"");
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.cigar, "4I");
    }

    #[test]
    fn both_empty_yields_an_empty_cigar() {
        let a = global(b"", b"");
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.cigar, "");
    }

    #[test]
    fn cigar_operation_counts_cover_both_sequences() {
        // Query length = = + X + I; target length = = + X + D.
        let (q, t) = ("ACGTACGTAA", "ACGTTACGT");
        let c = cigar(q, t);
        let (mut qn, mut tn, mut num) = (0usize, 0usize, String::new());
        for ch in c.chars() {
            if ch.is_ascii_digit() {
                num.push(ch);
            } else {
                let n: usize = num.parse().unwrap();
                num.clear();
                match ch {
                    '=' | 'X' => {
                        qn += n;
                        tn += n;
                    }
                    'I' => qn += n,
                    'D' => tn += n,
                    _ => panic!("unexpected op {ch}"),
                }
            }
        }
        assert_eq!(qn, q.len(), "query coverage");
        assert_eq!(tn, t.len(), "target coverage");
    }

    fn expand(cigar: &str, q: &str, r: &str) -> (String, String) {
        let (qa, ra) = cigar_to_seq(cigar, q.as_bytes(), r.as_bytes()).unwrap();
        (
            String::from_utf8(qa).unwrap(),
            String::from_utf8(ra).unwrap(),
        )
    }

    #[test]
    fn cigar_parses_into_runs() {
        assert_eq!(
            parse_cigar("2=1X10I").unwrap(),
            vec![(2, b'='), (1, b'X'), (10, b'I')]
        );
        assert_eq!(parse_cigar("").unwrap(), vec![]);
        assert!(parse_cigar("3M").is_none(), "M is not in the pipeline");
        assert!(parse_cigar("3").is_none(), "length with no operation");
        assert!(parse_cigar("=").is_none(), "operation with no length");
    }

    #[test]
    fn insertion_gaps_the_reference_and_deletion_gaps_the_query() {
        let (q, r) = expand("2=2I2=", "ACGTAC", "ACAC");
        assert_eq!((q.as_str(), r.as_str()), ("ACGTAC", "AC--AC"));
        let (q, r) = expand("2=2D2=", "ACAC", "ACGTAC");
        assert_eq!((q.as_str(), r.as_str()), ("AC--AC", "ACGTAC"));
    }

    #[test]
    fn an_expanded_cigar_round_trips_its_own_alignment() {
        let (q, r) = ("ACGTACGTAA", "ACGTTACGT");
        let a = global(q.as_bytes(), r.as_bytes());
        let (qa, ra) = expand(&a.cigar, q, r);
        assert_eq!(qa.len(), ra.len());
        assert_eq!(qa.replace('-', ""), q);
        assert_eq!(ra.replace('-', ""), r);
        // Gapped columns account for exactly the edit distance.
        let diffs = qa.bytes().zip(ra.bytes()).filter(|(a, b)| a != b).count();
        assert_eq!(diffs, a.edit_distance);
    }

    #[test]
    fn a_cigar_that_does_not_cover_both_sequences_is_rejected() {
        assert!(cigar_to_seq("2=", b"ACGT", b"ACGT").is_none());
        assert!(cigar_to_seq("4=", b"ACGT", b"ACG").is_none());
    }

    /// Tie cases verified directly against Python edlib, not assumed.
    ///
    /// The preference is applied walking *backwards*, so the gap lands at the
    /// end of the emitted CIGAR rather than the start --- which is the opposite
    /// of what "insertion first" suggests when read forwards.
    #[test]
    fn tie_cases_match_edlib_ground_truth() {
        assert_eq!(cigar("AA", "A"), "1=1I");
        assert_eq!(cigar("A", "AA"), "1=1D");
        assert_eq!(cigar("ACGT", "ACT"), "2=1I1=");
        assert_eq!(cigar("AAA", "A"), "1=2I");
    }
}

/// Differential oracle against Python's edlib.
///
/// ```bash
/// bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/cg
/// CIGAR_CASES=/tmp/cg/cigar_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml align::oracle -- --nocapture
/// ```
#[cfg(test)]
mod oracle {
    use super::*;

    #[test]
    fn matches_edlib_cigars_exactly() {
        let Ok(path) = std::env::var("CIGAR_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("CIGAR_CASES={path} unreadable: {e}"));

        let (mut n, mut bad_cigar, mut bad_ed, mut shown) = (0usize, 0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 4 {
                continue;
            }
            let got = global(f[0].as_bytes(), f[1].as_bytes());
            let want_ed: usize = f[2].parse().unwrap();
            if got.edit_distance != want_ed {
                bad_ed += 1;
            }
            if got.cigar != f[3] {
                if shown < 3 {
                    shown += 1;
                    eprintln!("cigar mismatch:\n  edlib: {}\n  rust : {}", f[3], got.cigar);
                }
                bad_cigar += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(bad_ed, 0, "{bad_ed} of {n} edit distances differed");
        assert_eq!(bad_cigar, 0, "{bad_cigar} of {n} CIGARs differed");
        eprintln!("matched edlib on all {n} recorded path alignments");
    }

    /// Replay `help_functions.cigar_to_seq`: `#cigar query ref q_aln r_aln`.
    ///
    /// Separate from the CIGAR oracle above because the two can fail
    /// independently — expansion is deterministic given a CIGAR, so a mismatch
    /// here is this function's bug and not the aligner's.
    ///
    /// ```bash
    /// CIGAR_TO_SEQ_CASES=/tmp/cg/cigar_to_seq_cases.tsv \
    ///   cargo test --manifest-path rust/Cargo.toml align::oracle -- --nocapture
    /// ```
    #[test]
    fn expands_cigars_like_the_reference() {
        let Ok(path) = std::env::var("CIGAR_TO_SEQ_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("CIGAR_TO_SEQ_CASES={path} unreadable: {e}"));

        let (mut n, mut bad) = (0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 5, "short row: {line}");
            let (q_aln, r_aln) = cigar_to_seq(f[0], f[1].as_bytes(), f[2].as_bytes())
                .unwrap_or_else(|| panic!("rejected a CIGAR the reference expanded: {}", f[0]));
            if q_aln != f[3].as_bytes() || r_aln != f[4].as_bytes() {
                if bad < 3 {
                    eprintln!(
                        "mismatch on {}:\n  python: {} / {}\n  rust  : {} / {}",
                        f[0],
                        f[3],
                        f[4],
                        String::from_utf8_lossy(&q_aln),
                        String::from_utf8_lossy(&r_aln)
                    );
                }
                bad += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(bad, 0, "{bad} of {n} expansions differed");
        eprintln!("expanded all {n} recorded CIGARs identically");
    }
}
