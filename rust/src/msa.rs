//! The multialignment matrix and its position frequency matrix, matching
//! `correct_seqs.create_multialignment_matrix` and the PFM built inline in
//! `get_best_corrections` (`isONcorrect.py:837`).
//!
//! `get_best_corrections` has a consensus and one pairwise edlib alignment per
//! supporting segment. Those alignments are pairwise-to-consensus, not a
//! multiple alignment: two reads can each insert two bases at the same place
//! without the insertions occupying the same columns. This module stacks them
//! into one matrix so a column means the same thing in every row.
//!
//! # The positioned vector
//!
//! Each pairwise alignment is first re-expressed against a vector of
//! `2 * len(consensus) + 1` slots: **odd slots hold the base aligned to
//! consensus position `i`** (possibly `-`), **even slots hold the insertion
//! before it** (a string, possibly multi-character, or `-` for none). That is
//! `position_query_to_alignment`, and it makes the consensus coordinate system
//! shared across rows for free — every row has exactly the same slot count.
//!
//! Widening then happens per slot: a slot whose longest insertion is one
//! character stays one column wide, and a longer one is padded out to
//! `"-" + longest + "-"` and every other insertion is fitted into it by
//! [`best_solution`].
//!
//! # What decides the answer
//!
//! * **the widest insertion is chosen by `sorted(...)[0]` among the longest** —
//!   lexicographically smallest, not first-seen, so a `BTreeSet` is the natural
//!   container and iteration order does not leak;
//! * `position_solutions` is keyed by the padded insertion **string**, so the
//!   cache is shared across slots that happen to widen to the same thing. Pure
//!   function of `(max_insertion, insertion)`, so this is a cache and nothing
//!   more — but it means the same insertion always resolves identically;
//! * [`best_solution`] falls through four strategies in order, and the third
//!   calls edlib with `task="path"` again ([`crate::align::global`]), so its
//!   tie-break reaches the matrix too;
//! * rows keep the partition's insertion order. The PFM does not care, but
//!   later stages iterate the matrix and the corrected output does.

use rustc_hash::FxHashMap;

use crate::align;

/// The alphabet the reference's PFM counts, in its literal dict order:
/// `{"A": 0, "C": 0, "G": 0, "T": 0, "U": 0, "-": 0, "N": 0}`.
///
/// Note this is the seven-symbol dict from `get_best_corrections`, which
/// includes `N`. `correct_seqs` has a six-symbol variant without it; that one
/// belongs to `correct_to_consensus`, which isONcorrect never calls.
pub const ALPHABET: [u8; 7] = [b'A', b'C', b'G', b'T', b'U', b'-', b'N'];

/// Column of the PFM: one count per [`ALPHABET`] symbol.
///
/// A `HashMap` per column would match the reference's shape, but the alphabet
/// is seven fixed symbols and there is one column per alignment position, so an
/// array is the same data in a fraction of the space. See *Working agreements*
/// in PORTING.md.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Column([u32; ALPHABET.len()]);

impl Column {
    /// Count for `nucl`, or 0 for a symbol outside the alphabet.
    ///
    /// The reference would raise `KeyError` instead; nothing reaches it with a
    /// symbol the matrix does not contain, and the matrix can only contain what
    /// [`multialignment_matrix`] admits.
    pub fn get(&self, nucl: u8) -> u32 {
        symbol_index(nucl).map_or(0, |i| self.0[i])
    }

    fn add(&mut self, nucl: u8, n: u32) {
        if let Some(i) = symbol_index(nucl) {
            self.0[i] += n;
        }
    }
}

fn symbol_index(nucl: u8) -> Option<usize> {
    ALPHABET.iter().position(|&c| c == nucl)
}

/// Position frequency matrix over `rows`, one column per alignment position.
///
/// Every row counts once — the reference does `PFM[j][n] += 1`, with no
/// weighting by depth. **The consensus row is excluded by the caller**:
/// `get_best_corrections` skips the `"ref"` key, and since
/// [`multialignment_matrix`] returns rows positionally, the caller passes the
/// supporting rows only.
pub fn pfm(rows: &[Vec<u8>]) -> Vec<Column> {
    let width = rows.first().map_or(0, Vec::len);
    let mut pfm = vec![Column::default(); width];
    for row in rows {
        for (j, &n) in row.iter().enumerate() {
            pfm[j].add(n, 1);
        }
    }
    pfm
}

/// A positioned vector: `2 * t + 1` slots laid out in one flat buffer.
///
/// The obvious representation is `Vec<Vec<u8>>`, one allocation per slot, which
/// is the shape the reference's list of strings has. That is
/// `2 * len(consensus) + 1` tiny allocations per supporting segment per
/// correction interval, and nearly all of them hold the single byte `"-"`. On a
/// profile of ten real clusters `create_multialignment_matrix` was **14.3%** of
/// runtime, most of it allocator traffic.
///
/// Every slot is either one character or a **contiguous run of the aligned
/// query** — an insertion is exactly the query characters opposite a run of
/// consensus gaps — so all of them fit in one buffer with an offset table. Two
/// allocations per segment instead of `2t + 1`.
///
/// This is a representation change only. The slot *contents* are identical, so
/// the 192 252-row matrix oracle must not move.
#[derive(Debug, Clone, Default)]
pub struct Positioned {
    buf: Vec<u8>,
    slots: Vec<(u32, u32)>,
}

impl Positioned {
    /// Number of slots, always `2 * t + 1`.
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// Contents of slot `i`.
    #[inline]
    pub fn slot(&self, i: usize) -> &[u8] {
        let (off, len) = self.slots[i];
        &self.buf[off as usize..off as usize + len as usize]
    }

    /// The slots in order.
    pub fn iter(&self) -> impl Iterator<Item = &[u8]> + '_ {
        self.slots
            .iter()
            .map(move |&(off, len)| &self.buf[off as usize..off as usize + len as usize])
    }

    fn push(&mut self, bytes: &[u8]) {
        let off = self.buf.len() as u32;
        self.buf.extend_from_slice(bytes);
        self.slots.push((off, bytes.len() as u32));
    }
}

/// Re-express a pairwise alignment against the consensus slot vector, matching
/// `position_query_to_alignment(query_aligned, target_aligned, 0)`.
///
/// Returns `2 * t + 1` slots, where `t` is the number of non-gap characters in
/// `target_aligned`. Even slots are insertions (`-` when there is none), odd
/// slots are the aligned query character.
///
/// The reference also returns the start and end vector positions, but they are
/// `0` and `2 * t` by construction and it immediately asserts as much.
pub fn positioned(query_aligned: &[u8], target_aligned: &[u8]) -> Positioned {
    assert_eq!(
        query_aligned.len(),
        target_aligned.len(),
        "the two rows of a pairwise alignment must have equal length"
    );
    let n = target_aligned.len();
    let mut out = Positioned {
        buf: Vec::with_capacity(n + n / 2 + 2),
        slots: Vec::with_capacity(n * 2 + 1),
    };

    // An insertion is a contiguous run of query characters opposite consensus
    // gaps, so it is tracked as a range into `query_aligned` and copied once.
    let mut run_start: Option<usize> = None;
    for (p, &t) in target_aligned.iter().enumerate() {
        if t == b'-' {
            run_start.get_or_insert(p);
        } else {
            match run_start.take() {
                Some(start) => out.push(&query_aligned[start..p]),
                None => out.push(b"-"),
            }
            out.push(&query_aligned[p..p + 1]);
        }
    }
    match run_start {
        Some(start) => out.push(&query_aligned[start..]),
        None => out.push(b"-"),
    }
    out
}

/// Fit `q_ins` into the padded slot `max_insertion`, matching
/// `correct_seqs.get_best_solution`.
///
/// Four strategies, tried in order, and the order is the behaviour:
///
/// 1. **no insertion** (`q_ins == "-"`) — the slot is all gaps;
/// 2. **substring** — if `q_ins` occurs in `max_insertion`, place it there.
///    `str.find` returns the *leftmost* occurrence;
/// 3. **alignment** — align `max_insertion` to `q_ins` and thread `q_ins`
///    through, gapping where `max_insertion` has extra characters. Rejected if
///    the alignment deletes from `q_ins` (`"D" in cigar`), since that would
///    drop bases the read actually has;
/// 4. **shift** — the offset with the most matching characters, and failing
///    that, flush left.
///
/// The returned slot always has `max_insertion.len()` characters.
pub fn best_solution(max_insertion: &[u8], q_ins: &[u8]) -> Vec<u8> {
    if q_ins == b"-" {
        return vec![b'-'; max_insertion.len()];
    }

    // 2. Leftmost occurrence.
    if let Some(pos) = find_subslice(max_insertion, q_ins) {
        let mut out = vec![b'-'; pos];
        out.extend_from_slice(q_ins);
        out.resize(max_insertion.len(), b'-');
        return out;
    }

    // 3. Thread through an alignment, if it does not delete from `q_ins`.
    if let Some(threaded) = min_ed(max_insertion, q_ins) {
        return threaded;
    }

    // 4. Best-matching offset, then flush left.
    //
    // `max_p` starts at 0 and the comparison is strict, so an offset of 0 is
    // never *selected* here — it falls through to the flush-left branch, which
    // produces the same string anyway.
    let (mut max_p, mut max_matches) = (0usize, 0usize);
    for p in 0..=max_insertion.len().saturating_sub(q_ins.len()) {
        let matches = q_ins
            .iter()
            .zip(&max_insertion[p..])
            .filter(|(a, b)| a == b)
            .count();
        if matches > max_matches {
            max_p = p;
            max_matches = matches;
        }
    }

    let mut out = vec![b'-'; max_p];
    out.extend_from_slice(&q_ins[..q_ins.len().min(max_insertion.len() - max_p)]);
    out.resize(max_insertion.len(), b'-');
    out
}

/// Thread `q_ins` through an alignment to `max_insertion`, matching
/// `correct_seqs.min_ed`.
///
/// The alignment runs `max_insertion` as the *query*, so `I` consumes
/// `max_insertion` (a position `q_ins` does not reach, hence a gap) and `D`
/// would consume `q_ins`. Any `D` means the alignment drops a character of
/// `q_ins`, which is not allowed — the reference returns `""`, falsy, and the
/// caller moves on.
fn min_ed(max_insertion: &[u8], q_ins: &[u8]) -> Option<Vec<u8>> {
    let cigar = align::global(max_insertion, q_ins).cigar;
    if cigar.contains('D') {
        return None;
    }
    let ops = align::parse_cigar(&cigar)?;

    let mut out = Vec::with_capacity(max_insertion.len());
    let mut q_pos = 0usize;
    for (len, op) in ops {
        if op == b'I' {
            out.extend(std::iter::repeat_n(b'-', len));
        } else {
            out.extend_from_slice(&q_ins[q_pos..q_pos + len]);
            q_pos += len;
        }
    }
    Some(out)
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() > haystack.len() {
        return None;
    }
    (0..=haystack.len() - needle.len()).find(|&i| &haystack[i..i + needle.len()] == needle)
}

/// The multialignment matrix, matching `create_multialignment_matrix` composed
/// with `create_multialignment_format_NEW(dict, 0, 2 * len(consensus))`.
///
/// `rows` are `(query_aligned, target_aligned)` pairs in partition order — the
/// consensus row included, as `get_best_corrections` puts it there under the
/// `"ref"` key. Every `target_aligned` must expand to the same consensus, so
/// every positioned vector has the same length; the reference asserts this.
///
/// The reference's window arguments are always the full vector, so the
/// "reads covering the region" filter admits every row and is not reproduced.
///
/// # Panics
///
/// On an empty `rows`, on rows that disagree about the consensus length, or on
/// a character outside `ACGTU-N` — all three are `assert`s or a `KeyError` in
/// the reference, and none is reachable from a well-formed alignment.
pub fn multialignment_matrix(rows: &[(&[u8], &[u8])]) -> Vec<Vec<u8>> {
    assert!(!rows.is_empty(), "multialignment of an empty partition");

    let segments: Vec<Positioned> = rows.iter().map(|&(q, t)| positioned(q, t)).collect();

    let nr_pos = segments[0].len();
    assert!(
        segments.iter().all(|s| s.len() == nr_pos),
        "rows disagree about the consensus length"
    );

    // How wide each slot becomes: `-` for a one-character slot, else the
    // lexicographically smallest of the longest insertions, padded either side.
    //
    // The reference takes `sorted(set(...))[0]` of the longest, and this used to
    // mirror that with a `BTreeSet` per slot — `2t + 1` trees per matrix, each
    // paying a string compare per insert per segment. But the *only* things the
    // set was ever asked for were the maximum length and the lexicographically
    // smallest string at that length, and a single pass tracks both directly.
    let mut max_insertions: Vec<Vec<u8>> = Vec::with_capacity(nr_pos);
    for p in 0..nr_pos {
        let mut best: Option<&[u8]> = None;
        for segment in &segments {
            let slot = segment.slot(p);
            // Longer wins outright; equal length breaks lexicographically, which
            // is exactly `sorted(...)[0]` over the longest.
            let better = match best {
                None => true,
                // Longer wins outright; at equal length the *smaller* string
                // wins, because the reference takes `sorted(...)[0]`.
                Some(b) => slot.len() > b.len() || (slot.len() == b.len() && slot < b),
            };
            if better {
                best = Some(slot);
            }
        }
        let widest = best.expect("a non-empty partition has slots");
        if widest.len() <= 1 {
            max_insertions.push(vec![b'-']);
            continue;
        }
        assert_eq!(p % 2, 0, "an insertion in a consensus-base slot");
        let mut padded = Vec::with_capacity(widest.len() + 2);
        padded.push(b'-');
        padded.extend_from_slice(widest);
        padded.push(b'-');
        max_insertions.push(padded);
    }

    // Solutions, keyed by `(padded slot, insertion)` — the reference's
    // `position_solutions`, shared across slots exactly as there, and filled on
    // demand rather than pre-enumerated. It is a cache, so *when* an entry is
    // computed cannot reach the output.
    let mut solutions: FxHashMap<(&[u8], &[u8]), Vec<u8>> = FxHashMap::default();

    segments
        .iter()
        .map(|segment| {
            let mut row = Vec::new();
            for (p, slot) in segment.iter().enumerate() {
                let max_ins = &max_insertions[p];
                if max_ins.len() <= 1 {
                    // A one-character slot is the identity: `position_solutions`
                    // maps each of `ACGTU-N` to itself.
                    assert!(
                        slot.len() == 1 && symbol_index(slot[0]).is_some(),
                        "character outside ACGTU-N in the alignment"
                    );
                    row.push(slot[0]);
                } else {
                    let key = (max_ins.as_slice(), slot);
                    let sol = solutions
                        .entry(key)
                        .or_insert_with(|| best_solution(max_ins, slot));
                    row.extend_from_slice(sol);
                }
            }
            row
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[u8]) -> String {
        String::from_utf8(v.to_vec()).unwrap()
    }

    fn slots(query: &str, target: &str) -> Vec<String> {
        positioned(query.as_bytes(), target.as_bytes())
            .iter()
            .map(s)
            .collect()
    }

    fn matrix(rows: &[(&str, &str)]) -> Vec<String> {
        let pairs: Vec<(&[u8], &[u8])> = rows
            .iter()
            .map(|(q, t)| (q.as_bytes(), t.as_bytes()))
            .collect();
        multialignment_matrix(&pairs).iter().map(|v| s(v)).collect()
    }

    #[test]
    fn a_gapless_alignment_interleaves_empty_insertion_slots() {
        assert_eq!(
            slots("ACGT", "ACGT"),
            ["-", "A", "-", "C", "-", "G", "-", "T", "-"]
        );
    }

    #[test]
    fn an_insertion_lands_in_the_slot_before_its_consensus_base() {
        // Query has TT between consensus C and G.
        assert_eq!(
            slots("ACTTGT", "AC--GT"),
            ["-", "A", "-", "C", "TT", "G", "-", "T", "-"]
        );
    }

    #[test]
    fn a_trailing_insertion_lands_in_the_final_slot() {
        assert_eq!(
            slots("ACGTAA", "ACGT--"),
            ["-", "A", "-", "C", "-", "G", "-", "T", "AA"]
        );
    }

    #[test]
    fn a_deletion_shows_as_a_gap_in_its_base_slot() {
        assert_eq!(
            slots("A-GT", "ACGT"),
            ["-", "A", "-", "-", "-", "G", "-", "T", "-"]
        );
    }

    #[test]
    fn the_slot_count_is_two_per_consensus_base_plus_one() {
        for (q, t) in [("ACGT", "ACGT"), ("ACTTGT", "AC--GT"), ("A-GT", "ACGT")] {
            let bases = t.bytes().filter(|&c| c != b'-').count();
            assert_eq!(positioned(q.as_bytes(), t.as_bytes()).len(), 2 * bases + 1);
        }
    }

    #[test]
    fn rows_without_insertions_widen_to_one_column_per_slot() {
        let m = matrix(&[("ACGT", "ACGT"), ("AGGT", "ACGT")]);
        assert_eq!(m, ["-A-C-G-T-", "-A-G-G-T-"]);
    }

    #[test]
    fn a_longer_insertion_widens_the_slot_for_every_row() {
        // Row 0 inserts TT, row 1 inserts nothing, row 2 inserts T.
        let m = matrix(&[("ACTTGT", "AC--GT"), ("ACGT", "ACGT"), ("ACTGT", "AC-GT")]);
        // The slot becomes "-TT-": four columns in every row.
        assert_eq!(m[0], "-A-C-TT-G-T-");
        assert_eq!(m[1], "-A-C----G-T-");
        assert_eq!(m[2], "-A-C-T--G-T-");
        assert!(m.iter().all(|r| r.len() == m[0].len()));
    }

    #[test]
    fn the_widest_insertion_is_the_lexicographically_smallest_of_the_longest() {
        // Two insertions of length 2: "TT" and "GG". "GG" sorts first.
        let m = matrix(&[("ACTTGT", "AC--GT"), ("ACGGGT", "AC--GT")]);
        assert_eq!(
            m[1], "-A-C-GG-G-T-",
            "the widest slot is GG, matched exactly"
        );
        // TT is not a substring of "-GG-", so it is threaded through instead.
        assert_eq!(m[0].replace('-', ""), "ACTTGT");
    }

    #[test]
    fn no_insertion_fills_the_widened_slot_with_gaps() {
        assert_eq!(s(&best_solution(b"-TT-", b"-")), "----");
    }

    #[test]
    fn a_substring_insertion_is_placed_at_its_leftmost_occurrence() {
        assert_eq!(s(&best_solution(b"-ACGT-", b"CG")), "--CG--");
        assert_eq!(s(&best_solution(b"-AA-", b"A")), "-A--");
    }

    #[test]
    fn a_non_substring_insertion_is_threaded_through_an_alignment() {
        // AG against GACG: no substring match, but A and G can both be placed.
        let got = s(&best_solution(b"-GACG-", b"AG"));
        assert_eq!(got.len(), 6);
        assert_eq!(got.replace('-', ""), "AG");
    }

    #[test]
    fn a_solution_always_fills_the_slot_exactly() {
        for max_ins in [&b"-A-"[..], b"-AC-", b"-ACGT-", b"-TTTT-"] {
            for q in [&b"-"[..], b"A", b"C", b"AC", b"CA", b"GGG", b"ACGT"] {
                let sol = best_solution(max_ins, q);
                assert_eq!(sol.len(), max_ins.len(), "{:?} into {:?}", s(q), s(max_ins));
            }
        }
    }

    #[test]
    fn the_pfm_counts_every_row_once_over_the_seven_symbol_alphabet() {
        let rows: Vec<Vec<u8>> = ["A-CN", "AGCU", "A-CN"]
            .iter()
            .map(|r| r.as_bytes().to_vec())
            .collect();
        let p = pfm(&rows);
        assert_eq!(p.len(), 4);
        assert_eq!(p[0].get(b'A'), 3);
        assert_eq!(p[1].get(b'-'), 2);
        assert_eq!(p[1].get(b'G'), 1);
        assert_eq!(p[3].get(b'N'), 2);
        assert_eq!(p[3].get(b'U'), 1);
        // Every column sums to the row count.
        for col in &p {
            assert_eq!(ALPHABET.iter().map(|&c| col.get(c)).sum::<u32>(), 3);
        }
    }

    #[test]
    fn the_pfm_of_no_rows_has_no_columns() {
        assert!(pfm(&[]).is_empty());
    }
}

/// Differential oracle against `correct_seqs.create_multialignment_matrix`.
///
/// The matrix is where several non-unique choices land at once — edlib's CIGAR
/// tie-break, the widest-insertion pick, and [`best_solution`]'s four-way
/// fallthrough — so real cases are worth far more than constructed ones.
///
/// ```bash
/// bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/msa
/// MSA_CASES=/tmp/msa/msa_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml msa::oracle -- --nocapture
/// ```
#[cfg(test)]
mod oracle {
    use super::*;

    /// `#case\trow\tquery_aligned\ttarget_aligned\taln_row`, rows of a case
    /// contiguous and in partition order.
    #[test]
    fn matches_the_reference_matrix_exactly() {
        let Ok(path) = std::env::var("MSA_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("MSA_CASES={path} unreadable: {e}"));

        let mut cases: Vec<Vec<(String, String, String)>> = Vec::new();
        let mut current = String::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 5, "short row: {line}");
            if f[0] != current {
                current = f[0].to_string();
                cases.push(Vec::new());
            }
            cases
                .last_mut()
                .unwrap()
                .push((f[2].to_string(), f[3].to_string(), f[4].to_string()));
        }
        assert!(!cases.is_empty(), "no cases in {path}");

        let (mut rows, mut bad, mut shown) = (0usize, 0usize, 0usize);
        for (c, case) in cases.iter().enumerate() {
            let pairs: Vec<(&[u8], &[u8])> = case
                .iter()
                .map(|(q, t, _)| (q.as_bytes(), t.as_bytes()))
                .collect();
            let got = multialignment_matrix(&pairs);
            for (r, (want, got)) in case.iter().map(|c| &c.2).zip(&got).enumerate() {
                rows += 1;
                if want.as_bytes() != got.as_slice() {
                    if shown < 3 {
                        shown += 1;
                        eprintln!(
                            "case {c} row {r}:\n  python: {want}\n  rust  : {}",
                            String::from_utf8_lossy(got)
                        );
                    }
                    bad += 1;
                }
            }
        }
        assert_eq!(bad, 0, "{bad} of {rows} matrix rows differed");
        eprintln!(
            "matched the reference on all {rows} rows across {} matrices",
            cases.len()
        );
    }
}
