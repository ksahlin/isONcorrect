//! Stage timing, mirroring `bench/profile_python.py`'s output.
//!
//! Sampling profilers are the wrong tool here for two reasons: release-build
//! inlining smears attribution (a `sample` run of this binary was full of `???`
//! frames), and the numbers cannot be lined up against the reference's own
//! profile. This instruments the same stages the Python profiler does, reports
//! the same columns, and therefore answers the question that matters — *which
//! stage is the port slower at than Python* — rather than just where its own
//! time goes.
//!
//! **Self vs inclusive.** `get_best_corrections` contains the POA, the pairwise
//! alignments and the MSA, so inclusive time double-counts nesting. Each scope
//! subtracts the inclusive time of its children, so the `self` column sums to at
//! most the wall clock and the shortfall is honest: I/O, parsing, glue.
//!
//! Off unless `ISONCORRECT_PROFILE` is set, and instrumented at *stage*
//! granularity only. A timer pair is tens of nanoseconds, which is nothing
//! against a POA call but would distort `editdist`, which runs hundreds of
//! thousands of times per cluster.
//!
//! ```bash
//! ISONCORRECT_PROFILE=1 isONcorrect --fastq cluster.fastq --outfolder /tmp/x
//! ```

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Mutex;
use std::time::{Duration, Instant};

/// One open scope: when it drops, its inclusive time is credited to the parent's
/// child total so the parent can subtract it.
struct Frame {
    label: &'static str,
    start: Instant,
    children: Duration,
}

thread_local! {
    static STACK: RefCell<Vec<Frame>> = const { RefCell::new(Vec::new()) };
}

/// `label -> (calls, self, inclusive)`.
type Totals = HashMap<&'static str, (u64, Duration, Duration)>;
static TOTALS: Mutex<Option<Totals>> = Mutex::new(None);

fn enabled() -> bool {
    use std::sync::atomic::{AtomicU8, Ordering};
    static ON: AtomicU8 = AtomicU8::new(u8::MAX);
    let cached = ON.load(Ordering::Relaxed);
    if cached != u8::MAX {
        return cached == 1;
    }
    let on = std::env::var("ISONCORRECT_PROFILE")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false);
    ON.store(u8::from(on), Ordering::Relaxed);
    on
}

/// Time a stage until the returned guard drops. Does nothing when profiling is
/// off, so leaving these in the hot path costs one relaxed load.
pub fn scope(label: &'static str) -> Option<Guard> {
    if !enabled() {
        return None;
    }
    STACK.with(|s| {
        s.borrow_mut().push(Frame {
            label,
            start: Instant::now(),
            children: Duration::ZERO,
        })
    });
    Some(Guard)
}

/// Drop guard; see [`scope`].
pub struct Guard;

impl Drop for Guard {
    fn drop(&mut self) {
        let done = STACK.with(|s| s.borrow_mut().pop());
        let Some(frame) = done else { return };
        let inclusive = frame.start.elapsed();
        let self_time = inclusive.saturating_sub(frame.children);

        // Credit the parent so it can exclude us from its own self time.
        STACK.with(|s| {
            if let Some(parent) = s.borrow_mut().last_mut() {
                parent.children += inclusive;
            }
        });

        let mut guard = TOTALS.lock().expect("profile totals");
        let totals = guard.get_or_insert_with(HashMap::new);
        let entry = totals
            .entry(frame.label)
            .or_insert((0, Duration::ZERO, Duration::ZERO));
        entry.0 += 1;
        entry.1 += self_time;
        entry.2 += inclusive;
    }
}

/// Print the table, in the same shape `bench/profile_python.py` prints.
pub fn report(total_wall: Duration) {
    if !enabled() {
        return;
    }
    let guard = TOTALS.lock().expect("profile totals");
    let Some(totals) = guard.as_ref() else { return };
    if totals.is_empty() {
        return;
    }

    let mut rows: Vec<(&&str, &(u64, Duration, Duration))> = totals.iter().collect();
    rows.sort_by(|a, b| b.1 .1.cmp(&a.1 .1));

    let wall = total_wall.as_secs_f64();
    println!();
    println!("self column sums to at most the total. 'inclusive' double-counts nesting.");
    println!();
    println!(
        "{:<46}{:>10}{:>10}{:>8}{:>10}",
        "stage", "calls", "self_s", "self_%", "incl_s"
    );
    println!("{}", "-".repeat(84));
    let mut accounted = 0.0;
    for (label, (calls, self_t, incl)) in &rows {
        let self_s = self_t.as_secs_f64();
        accounted += self_s;
        println!(
            "{:<46}{:>10}{:>10.2}{:>8.1}{:>10.2}",
            label,
            calls,
            self_s,
            if wall > 0.0 {
                100.0 * self_s / wall
            } else {
                0.0
            },
            incl.as_secs_f64()
        );
    }
    println!("{}", "-".repeat(84));
    println!(
        "{:<46}{:>10}{:>10.2}{:>8.1}",
        "accounted for",
        "",
        accounted,
        if wall > 0.0 {
            100.0 * accounted / wall
        } else {
            0.0
        }
    );
    println!(
        "{:<46}{:>10}{:>10.2}{:>8.1}",
        "unattributed (I/O, parsing, glue)",
        "",
        wall - accounted,
        if wall > 0.0 {
            100.0 * (wall - accounted) / wall
        } else {
            0.0
        }
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scopes_are_free_when_profiling_is_off() {
        // The default in tests: no env var, so `scope` hands back None and
        // nothing is recorded.
        assert!(scope("anything").is_none());
    }
}
