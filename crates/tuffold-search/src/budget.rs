//! The memory Sage's index will take, counted before it is built.
//!
//! Sage builds its index in one go: every digest of every protein, times
//! every modification combination, times two for the reversed decoys, is
//! collected into one list before duplicates are merged, and then every
//! peptide's fragments into another. A semi-specific search, or many
//! missed cleavages or variable modifications, can ask for more memory than
//! the computer has, and the operating system swaps or kills the process
//! instead of saying so. So the same digest is first only counted, protein
//! by protein, with Sage's own code and nothing kept; the search is refused
//! when the estimate exceeds its memory budget. Counting stops as soon as
//! the budget is passed, so a hopeless search is refused quickly.

use std::mem::size_of;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;
use sage_core::database::{Parameters, Theoretical};
use sage_core::enzyme::{Digest, EnzymeParameters};
use sage_core::peptide::Peptide;

/// Allocator bookkeeping and rounding per heap allocation.
const ALLOCATION: u64 = 16;
/// A reference-counted allocation's two counters.
const ARC_HEADER: u64 = 16;
/// Rayon collects a parallel iterator into pieces and then joins them, so
/// for a while both exist.
const COLLECT_COPIES: u64 = 2;
/// Bytes per entry of the set of target sequences Sage builds to drop
/// decoys that equal a target: the bytes, the hash table slot, the lock
/// shards' share.
const SET_ENTRY: u64 = 48;

/// The share of the memory free when a search starts that its index may
/// take; the rest is for one run's spectra, the project, and everything else
/// running.
const FREE_SHARE: f64 = 0.75;
/// The budget when the system does not say how much memory it has.
const UNKNOWN_SYSTEM_BUDGET: u64 = 4_000_000_000;

/// A search's index budget when none is set: three quarters of the memory
/// free right now, else half of all memory if the free amount is unknown.
pub fn default_memory_budget() -> u64 {
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    let free = match (free_memory().unwrap_or(0), system.total_memory()) {
        (0, 0) => return UNKNOWN_SYSTEM_BUDGET,
        (0, total) => total / 2,
        (free, _) => free,
    };
    (free as f64 * FREE_SHARE) as u64
}

/// Memory the system can hand out without swapping. sysinfo's figure on
/// macOS subtracts the pages the compressor occupies from free, inactive,
/// and purgeable pages, but those pages are in none of the three, so they
/// were taken off twice: with 5.8 GB compressed it said 0 and the budget
/// fell back to half of all memory. The counts come from `vm_stat`, since
/// the kernel call needs `unsafe` (forbidden here).
#[cfg(target_os = "macos")]
fn free_memory() -> Option<u64> {
    let output = std::process::Command::new("/usr/bin/vm_stat")
        .output()
        .ok()
        .filter(|output| output.status.success())?;
    let (counts, page_size) = parse_vm_stat(&String::from_utf8_lossy(&output.stdout))?;
    Some(reclaimable_pages(&counts).saturating_mul(page_size))
}

#[cfg(not(target_os = "macos"))]
fn free_memory() -> Option<u64> {
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    Some(system.available_memory())
}

/// The page counts the budget uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(not(target_os = "macos"), allow(dead_code))]
struct VmCounts {
    /// Free pages, speculative ones (read ahead, not yet used) included.
    free: u64,
    /// Not touched lately: file pages are dropped, the rest compressed.
    inactive: u64,
    /// Volatile caches the kernel may empty.
    purgeable: u64,
}

/// Free, inactive, and purgeable pages. Active pages could be compressed
/// too, and the kernel's own level (`kern.memorystatus_level`) counts them,
/// but taking them means compressing what other programs are using now.
#[cfg_attr(not(target_os = "macos"), allow(dead_code))]
fn reclaimable_pages(counts: &VmCounts) -> u64 {
    counts
        .free
        .saturating_add(counts.inactive)
        .saturating_add(counts.purgeable)
}

/// The counts and page size in `vm_stat`'s output, or `None` if any is
/// missing. Its "Pages free" leaves out the speculative pages, which are
/// just as free.
#[cfg_attr(not(target_os = "macos"), allow(dead_code))]
fn parse_vm_stat(text: &str) -> Option<(VmCounts, u64)> {
    let page_size = text
        .split("page size of ")
        .nth(1)?
        .split_whitespace()
        .next()?
        .parse()
        .ok()?;
    let count = |label: &str| -> Option<u64> {
        text.lines()
            .find_map(|line| line.strip_prefix(label)?.strip_prefix(':'))?
            .trim()
            .trim_end_matches('.')
            .parse()
            .ok()
    };
    let counts = VmCounts {
        free: count("Pages free")?.saturating_add(count("Pages speculative")?),
        inactive: count("Pages inactive")?,
        purgeable: count("Pages purgeable")?,
    };
    Some((counts, page_size))
}

/// Counting goes on past the budget up to this multiple of it, so a refusal
/// can say how much the search would really need, yet a search that could
/// never fit (non-specific, say) is still refused in seconds.
pub(crate) const COUNT_PAST_BUDGET: u64 = 4;

/// What the index will hold and the most memory building it takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct IndexEstimate {
    /// Target digests (one per protein a sequence occurs in).
    pub digests: u64,
    /// Target and decoy peptide forms before duplicates are merged.
    pub forms: u64,
    /// Theoretical fragments, at most.
    pub fragments: u64,
    /// The largest amount of memory building the index holds at once.
    pub peak_bytes: u64,
    /// False when counting stopped early because the budget was passed:
    /// the numbers are then only what was counted so far.
    pub complete: bool,
}

/// Running totals, shared by the counting threads.
#[derive(Default)]
struct Totals {
    digests: AtomicU64,
    digest_bytes: AtomicU64,
    forms: AtomicU64,
    form_bytes: AtomicU64,
    fragments: AtomicU64,
}

impl Totals {
    fn estimate(&self, complete: bool) -> IndexEstimate {
        let digests = self.digests.load(Ordering::Relaxed);
        let digest_bytes = self.digest_bytes.load(Ordering::Relaxed);
        let forms = self.forms.load(Ordering::Relaxed);
        let form_bytes = self.form_bytes.load(Ordering::Relaxed);
        let fragments = self.fragments.load(Ordering::Relaxed);
        // Digesting: the digests, the set of target sequences, and the
        // peptide forms being collected.
        let digesting = digest_bytes
            .saturating_add(digests.saturating_mul(SET_ENTRY))
            .saturating_add(form_bytes.saturating_mul(COLLECT_COPIES));
        // Indexing: the (at most this many) peptides and their fragments
        // being collected.
        let fragment_bytes = fragments
            .saturating_mul(size_of::<Theoretical>() as u64)
            .saturating_mul(COLLECT_COPIES);
        let indexing = form_bytes.saturating_add(fragment_bytes);
        IndexEstimate {
            digests,
            forms,
            fragments,
            peak_bytes: digesting.max(indexing),
            complete,
        }
    }
}

/// Count what `params.build(fasta)` would hold for these target proteins,
/// stopping once the peak passes `limit` bytes (then `complete` is false).
pub(crate) fn estimate(
    proteins: &[(Arc<String>, String)],
    params: &Parameters,
    limit: u64,
) -> IndexEstimate {
    let enzyme: EnzymeParameters = params.enzyme.clone().into();
    let mods: Vec<_> = params
        .variable_mods
        .iter()
        .flat_map(|(site, masses)| masses.iter().map(move |mass| (*site, *mass)))
        .collect();
    let totals = Totals::default();
    let within = proteins.par_iter().try_for_each(|(accession, sequence)| {
        for digest in enzyme.digest(sequence, accession.clone()) {
            count_digest(digest, &mods, params, &totals);
        }
        // Totals only grow, so once over the limit, always over it.
        if totals.estimate(false).peak_bytes > limit {
            Err(())
        } else {
            Ok(())
        }
    });
    totals.estimate(within.is_ok())
}

fn count_digest(
    digest: Digest,
    mods: &[(sage_core::modification::ModificationSpecificity, f32)],
    params: &Parameters,
    totals: &Totals,
) {
    let length = digest.sequence.len() as u64;
    totals.digests.fetch_add(1, Ordering::Relaxed);
    totals.digest_bytes.fetch_add(
        size_of::<Digest>() as u64 + length + ALLOCATION,
        Ordering::Relaxed,
    );
    let Ok(peptide) = Peptide::try_from(digest) else {
        return;
    };
    let kept = peptide
        .apply(mods, &params.static_mods, params.max_variable_mods)
        .into_iter()
        .filter(|form| {
            form.monoisotopic >= params.peptide_min_mass
                && form.monoisotopic <= params.peptide_max_mass
        })
        .count() as u64;
    if kept == 0 {
        return;
    }
    // Each kept form and its reversed decoy.
    let forms = kept * 2;
    // A form: the struct, its per-residue modification masses, its protein
    // list (one entry before duplicates merge); the target forms of one
    // digest share a sequence, every decoy has its own.
    let per_form = size_of::<Peptide>() as u64
        + (4 * length + ALLOCATION)
        + (size_of::<Arc<String>>() as u64 + ALLOCATION);
    let sequences = (1 + kept) * (ARC_HEADER + length + ALLOCATION);
    totals.forms.fetch_add(forms, Ordering::Relaxed);
    totals
        .form_bytes
        .fetch_add(forms * per_form + sequences, Ordering::Relaxed);
    // b and y ions past the first `min_ion_index`, before the m/z filter.
    let ions = length
        .saturating_sub(1)
        .saturating_sub(params.min_ion_index as u64)
        * params.ion_kinds.len() as u64;
    totals.fragments.fetch_add(forms * ions, Ordering::Relaxed);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `vm_stat` on an 18 GB Mac (abridged; the full output has more
    /// lines between these).
    const VM_STAT: &str = "\
Mach Virtual Memory Statistics: (page size of 16384 bytes)
Pages free:                                     9864.
Pages active:                                 284254.
Pages inactive:                               281979.
Pages speculative:                              1598.
Pages wired down:                             193845.
Pages purgeable:                                4825.
Pages stored in compressor:                  1443504.
Pages occupied by compressor:                 354717.
";

    #[test]
    fn compressed_pages_are_not_taken_off_free_memory() {
        // free = 9,864 + 1,598 speculative = 11,462; reclaimable = 11,462 +
        // 281,979 inactive + 4,825 purgeable = 298,266 pages × 16,384 =
        // 4,886,790,144 bytes. sysinfo took the 354,717 pages occupied by
        // the compressor off that too, leaving 0.
        let (counts, page_size) = parse_vm_stat(VM_STAT).unwrap();
        assert_eq!(
            counts,
            VmCounts {
                free: 11_462,
                inactive: 281_979,
                purgeable: 4_825
            }
        );
        assert_eq!(page_size, 16_384);
        assert_eq!(reclaimable_pages(&counts) * page_size, 4_886_790_144);
    }

    #[test]
    fn vm_stat_missing_a_count_gives_nothing() {
        let without_inactive = VM_STAT.replace("Pages inactive", "Pages idle");
        assert_eq!(parse_vm_stat(&without_inactive), None);
        assert_eq!(parse_vm_stat("Pages free: 1."), None);
    }

    #[test]
    fn free_memory_is_known_and_below_the_total() {
        let mut system = sysinfo::System::new();
        system.refresh_memory();
        let free = free_memory().unwrap();
        assert!(free > 0 && free <= system.total_memory(), "{free}");
        assert!(default_memory_budget() <= free);
    }
}
