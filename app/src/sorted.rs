//! Rows as loaded plus the order they are shown in, for sortable tables.
//! Sorting permutes indices (rows are never copied or dropped), ties keep
//! the load order, and clearing the sort restores it. New rows keep the
//! current sort.

use std::cmp::Ordering;

/// Which way a column is sorted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Ascending,
    Descending,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SortedRows<T> {
    loaded: Vec<T>,
    /// `order[shown] = loaded index`.
    order: Vec<usize>,
    sort: Option<(usize, Direction)>,
}

impl<T> Default for SortedRows<T> {
    fn default() -> Self {
        Self {
            loaded: Vec::new(),
            order: Vec::new(),
            sort: None,
        }
    }
}

impl<T> SortedRows<T> {
    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// The row shown at `ix`.
    pub fn get(&self, ix: usize) -> Option<&T> {
        self.order.get(ix).map(|&loaded| &self.loaded[loaded])
    }

    /// Rows in the order shown.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.order.iter().map(|&loaded| &self.loaded[loaded])
    }

    /// Where the first row matching `predicate` is shown.
    pub fn position(&self, predicate: impl Fn(&T) -> bool) -> Option<usize> {
        self.iter().position(predicate)
    }

    /// The direction `col` is sorted in, if it is the sorted column.
    pub fn direction(&self, col: usize) -> Option<Direction> {
        self.sort.filter(|(c, _)| *c == col).map(|(_, d)| d)
    }

    /// Replace the rows, keeping the current sort.
    pub fn replace(&mut self, rows: Vec<T>, compare: impl Fn(usize, &T, &T) -> Compared) {
        self.loaded = rows;
        self.apply(compare);
    }

    /// Sort by `col` (`compare` orders two rows ascending by a column), or
    /// back to the load order with `None`.
    pub fn sort_by(
        &mut self,
        sort: Option<(usize, Direction)>,
        compare: impl Fn(usize, &T, &T) -> Compared,
    ) {
        self.sort = sort;
        self.apply(compare);
    }

    /// Sort as [`Self::sort_by`] and return where the row shown at
    /// `selected` moved to, so a selection follows its row.
    pub fn sort_keeping(
        &mut self,
        sort: Option<(usize, Direction)>,
        compare: impl Fn(usize, &T, &T) -> Compared,
        selected: Option<usize>,
    ) -> Option<usize> {
        let loaded = selected.and_then(|ix| self.order.get(ix).copied());
        self.sort_by(sort, compare);
        loaded.and_then(|loaded| self.order.iter().position(|&l| l == loaded))
    }

    fn apply(&mut self, compare: impl Fn(usize, &T, &T) -> Compared) {
        let mut order: Vec<usize> = (0..self.loaded.len()).collect();
        if let Some((col, direction)) = self.sort {
            // Stable: ties keep the load order in both directions.
            order.sort_by(|&a, &b| {
                let (a, b) = (&self.loaded[a], &self.loaded[b]);
                let compared = compare(col, a, b);
                let value = match direction {
                    Direction::Ascending => compared.value,
                    Direction::Descending => compared.value.reverse(),
                };
                // Missing values stay last whichever way the column sorts.
                compared.missing.then(value)
            });
        }
        self.order = order;
    }
}

/// How two rows compare on a column, ascending: `missing` puts rows
/// without a value after those with one (in either direction), `value`
/// orders the known values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Compared {
    pub missing: Ordering,
    pub value: Ordering,
}

impl From<Ordering> for Compared {
    fn from(value: Ordering) -> Self {
        Self {
            missing: Ordering::Equal,
            value,
        }
    }
}

/// Compare optional values; missing ones sort after every known one.
pub fn missing_last<V>(
    a: Option<V>,
    b: Option<V>,
    compare: impl Fn(&V, &V) -> Ordering,
) -> Compared {
    match (a, b) {
        (Some(a), Some(b)) => compare(&a, &b).into(),
        (a, b) => Compared {
            missing: a.is_none().cmp(&b.is_none()),
            value: Ordering::Equal,
        },
    }
}

/// Compare text with a trailing number numerically, so `scan=9` comes
/// before `scan=10` and `P2` before `P10`.
pub fn natural(a: &str, b: &str) -> Ordering {
    let ((head_a, number_a), (head_b, number_b)) = (number_tail(a), number_tail(b));
    head_a
        .cmp(head_b)
        .then_with(|| number_a.cmp(&number_b))
        .then_with(|| a.cmp(b))
}

fn number_tail(s: &str) -> (&str, Option<u128>) {
    let head = s.trim_end_matches(|c: char| c.is_ascii_digit());
    (head, s[head.len()..].parse().ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn by_value(col: usize, a: &(i32, &str), b: &(i32, &str)) -> Compared {
        match col {
            0 => a.0.cmp(&b.0).into(),
            _ => a.1.cmp(b.1).into(),
        }
    }

    fn shown(rows: &SortedRows<(i32, &'static str)>) -> Vec<&'static str> {
        rows.iter().map(|r| r.1).collect()
    }

    #[test]
    fn sorting_permutes_and_clearing_restores_the_load_order() {
        let mut rows = SortedRows::default();
        rows.replace(vec![(2, "b"), (1, "c"), (2, "a"), (0, "d")], by_value);
        assert_eq!(shown(&rows), ["b", "c", "a", "d"]);

        rows.sort_by(Some((0, Direction::Ascending)), by_value);
        // Ties (b and a, both 2) keep the load order.
        assert_eq!(shown(&rows), ["d", "c", "b", "a"]);
        rows.sort_by(Some((0, Direction::Descending)), by_value);
        assert_eq!(shown(&rows), ["b", "a", "c", "d"]);
        assert_eq!(rows.direction(0), Some(Direction::Descending));
        assert_eq!(rows.direction(1), None);
        assert_eq!(rows.get(3), Some(&(0, "d")));
        assert_eq!(rows.position(|r| r.1 == "c"), Some(2));

        rows.sort_by(None, by_value);
        assert_eq!(shown(&rows), ["b", "c", "a", "d"]);
    }

    #[test]
    fn new_rows_keep_the_current_sort() {
        let mut rows = SortedRows::default();
        rows.sort_by(Some((1, Direction::Ascending)), by_value);
        rows.replace(vec![(1, "z"), (2, "x"), (3, "y")], by_value);
        assert_eq!(shown(&rows), ["x", "y", "z"]);
        assert_eq!(rows.len(), 3);
        rows.replace(Vec::new(), by_value);
        assert!(rows.is_empty() && rows.get(0).is_none());
    }

    #[test]
    fn a_selection_follows_its_row() {
        let mut rows = SortedRows::default();
        rows.replace(vec![(3, "a"), (1, "b"), (2, "c")], by_value);
        // "a" is shown first; sorted by value it moves to the end.
        let moved = rows.sort_keeping(Some((0, Direction::Ascending)), by_value, Some(0));
        assert_eq!(moved, Some(2));
        assert_eq!(rows.get(2), Some(&(3, "a")));
        assert_eq!(rows.sort_keeping(None, by_value, None), None);
        assert_eq!(rows.sort_keeping(None, by_value, Some(9)), None);
    }

    #[test]
    fn trailing_numbers_compare_as_numbers() {
        assert_eq!(natural("scan=9", "scan=10"), Ordering::Less);
        assert_eq!(natural("scan=10", "scan=10"), Ordering::Equal);
        assert_eq!(natural("index=5", "scan=1"), Ordering::Less);
        assert_eq!(
            natural("P02", "P2"),
            Ordering::Less,
            "equal numbers fall back to text"
        );
        assert_eq!(natural("run", "run1"), Ordering::Less);
    }

    #[test]
    fn missing_values_sort_last_in_both_directions() {
        let by_q = |_: usize, a: &(Option<f64>, &str), b: &(Option<f64>, &str)| {
            missing_last(a.0, b.0, f64::total_cmp)
        };
        let names = |rows: &SortedRows<(Option<f64>, &'static str)>| {
            rows.iter().map(|r| r.1).collect::<Vec<_>>()
        };
        let mut rows = SortedRows::default();
        rows.replace(
            vec![(None, "x"), (Some(0.2), "b"), (Some(0.1), "a"), (None, "y")],
            by_q,
        );
        rows.sort_by(Some((0, Direction::Ascending)), by_q);
        assert_eq!(names(&rows), ["a", "b", "x", "y"]);
        rows.sort_by(Some((0, Direction::Descending)), by_q);
        assert_eq!(names(&rows), ["b", "a", "x", "y"]);
    }
}
