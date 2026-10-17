# automatic trait selection preserves candidate priority

    Code
      .infer_gc_trait_col(df)
    Message
      i Multiple glycan trait columns found: glycan_structure, glycan_composition, trait, and motif.
      i Using glycan_structure. Set `by = "<column>"` to select another column.
    Output
      [1] "glycan_structure"
    Code
      .infer_gc_trait_col(df, prefer_trait = TRUE)
    Message
      i Multiple glycan trait columns found: glycan_structure, glycan_composition, trait, and motif.
      i Using trait. Set `by = "<column>"` to select another column.
    Output
      [1] "trait"

# trait selection rejects invalid names and grouping columns clearly

    Code
      .infer_gc_trait_col(df, by = "missing")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! The column "missing" selected by `by` is missing from `dea_res`.

---

    Code
      .infer_gc_trait_col(df, by = character())
    Condition
      Error in `.infer_gc_trait_col()`:
      ! `by` must be NULL or a single, non-empty column name.

---

    Code
      .infer_gc_trait_col(df, by = c("trait", "motif"))
    Condition
      Error in `.infer_gc_trait_col()`:
      ! `by` must be NULL or a single, non-empty column name.

---

    Code
      .infer_gc_trait_col(df, by = NA_character_)
    Condition
      Error in `.infer_gc_trait_col()`:
      ! `by` must be NULL or a single, non-empty column name.

---

    Code
      .infer_gc_trait_col(df, by = "")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! `by` must be NULL or a single, non-empty column name.

---

    Code
      .infer_gc_trait_col(df, by = 1)
    Condition
      Error in `.infer_gc_trait_col()`:
      ! `by` must be NULL or a single, non-empty column name.

---

    Code
      .infer_gc_trait_col(tibble::tibble(protein = "P1"))
    Condition
      Error in `.infer_gc_trait_col()`:
      ! Cannot determine glycan traits.
      i At least one of these columns is needed: glycan_structure, glycan_composition, trait, and motif
      i Did you accidentally set `add_info` to `FALSE` when performing DEA with glystats?

---

    Code
      .infer_gc_trait_col(df, by = "trait")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! The `by` column trait cannot be used to group glycan traits.
      i Use a one-dimensional atomic vector, such as character or factor labels.

---

    Code
      .infer_gc_trait_col(df, by = "trait")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! The `by` column trait cannot be used to group glycan traits.
      i Use a one-dimensional atomic vector, such as character or factor labels.

---

    Code
      .infer_gc_trait_col(df, by = "trait")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! The `by` column trait must have a unique name.

