# glycan-centric ORA rejects invalid by before enrichment

    Code
      enrich_gc_ora_kegg(df, by = "missing")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! The column "missing" selected by `by` is missing from `dea_res`.

---

    Code
      enrich_gc_ora_kegg(df, by = "trait")
    Condition
      Error in `.infer_gc_trait_col()`:
      ! The `by` column trait cannot be used to group glycan traits.
      i Use a one-dimensional atomic vector, such as character or factor labels.

# enrich_gc_ora_go errors on data.frame with missing columns

    Code
      suppressMessages(enrich_gc_ora_go(dea_res_missing_trait))
    Condition
      Error in `.infer_gc_trait_col()`:
      ! Cannot determine glycan traits.
      i At least one of these columns is needed: glycan_structure, glycan_composition, trait, and motif
      i Did you accidentally set `add_info` to `FALSE` when performing DEA with glystats?

