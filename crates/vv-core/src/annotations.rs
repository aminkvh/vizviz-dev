//! File header metadata (title, method, resolution, citation, ...).
//!
//! Stored generically as mmCIF-shaped categories, so keeping one more field
//! is a parser change, not a data-model change. PDB header records are
//! mapped onto the same category/item names by the reader, so consumers
//! see one vocabulary regardless of the source format.

/// One mmCIF category (or a PDB header record mapped onto mmCIF names):
/// column names plus rows of values, in file order. Single-value categories
/// have one row. Null values (`?`, `.`) are stored as empty strings.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AnnotationCategory {
    /// Category name without the leading underscore, e.g. `"struct"`.
    pub name: String,
    /// Item names without the category prefix, e.g. `["title"]`.
    pub items: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl AnnotationCategory {
    /// Column index of `item`; case-insensitive like the CIF reader.
    pub fn column(&self, item: &str) -> Option<usize> {
        self.items.iter().position(|i| i.eq_ignore_ascii_case(item))
    }

    /// Value at (`item`, `row`), or `None` when missing or null (empty).
    pub fn get(&self, item: &str, row: usize) -> Option<&str> {
        let col = self.column(item)?;
        let v = self.rows.get(row)?.get(col)?.as_str();
        (!v.is_empty()).then_some(v)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Annotations {
    pub categories: Vec<AnnotationCategory>,
}

impl Annotations {
    pub fn category(&self, name: &str) -> Option<&AnnotationCategory> {
        self.categories
            .iter()
            .find(|c| c.name.eq_ignore_ascii_case(name))
    }

    /// First-row value of `category.item`, e.g. `get("struct", "title")`.
    /// Missing and null both read as `None`.
    pub fn get(&self, category: &str, item: &str) -> Option<&str> {
        self.category(category)?.get(item, 0)
    }

    pub fn is_empty(&self) -> bool {
        self.categories.is_empty()
    }

    // Convenience accessors that know the usual places.

    /// `struct.title`
    pub fn title(&self) -> Option<&str> {
        self.get("struct", "title")
    }

    /// `exptl.method`
    pub fn method(&self) -> Option<&str> {
        self.get("exptl", "method")
    }

    /// `refine.ls_d_res_high` (crystallography), else
    /// `em_3d_reconstruction.resolution` (cryo-EM). Angstroms.
    pub fn resolution(&self) -> Option<f32> {
        self.get("refine", "ls_d_res_high")
            .or_else(|| self.get("em_3d_reconstruction", "resolution"))
            .and_then(|s| s.trim().parse().ok())
    }

    /// `pdbx_database_status.recvd_initial_deposition_date`, as written in
    /// the file (`1984-03-07` in mmCIF, `07-MAR-84` in PDB).
    pub fn deposition_date(&self) -> Option<&str> {
        self.get("pdbx_database_status", "recvd_initial_deposition_date")
    }

    /// `struct_keywords.pdbx_keywords`
    pub fn keywords(&self) -> Option<&str> {
        self.get("struct_keywords", "pdbx_keywords")
    }

    /// Source organism of the first entity: `entity_src_gen` (expressed)
    /// or `entity_src_nat` (natural source).
    pub fn organism(&self) -> Option<&str> {
        self.get("entity_src_gen", "pdbx_gene_src_scientific_name")
            .or_else(|| self.get("entity_src_nat", "pdbx_organism_scientific"))
    }

    /// `citation.title` of the primary citation (row 0).
    pub fn citation_title(&self) -> Option<&str> {
        self.get("citation", "title")
    }

    /// `citation.pdbx_database_id_DOI` of the primary citation (row 0).
    pub fn doi(&self) -> Option<&str> {
        self.get("citation", "pdbx_database_id_DOI")
    }

    /// (entity id, description) pairs from `entity.id` /
    /// `entity.pdbx_description`. Rows without an id are skipped.
    pub fn entities(&self) -> Vec<(String, String)> {
        let Some(cat) = self.category("entity") else {
            return Vec::new();
        };
        (0..cat.rows.len())
            .filter_map(|row| {
                let id = cat.get("id", row)?;
                let desc = cat.get("pdbx_description", row).unwrap_or("");
                Some((id.to_string(), desc.to_string()))
            })
            .collect()
    }

    /// UniProt accessions from `struct_ref` rows whose `db_name` is `UNP`,
    /// one per row (a PDB file lists one DBREF per chain, so repeats happen).
    pub fn uniprot_accessions(&self) -> Vec<String> {
        let Some(cat) = self.category("struct_ref") else {
            return Vec::new();
        };
        (0..cat.rows.len())
            .filter(|&row| cat.get("db_name", row).is_some_and(|db| db == "UNP"))
            .filter_map(|row| cat.get("pdbx_db_accession", row))
            .map(str::to_string)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cat(name: &str, items: &[&str], rows: &[&[&str]]) -> AnnotationCategory {
        AnnotationCategory {
            name: name.to_string(),
            items: items.iter().map(|s| s.to_string()).collect(),
            rows: rows
                .iter()
                .map(|r| r.iter().map(|s| s.to_string()).collect())
                .collect(),
        }
    }

    #[test]
    fn get_reads_first_row_and_treats_null_as_absent() {
        let a = Annotations {
            categories: vec![cat("struct", &["entry_id", "title"], &[&["1ABC", "Hello"]])],
        };
        assert_eq!(a.get("struct", "title"), Some("Hello"));
        assert_eq!(a.get("STRUCT", "TITLE"), Some("Hello"));
        assert_eq!(a.title(), Some("Hello"));
        assert_eq!(a.get("struct", "missing"), None);
        assert_eq!(a.get("nope", "title"), None);

        let nulls = Annotations {
            categories: vec![cat("exptl", &["method"], &[&[""]])],
        };
        assert_eq!(nulls.method(), None);
        assert!(!nulls.is_empty());
    }

    #[test]
    fn resolution_prefers_refine_then_falls_back_to_em() {
        let xray = Annotations {
            categories: vec![
                cat("refine", &["ls_d_res_high"], &[&["1.74"]]),
                cat("em_3d_reconstruction", &["resolution"], &[&["3.2"]]),
            ],
        };
        assert_eq!(xray.resolution(), Some(1.74));

        let em = Annotations {
            categories: vec![
                cat("refine", &["ls_d_res_high"], &[&[""]]),
                cat("em_3d_reconstruction", &["resolution"], &[&["3.2"]]),
            ],
        };
        assert_eq!(em.resolution(), Some(3.2));

        let none = Annotations {
            categories: vec![cat("refine", &["ls_d_res_high"], &[&["n/a"]])],
        };
        assert_eq!(none.resolution(), None);
    }

    #[test]
    fn entities_pair_ids_with_descriptions() {
        let a = Annotations {
            categories: vec![cat(
                "entity",
                &["id", "type", "pdbx_description"],
                &[
                    &["1", "polymer", "Hemoglobin subunit alpha"],
                    &["2", "water", ""],
                    &["", "polymer", "orphan row"],
                ],
            )],
        };
        assert_eq!(
            a.entities(),
            vec![
                ("1".to_string(), "Hemoglobin subunit alpha".to_string()),
                ("2".to_string(), String::new()),
            ]
        );
        assert!(Annotations::default().entities().is_empty());
    }

    #[test]
    fn uniprot_accessions_filter_on_db_name() {
        let a = Annotations {
            categories: vec![cat(
                "struct_ref",
                &["id", "db_name", "pdbx_db_accession"],
                &[
                    &["1", "UNP", "P69905"],
                    &["2", "PDB", "4HHB"],
                    &["3", "UNP", "P68871"],
                    &["4", "UNP", ""],
                ],
            )],
        };
        assert_eq!(a.uniprot_accessions(), vec!["P69905", "P68871"]);
    }

    #[test]
    fn empty_annotations() {
        let a = Annotations::default();
        assert!(a.is_empty());
        assert_eq!(a.title(), None);
        assert_eq!(a.resolution(), None);
        assert!(a.uniprot_accessions().is_empty());
    }
}
