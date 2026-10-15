## Drug-name normalisation shared by steps 8 and 10 (sourced, not run on its own).
## Key = upper case, trimmed, trailing salt / hydrate / counter-ion words removed, so DGIdb "AMLODIPINE BESYLATE"
## and iLINCS "Amlodipine" share the key AMLODIPINE. Only counter-ions and hydration states are listed: esters,
## prodrugs, formulations (liposomal, conjugates) and stereo descriptors are different entities and are NOT removed.
SALT_WORDS <- paste0("(DIHYDROCHLORIDE|HYDROCHLORIDE|HYDROBROMIDE|HCL|MESYLATE|MESILATE|BESYLATE|BESILATE|",
                     "MALEATE|FUMARATE|CITRATE|SULFATE|SULPHATE|PHOSPHATE|TARTRATE|ACETATE|SUCCINATE|",
                     "TOSYLATE|TOSILATE|LACTATE|GLUCONATE|NITRATE|SODIUM|POTASSIUM|CALCIUM|MAGNESIUM|",
                     "CHLORIDE|BROMIDE|ISETHIONATE|ANHYDROUS|MONOHYDRATE|DIHYDRATE|TRIHYDRATE|HEMIHYDRATE|HYDRATE)")
drug_key <- function(x) {
  x <- toupper(trimws(x))
  repeat {
    y <- trimws(sub(paste0("[[:space:],()-]+", SALT_WORDS, "$"), "", x))
    if (identical(y, x)) break
    x <- y
  }
  x
}
## DGIdb returns its own entity names, which can be a related but different entity (a liposomal form, a
## conjugate, a prodrug ester). Each DGIdb name is accepted only if its key equals the key of an iLINCS query name;
## the matching query name(s) are kept. Returns one row per DGIdb name.
match_query_names <- function(dgidb_names, query_names) {
  dn <- unique(dgidb_names[!is.na(dgidb_names)])
  qk <- split(query_names, drug_key(query_names))
  q  <- vapply(drug_key(dn), function(k) if (k %in% names(qk)) paste(sort(unique(qk[[k]])), collapse = ";") else NA_character_, character(1))
  data.frame(Drug = dn, Name_key = drug_key(dn), Query_name = unname(q),
             Name_match = ifelse(is.na(q), "unmatched", "normalised name"), stringsAsFactors = FALSE)
}
