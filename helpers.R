library(data.table)

find_overlaps <- function(
    # Find genomic overlaps based on start, end columns

    # wrapper around data.table::foverlaps, should yield same result as
    # IRanges::findOverlaps(plyranges::as_iranges(tb1), plyranges::as_iranges(tb2))
    tb1,
    tb2,
    # contrast with "any", which allows partial overlaps
    overlap_type = "within"
) {
    tb1_dt <- data.table::as.data.table(tb1)
    data.table::setkey(tb1_dt, start, end)
    tb2_dt <- data.table::as.data.table(dplyr::distinct(tb2, start, end))
    data.table::setkey(tb2_dt, start, end)
    res <- data.table::foverlaps(
        tb1_dt,
        tb2_dt,
        nomatch = NULL,
        type = overlap_type
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::select(-c(start, end)) %>%
        dplyr::rename(
            start=i.start,
            end=i.end
            ) %>%
        dplyr::distinct()
    return(res)
}