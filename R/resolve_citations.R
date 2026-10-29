# Resolve each citation in the bioassay databases to the paper it refers to
# (supplement, "Bioassay data processing"), so that the prep can tell papers
# apart by more than first author and year. Run by hand when a database is
# re-downloaded: it reads the citations the prep lists
# (outputs/review/bioassay_citations.csv), looks up those not yet in
# review/bioassay_citations.csv, and adds them to it. Rows decided by manual
# review (`method` "manual", or "report" for a citation that is a report) are
# left alone.
#
# Each citation is looked up by the best evidence it gives, in this order: a
# DOI in its text or link; a PubMed or PMC ID; a DOI on the page its link
# points to; a search of Crossref and Europe PMC by its title; the first page
# of the PDF its link points to; or else a search by its first author and year
# with words of the field. A match must agree with the
# citation's first author and year (within one year, for online and print
# dates). A preprint is resolved to the paper it became, where there is one.

# Needs the httr2 package (R/packages.R) and the pdftotext program (poppler)
# to read linked PDFs.

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_prep_functions.R")

user_agent <- "ir_cube bioassay prep (https://github.com/idem-lab/ir_cube)"

# Lookups are cached on disk (outputs/cache/citations), so a re-run repeats
# none; a failed request isn't cached, and is tried again on the next run
cache_dir <- "outputs/cache/citations"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cached <- function(url, fetch) {
  file <- file.path(cache_dir, paste0(rlang::hash(url), ".rds"))
  if (file.exists(file)) return(readRDS(file)$value)
  value <- fetch(url)
  if (!identical(attr(value, "failed"), TRUE)) saveRDS(list(value = value), file)
  value
}
failed <- structure(list(), failed = TRUE)

get_json <- function(url) {
  value <- cached(url, get_json_now)
  if (identical(attr(value, "failed"), TRUE)) NULL else value
}
get_page <- function(url) {
  value <- cached(url, get_page_now)
  if (identical(attr(value, "failed"), TRUE)) NA_character_ else value
}

get_json_now <- function(url) {
  response <- tryCatch(
    request(url) %>%
      req_user_agent(user_agent) %>%
      req_timeout(30) %>%
      req_retry(max_tries = 3) %>%
      req_error(is_error = function(resp) FALSE) %>%
      req_perform(),
    error = function(e) NULL)
  if (is.null(response) || resp_status(response) >= 500 ||
        resp_status(response) == 429) return(failed)
  if (resp_status(response) != 200) return(NULL)
  tryCatch(resp_body_json(response), error = function(e) NULL)
}

get_page_now <- function(url) {
  response <- tryCatch(
    request(url) %>%
      req_user_agent(user_agent) %>%
      req_timeout(30) %>%
      req_error(is_error = function(resp) FALSE) %>%
      req_perform(),
    error = function(e) NULL)
  if (is.null(response)) return(failed)
  if (resp_status(response) != 200) return(NA_character_)
  type <- resp_content_type(response)
  if (!grepl("html", type)) return(NA_character_)
  tryCatch(resp_body_string(response), error = function(e) NA_character_)
}

# identifiers in a citation's text or link
clean_doi <- function(doi) {
  doi <- tolower(sub("[.,;)\\]]+$", "", doi))
  sub("(/full|/abstract|/pdf|\\.pdf|/epdf|[?&].*)$", "", doi)
}
find_doi <- function(text) {
  doi <- str_match(text, "10\\.\\d{4,9}/[^\\s,;\"<>]+")[, 1]
  # BioMed Central's older links give the DOI's suffix only
  bmc <- str_match(text, paste0("/content/(?:pdf/)?((?:14\\d\\d|17\\d\\d|20\\d\\d)",
                                "-\\d{3}[\\dx]-\\d+-\\d+|s\\d{5}-\\d{3}-\\d{4,5}-[0-9x])"))[, 2]
  # Nature's links give the DOI's suffix only, and Hindawi's its year and number
  nature <- str_match(text, "nature\\.com/articles/(s\\d{5}-\\d{3}-\\d{4,5}-[0-9x])")[, 2]
  hindawi <- str_match(text, "hindawi\\.com/.*?/((?:19|20)\\d{2})/(\\d{5,7})")
  clean_doi(coalesce(doi,
                     if_else(is.na(bmc), NA, paste0("10.1186/", bmc)),
                     if_else(is.na(nature), NA, paste0("10.1038/", nature)),
                     if_else(is.na(hindawi[, 2]), NA,
                             paste0("10.1155/", hindawi[, 2], "/", hindawi[, 3]))))
}
find_pmid <- function(text) {
  coalesce(str_match(text, "pubmed(?:\\.ncbi\\.nlm\\.nih\\.gov)?/(\\d{6,9})")[, 2],
           str_match(text, "PMID:? ?(\\d{6,9})")[, 2],
           str_match(text, "europepmc\\.org/(?:abstract|article)/MED/(\\d{6,9})")[, 2])
}
find_pmc <- function(text) str_match(text, "PMC\\d{5,8}")[, 1]

# the title a citation gives, where it gives one: after "(year)" near the start
# in the Vector Atlas and MTM, and after the authors in the IR Mapper extracts'
# Vancouver style (authors' initials, then a full stop or "(year)"); NA for
# "Author et al. year" citations. A title ends at a
# full stop, but not one of an abbreviation in it ("s. s.", "s.l.", "An.").
citation_title <- function(citation) {
  x <- str_squish(gsub("\u00c2", " ", plain_text(citation)))
  end <- "(?:\\.\\.|(?<!\\b(?:[A-Za-z]|An|Ae|Cx|sp|spp))\\.\\s|$)"
  tries <- cbind(
    str_match(x, paste0("^[^()]{0,80}\\((?:19|20)\\d{2}[a-z]?\\)\\.?\\s+(.+?)", end))[, 2],
    str_match(x, paste0("et al\\.?,?\\s+(.+?)", end))[, 2],
    str_match(x, paste0("^.*?[A-Z]{1,3}(?:\\.|\\s*\\((?:19|20)\\d{2}\\))\\s+(.+?)", end))[, 2])
  apply(tries, 1, function(t) {
    t <- t[!is.na(t) & str_count(t, "\\S+") >= 5]
    if (length(t)) t[1] else NA_character_
  })
}

words_of <- function(x) {
  w <- strsplit(gsub("[^a-z0-9 ]", " ",
                     tolower(stringi::stri_trans_general(x, "Latin-ASCII"))),
                " +")
  lapply(w, function(v) v[nchar(v) > 2])
}
# the share of one title's words in the other's: of the shorter title, or of
# the longer one (stricter: a short generic title doesn't match a longer one)
title_overlap <- function(a, b, of = min) {
  wa <- words_of(a)[[1]]
  wb <- words_of(b)[[1]]
  if (!length(wa) || !length(wb)) return(0)
  length(intersect(wa, wb)) / of(length(unique(wa)), length(unique(wb)))
}

# the first element of a list or vector, or NA
first_of <- function(x) if (length(x)) x[[1]] else NA

# a field of the first element of a list of records, or NA
first_field <- function(x, field) {
  if (length(x) && is.list(x[[1]]) && !is.null(x[[1]][[field]])) x[[1]][[field]] else NA
}

# a paper's metadata, as a one-row tibble
paper_row <- function(doi = NA, pmid = NA, title = NA, first_author = NA,
                      year = NA, journal = NA, author_name = NA) {
  tibble(doi = as.character(doi), pmid = as.character(pmid),
         title = plain_text(as.character(title)), first_author = as.character(first_author),
         author_name = as.character(coalesce(author_name, first_author)),
         year = as.integer(year), journal = as.character(journal))
}

from_crossref <- function(doi) {
  x <- get_json(paste0("https://api.crossref.org/works/", URLencode(doi, reserved = TRUE)))
  if (is.null(x)) return(NULL)
  m <- x$message
  published <- m$`is-preprint-of`
  if (!is.null(m$relation$`is-preprint-of`)) {
    target <- m$relation$`is-preprint-of`[[1]]
    if (identical(target$`id-type`, "doi")) {
      later <- from_crossref(tolower(target$id))
      if (!is.null(later)) return(later)
    }
  }
  year <- first_of(first_of(m$issued$`date-parts`))
  paper_row(doi = tolower(m$DOI), title = first_of(m$title),
            first_author = first_field(m$author, "family"),
            author_name = paste(first_field(m$author, "given"), first_field(m$author, "family")),
            year = year %||% NA, journal = first_of(m$`container-title`))
}

from_europepmc <- function(query) {
  x <- get_json(paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?",
                       "format=json&resultType=core&pageSize=25&query=",
                       URLencode(query, reserved = TRUE)))
  hits <- x$resultList$result
  if (is.null(hits) || !length(hits)) return(NULL)
  bind_rows(lapply(hits, function(h) {
    paper_row(doi = tolower(h$doi %||% NA), pmid = h$pmid %||% NA,
              title = h$title %||% NA,
              author_name = coalesce(first_field(h$authorList$author, "fullName"),
                                     sub(",.*", "", h$authorString %||% NA)),
              first_author = coalesce(first_field(h$authorList$author, "lastName"),
                                      sub(" .*", "", h$authorString %||% NA)),
              year = h$pubYear %||% NA,
              journal = h$journalInfo$journal$title %||% h$bookOrReportDetails$publisher %||% NA) %>%
      mutate(abstract = h$abstractText %||% NA_character_)
  }))
}

from_crossref_search <- function(citation, year) {
  x <- get_json(paste0("https://api.crossref.org/works?rows=10&query.bibliographic=",
                       URLencode(citation, reserved = TRUE),
                       if (!is.na(year)) paste0("&filter=from-pub-date:", year - 1,
                                                ",until-pub-date:", year + 1)))
  items <- x$message$items
  if (is.null(items) || !length(items)) return(NULL)
  bind_rows(lapply(items, function(m) {
    paper_row(doi = tolower(m$DOI), title = first_of(m$title),
              first_author = first_field(m$author, "family"),
              author_name = paste(first_field(m$author, "given"), first_field(m$author, "family")),
              year = first_of(first_of(m$issued$`date-parts`)) %||% NA,
              journal = first_of(m$`container-title`)) %>%
      mutate(abstract = m$abstract %||% NA_character_)
  }))
}

# DOIs of preprint servers (Research Square, bioRxiv and medRxiv,
# Preprints.org, Authorea, OSF)
is_preprint <- function(doi) grepl("^10\\.(21203|1101|20944|22541|31219)/", coalesce(doi, ""))

# the published paper a preprint became, found by its title and first author
published_version <- function(p) {
  hits <- from_crossref_search(paste(p$first_author, p$title), NA)
  if (is.null(hits)) return(NULL)
  hits <- hits %>%
    filter(!is_preprint(doi), name_key(first_author) == name_key(p$first_author)) %>%
    mutate(overlap = vapply(title, function(t) title_overlap(t, p$title), numeric(1))) %>%
    filter(overlap >= 0.8) %>%
    arrange(desc(overlap))
  if (nrow(hits)) hits[1, ] else NULL
}

# candidates that are one paper (Crossref and Europe PMC both list most papers,
# their titles written a little differently) kept once: the first of those
# whose titles share nine tenths of their words
one_per_paper <- function(candidates) {
  keep <- rep(TRUE, nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    if (!keep[i]) next
    for (j in seq_len(nrow(candidates))) {
      if (j > i && keep[j] &&
            title_overlap(candidates$title[i], candidates$title[j]) >= 0.9) {
        keep[j] <- FALSE
      }
    }
  }
  candidates[keep, ]
}

# text with HTML entities decoded and tags dropped (some citations give
# "&lt;I&gt;Anopheles gambiae&lt;/I&gt;" for the italic name)
plain_text <- function(x) {
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&quot;", "\"", x, fixed = TRUE)
  x <- gsub("&#0?39;|&apos;", "'", x)
  str_squish(gsub("</?[a-zA-Z][^>]*>", " ", x))
}

# a notice about a paper (a correction, erratum or retraction), which shares
# its title
is_notice <- function(title, doi) {
  grepl("^(correction|erratum|corrigendum|retraction|expression of concern)",
        coalesce(title, ""), ignore.case = TRUE) |
    grepl("/annotation/", coalesce(doi, ""), fixed = TRUE)
}

# words of the field, at least one of which a candidate found by author and
# year must have in its title or abstract
field_words <- "insecticid|resistan|anophel|pyrethroid|susceptib|mosquito|vector|malaria"

# text decoded twice (e.g. "Nâdri" for "N’dri"), decoded once more where that works
redecode <- function(x) {
  fixed <- suppressWarnings(iconv(x, from = "UTF-8", to = "latin1"))
  ok <- !is.na(fixed) & grepl("[\u00c2\u00c3\u00e2]", x)
  Encoding(fixed) <- "UTF-8"
  ok <- ok & validUTF8(fixed)
  if_else(ok, fixed, x)
}

# whether two surnames are the same, allowing for a misspelling (one letter
# different, two for names of six or more), or one part of a compound name
same_surname <- function(a, b) {
  a <- name_key(redecode(a))
  b <- name_key(redecode(b))
  ok <- !is.na(a) & !is.na(b) & a != "" & b != ""
  close <- vapply(seq_along(a), function(i) {
    x <- a[i]
    y <- b[i]
    if (is.na(x) || is.na(y) || x == "" || y == "") return(FALSE)
    x == y || grepl(x, y, fixed = TRUE) || grepl(y, x, fixed = TRUE) ||
      utils::adist(x, y)[1, 1] <= if (max(nchar(x), nchar(y)) >= 6) 2 else 1
  }, logical(1))
  ok & close
}

# whether a candidate agrees with the citation's first author and year
agrees <- function(candidate, surname, year) {
  names <- coalesce(candidate$author_name %||% candidate$first_author, candidate$first_author)
  same_author <- vapply(seq_len(nrow(candidate)), function(i) {
    parts <- c(candidate$first_author[i],
               strsplit(coalesce(names[i], ""), "[ ,.\\-]+")[[1]])
    parts <- parts[!is.na(parts) & nchar(parts) > 1]
    length(parts) > 0 && any(same_surname(parts, rep(surname, length(parts))))
  }, logical(1))
  same_year <- is.na(year) |
    (!is.na(candidate$year) & abs(as.integer(candidate$year) - year) <= 1)
  same_author & same_year
}

# a citation's first author and year (citation_author_year())
citation_surname <- function(citation) {
  sub(" \\d{4}$", "", citation_author_year(redecode(citation)))
}
citation_year <- function(citation) {
  suppressWarnings(as.integer(str_extract(citation_author_year(redecode(citation)),
                                          "\\d{4}$")))
}

# the paper a found record describes, checked against the citation, and
# followed from a preprint to its published paper
as_result <- function(p, citation, method, note = NA_character_) {
  if (!is.na(p$doi[1]) && is_preprint(p$doi[1])) {
    later <- published_version(p[1, ])
    if (!is.null(later)) {
      p <- later
      method <- paste(method, "(preprint, then its published paper)")
    }
  }
  surname <- citation_surname(citation)
  year <- citation_year(citation)
  own_title <- citation_title(citation)
  same_title <- !is.na(own_title) && !is.na(p$title[1]) &&
    title_overlap(p$title[1], own_title, of = max) >= 0.8
  same_year <- !is.na(year) && !is.na(p$year[1]) && abs(as.integer(p$year[1]) - year) <= 1
  ok <- agrees(p[1, ], surname, year) ||
    (same_title && (same_year || agrees(p[1, ], surname, NA)))
  p[1, ] %>%
    select(doi, pmid, title, first_author, year, journal) %>%
    mutate(method = method,
           check = if_else(ok, "ok", "author or year differs"),
           note = note)
}

# 1. by an identifier in the citation or its link, or on the linked page; a
# paper that disagrees with the citation's first author or year (a wrong link)
# isn't taken
by_identifier <- function(citation, link) {
  result <- by_identifier_any(citation, link)
  if (!is.null(result) && result$check %in% "ok") result else NULL
}
by_identifier_any <- function(citation, link) {
  text <- paste(citation, coalesce(URLdecode(coalesce(link, "")), ""))
  doi <- find_doi(text)
  if (!is.na(doi)) {
    p <- from_crossref(doi)
    # Europe PMC returns other papers for a DOI it doesn't know
    if (is.null(p)) p <- filter(bind_rows(paper_row(), from_europepmc(paste0('DOI:"', doi, '"'))),
                                doi %in% .env$doi)
    if (nrow(p)) return(as_result(p, citation, "doi"))
  }
  pmid <- find_pmid(text)
  if (!is.na(pmid)) {
    p <- from_europepmc(paste0("EXT_ID:", pmid, " AND SRC:MED"))
    if (!is.null(p)) return(as_result(p, citation, "pubmed id"))
  }
  pmc <- find_pmc(text)
  if (!is.na(pmc)) {
    p <- from_europepmc(paste0("PMCID:", pmc))
    if (!is.null(p)) return(as_result(p, citation, "pmc id"))
  }
  if (!is.na(link) && grepl("^https?://", link) && !grepl("\\.pdf($|\\?)", link)) {
    page <- get_page(link)
    page_doi <- find_doi(str_match(page, paste0(
      "(?i)(?:citation_doi|dc\\.identifier|doi\\.org/)[^0-9]{0,40}",
      "(10\\.\\d{4,9}/[^\"'<>\\s]+)"))[, 2])
    if (!is.na(page_doi)) {
      p <- from_crossref(page_doi)
      if (!is.null(p)) return(as_result(p, citation, "doi on linked page"))
    }
    # a journal page without a DOI often gives the paper's title, authors and
    # date in its metadata
    meta <- function(name) {
      str_match(page, paste0('(?i)<meta[^>]+name="', name, '"[^>]+content="([^"]+)"'))[, 2]
    }
    page_title <- meta("citation_title")
    if (!is.na(page_title)) {
      p <- paper_row(title = page_title,
                     first_author = sub(",.*| .*", "", meta("citation_author")),
                     author_name = meta("citation_author"),
                     year = str_extract(coalesce(meta("citation_publication_date"),
                                                 meta("citation_date")), "(19|20)\\d{2}"),
                     journal = meta("citation_journal_title"))
      return(as_result(p, citation, "title on linked page"))
    }
  }
  NULL
}

# 2. by the citation's own title, searched in Crossref and Europe PMC: a paper
# by its first author with its title (four fifths of the shorter title's
# words, as some citations shorten titles). Candidates are ranked by the share
# of the longer title's words, so a paper whose title merely contains the
# citation's words ranks below one with the same title, then by year: the year
# only ranks, as a citation's year can be wrong (e.g. the IR Mapper extracts'
# Messenger 2017, whose title holds "(2012-2016)").
by_title <- function(citation) {
  title <- citation_title(citation)
  surname <- citation_surname(citation)
  year <- citation_year(citation)
  if (is.na(title) || is.na(surname)) return(NULL)
  candidates <- bind_rows(
    from_crossref_search(plain_text(citation), NA),
    from_europepmc(paste0("TITLE:(", paste(head(words_of(title)[[1]], 8), collapse = " "),
                          ') AND AUTH:"', surname, '"')))
  if (!nrow(candidates)) return(NULL)
  candidates <- candidates %>%
    filter(!is_notice(title, doi), agrees(., .env$surname, NA)) %>%
    mutate(overlap = vapply(title, function(t) title_overlap(t, .env$title), numeric(1)),
           same = vapply(title, function(t) title_overlap(t, .env$title, of = max), numeric(1))) %>%
    filter(overlap >= 0.8) %>%
    arrange(desc(same), abs(year - coalesce(.env$year, year)), is.na(doi))
  if (!nrow(candidates)) return(NULL)
  as_result(candidates, citation, "title search")
}

# 3. from the citation's own text, where it gives a title but no paper was
# found
by_own_text <- function(citation) {
  title <- citation_title(citation)
  if (is.na(title) || is.na(citation_author_year(citation))) return(NULL)
  paper_row(title = title, first_author = citation_surname(citation),
            year = citation_year(citation)) %>%
    as_result(citation, "citation text")
}

# 4. by first author and year, with words of the field, in Crossref and Europe
# PMC, preferring the citation's own year: the candidate whose title the
# citation's linked PDF or page holds, or else, of several, the one naming the
# countries and sites of the citation's records. A single candidate with
# neither is left for manual review (e.g. "Hunt 2013" finds only "Cerebral
# malaria: gamma-interferon redux").
by_author_year <- function(citation, countries, sites, text = NA_character_) {
  surname <- citation_surname(citation)
  year <- citation_year(citation)
  if (is.na(surname)) return(NULL)
  candidates <- bind_rows(
    from_europepmc(paste0('AUTH:"', surname, '" AND PUB_YEAR:[', year - 1, " TO ", year + 1,
                          "] AND (anopheles OR insecticide OR mosquito OR malaria)")),
    from_crossref_search(paste(surname, "insecticide resistance Anopheles",
                               gsub("\\+", " ", countries)), year))
  if (!nrow(candidates)) return(NULL)
  candidates <- candidates %>%
    filter(!is_notice(title, doi), agrees(., .env$surname, .env$year),
           grepl(field_words, tolower(paste(title, coalesce(abstract, ""))))) %>%
    arrange(is.na(doi)) %>%
    one_per_paper()
  if (any(candidates$year %in% year)) candidates <- filter(candidates, year %in% .env$year)
  if (nrow(candidates) == 0) return(NULL)
  # the one whose title the citation's linked PDF or page holds
  if (!is.na(text)) {
    in_text <- vapply(candidates$title, words_found, numeric(1), text = text) >= 0.8
    if (sum(in_text) == 1) {
      return(as_result(candidates[in_text, ], citation, "author and year search",
                       paste(nrow(candidates), "candidates; chosen by the title in its linked document")))
    }
  }
  places <- unique(tolower(c(strsplit(countries, "\\+")[[1]], strsplit(sites, "\\|")[[1]])))
  places <- places[nchar(places) > 3]
  score <- vapply(seq_len(nrow(candidates)), function(i) {
    text_i <- tolower(paste(candidates$title[i], coalesce(candidates$abstract[i], "")))
    sum(vapply(places, function(w) grepl(w, text_i, fixed = TRUE), logical(1)))
  }, numeric(1))
  best <- which(score == max(score))
  if (nrow(candidates) > 1 && length(best) == 1 && max(score) >= 1) {
    return(as_result(candidates[best, ], citation, "author and year search",
                     paste(nrow(candidates), "candidates; chosen by places named")))
  }
  attr(candidates, "ambiguous") <- TRUE
  candidates
}

# a link written one way, for comparing links: without its protocol, "www.",
# case, escapes and a trailing slash; NA for anything that isn't a web address
link_key <- function(link) {
  x <- tolower(URLdecode(coalesce(link, "")))
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- sub("^https?://(www\\.)?", "", x)
  x <- sub("/+$", "", x)
  if_else(grepl("^https?://", coalesce(link, ""), ignore.case = TRUE), x, NA_character_)
}

# the text of the first page of a linked PDF, line by line (cached), or NA
pdf_text <- function(link) {
  if (is.na(link) || !grepl("^https?://", link)) return(NA_character_)
  value <- pdf_text_cached(link)
  if (identical(attr(value, "failed"), TRUE) || !length(value)) NA_character_ else value
}
pdf_text_cached <- function(link) {
  cached(paste0("pdf lines:", link), function(url) {
    file <- tempfile(fileext = ".pdf")
    response <- tryCatch(
      request(link) %>% req_user_agent(user_agent) %>% req_timeout(60) %>%
        req_error(is_error = function(resp) FALSE) %>% req_perform(path = file),
      error = function(e) NULL)
    if (is.null(response)) return(failed)
    if (resp_status(response) != 200 || !grepl("pdf", resp_content_type(response))) return(NA_character_)
    text <- tryCatch(system2("pdftotext", c("-f", "1", "-l", "1", "-q", file, "-"), stdout = TRUE),
                     error = function(e) character())
    unlink(file)
    if (!length(text)) NA_character_ else paste(text, collapse = "\n")
  })
}

# a paper's details from the first page of its PDF: its title, the lines
# above its author list (the first line naming the cited first author among
# other names); its first author, the part of the first name in that list
# that matches the cited one, with the full list in the note; its year, the
# first printed on the page; and its journal, from the page's first line (e.g.
# "Journal of Entomology and Zoology Studies 2017; 5(1): 701-707")
from_first_page <- function(text, surname) {
  lines <- str_squish(strsplit(text, "\n")[[1]])
  authors <- which(grepl(name_key(surname), name_key(lines), fixed = TRUE) &
                     grepl(",| & | and ", lines))[1]
  if (is.na(authors) || authors < 2) return(NULL)
  start <- max(c(0, which(lines[seq_len(authors - 1)] == ""))) + 1
  title <- paste(lines[start:(authors - 1)], collapse = " ")
  if (start >= authors || !between(str_count(title, "\\S+"), 4, 40)) return(NULL)
  after <- lines[authors:length(lines)]
  end <- authors - 2 + coalesce(which(after == "" | grepl("^(\\d|\\*|abstract)", after,
                                                         ignore.case = TRUE))[1], 2L)
  names <- str_squish(gsub("\\d+(,\\d+)*\\*?", "", paste(lines[authors:end], collapse = " ")))
  first_name <- str_split(sub("(,| & | and ).*", "", names), " ")[[1]]
  header <- lines[lines != ""][1]
  paper_row(title = title,
            first_author = first_name[same_surname(first_name, rep(surname, length(first_name)))][1],
            author_name = paste(first_name, collapse = " "),
            year = str_extract(text, "(19|20)\\d{2}"),
            journal = str_squish(sub("\\s*\\(?((19|20)\\d{2}|\\d+\\s*\\(\\d+\\)).*", "", header))) %>%
    mutate(authors = names)
}

# share of a title's words found in a text
words_found <- function(title, text) {
  w <- unique(words_of(title)[[1]])
  if (!length(w) || is.na(text)) return(0)
  mean(w %in% words_of(text)[[1]])
}

# the results two citations' records share: the same insecticide, number
# tested and number dead, between 0% and 100% mortality, whatever the year
# (databases give one paper's tests different years). `results` are
# ";"-separated "insecticide year tested dead", as the prep lists them.
result_keys <- function(results) {
  r <- strsplit(coalesce(results, ""), ";")[[1]]
  parts <- str_match(r, "^(.*) (\\d{4}) (\\d+) (\\d+)$")
  keep <- !is.na(parts[, 1]) & parts[, 5] != "0" & parts[, 4] != parts[, 5]
  unique(paste(parts[keep, 2], parts[keep, 4], parts[keep, 5]))
}
shared_results <- function(own, others) {
  own <- result_keys(own)
  vapply(others, function(o) length(intersect(own, result_keys(o))),
         numeric(1), USE.NAMES = FALSE)
}

# the one candidate whose citations share two or more results with this
# citation (one, for a citation with at most two), or none
by_shared_results <- function(candidates, own) {
  n_own <- length(result_keys(own))
  shared <- shared_results(own, candidates$results)
  pick <- which(shared >= if (n_own <= 2) 1 else 2)
  if (length(pick) == 1) candidates[pick, ] else candidates[0, ]
}

# a paper's identifier: its DOI, else its PubMed ID, else its first author,
# year and the first five words of its title
paper_id <- function(doi, pmid, first_author, year, title) {
  ref <- paste(name_key(first_author), year,
               vapply(words_of(title), function(w) paste(head(w, 5), collapse = " "),
                      character(1)))
  case_when(!is.na(doi) ~ paste0("doi:", doi),
            !is.na(pmid) ~ paste0("pmid:", pmid),
            !is.na(title) ~ paste0("ref:", ref),
            .default = NA_character_)
}

# ---- run ----

if (sys.nframe() == 0) {
  seen <- read_csv("outputs/review/bioassay_citations.csv",
                   col_types = cols(.default = "c"), na = "") %>%
    filter(is_paper_citation(citation))
  table_file <- "review/bioassay_citations.csv"
  columns <- c("database", "citation", "link", "doi", "pmid", "title",
               "first_author", "year", "journal", "method", "check", "note",
               "paper")
  done <- if (file.exists(table_file)) {
    read_csv(table_file, col_types = cols(.default = "c"), na = "")
  } else {
    as_tibble(setNames(rep(list(character()), length(columns)), columns))
  }
  todo <- anti_join(seen, filter(done, !method %in% c("unresolved", "ambiguous",
                                                     "link conflict")),
                    by = c("database", "citation", "link"))
  message(nrow(todo), " citations to resolve")

  # passes 1 to 3: identifiers, titles, own text
  pass1 <- parallel::mclapply(seq_len(nrow(todo)), function(i) {
    r <- todo[i, ]
    tryCatch(by_identifier(r$citation, r$link) %||% by_title(r$citation) %||%
               by_own_text(r$citation),
             error = function(e) NULL)
  }, mc.cores = 4)
  resolved <- bind_rows(lapply(seq_along(pass1), function(i) {
    if (is.null(pass1[[i]])) return(NULL)
    bind_cols(select(todo[i, ], database, citation, link), pass1[[i]])
  }))
  known <- bind_rows(filter(done, !method %in% c("unresolved", "ambiguous", "link conflict")) %>%
                       mutate(year = as.integer(year)),
                     resolved)

  # a paper found only from a citation's own text takes the DOI or PubMed ID
  # of another citation of it: the same first author and year, and a title
  # sharing four fifths of its words
  identified <- filter(known, !is.na(doi) | !is.na(pmid))
  for (i in which(resolved$method == "citation text")) {
    match <- identified %>%
      filter(agrees(., resolved$first_author[i], resolved$year[i])) %>%
      mutate(overlap = vapply(title, function(t) title_overlap(t, resolved$title[i]), numeric(1))) %>%
      filter(overlap >= 0.8) %>%
      arrange(desc(overlap))
    method <- "citation text, matched to another citation's paper"
    if (!nrow(match)) {
      # a paper whose title is given in another language (e.g. French in one
      # database, PubMed's English in another): the only paper by this first
      # author and year, or the one sharing its sites
      match <- identified %>%
        filter(agrees(., resolved$first_author[i], resolved$year[i])) %>%
        mutate(id = paper_id(doi, pmid, first_author, year, title)) %>%
        left_join(select(seen, database, citation, link, results),
                  by = c("database", "citation", "link")) %>%
        group_by(id) %>%
        summarise(across(c(doi, pmid, title, first_author, year, journal), first),
                  results = paste(results, collapse = ";"), .groups = "drop")
      # it must share results with the other citation's records: the only
      # paper found by this author and year may not be this one (e.g. a paper
      # in a journal with no DOI)
      own <- seen %>%
        filter(database == resolved$database[i], citation == resolved$citation[i])
      match <- by_shared_results(match, own$results[1])
      method <- "citation text, the same author, year and results as another citation's paper"
    }
    if (nrow(match)) {
      resolved[i, c("doi", "pmid", "title", "first_author", "year", "journal")] <-
        match[1, c("doi", "pmid", "title", "first_author", "year", "journal")]
      resolved$method[i] <- method
    }
  }
  known <- bind_rows(filter(done, !method %in% c("unresolved", "ambiguous", "link conflict")) %>%
                       mutate(year = as.integer(year)),
                     resolved)

  # pass 4: an "Author et al. Year" citation is the paper another database's
  # citation was resolved to, with the same first author and year, where there
  # is one; else it is searched for by author and year
  pending <- todo[vapply(pass1, is.null, logical(1)), ]
  papers <- known %>%
    mutate(id = paper_id(doi, pmid, first_author, year, title)) %>%
    filter(!is.na(id)) %>%
    left_join(select(seen, database, citation, link, countries, sites, results),
              by = c("database", "citation", "link"))
  resolve_pending <- function(r) {
    surname <- citation_surname(r$citation)
    year <- citation_year(r$citation)
    local <- papers %>%
      filter(agrees(., .env$surname, .env$year)) %>%
      group_by(id) %>%
      summarise(across(c(doi, pmid, title, first_author, journal), first),
                year = first(year),
                places = paste(c(countries, sites), collapse = "|"),
                results = paste(results, collapse = ";"), .groups = "drop") %>%
      # one paper under two identifiers (e.g. two DOIs) counts once
      arrange(is.na(doi)) %>%
      one_per_paper()
    result <- NULL
    # the same link as a resolved citation by the same first author
    same_link <- papers %>%
      filter(!is.na(link_key(link)), link_key(link) == link_key(r$link),
             same_surname(first_author, rep(surname, n())) |
               same_surname(citation_surname(citation), rep(surname, n()))) %>%
      distinct(id, .keep_all = TRUE)
    if (nrow(same_link) == 1) {
      return(bind_cols(select(r, database, citation, link),
                       as_result(same_link, r$citation, "same link as another citation")))
    }
    # a linked PDF whose first page holds the title of a resolved paper by the
    # same first author, or of a paper Crossref finds from that page
    is_pdf <- grepl("\\.pdf($|\\?)", coalesce(r$link, ""), ignore.case = TRUE)
    text <- if (is_pdf) {
      pdf_text(r$link)
    } else if (grepl("^https?://", coalesce(r$link, ""))) {
      plain_text(get_page(r$link))
    } else {
      NA_character_
    }
    first_page <- if (is_pdf && !is.na(text)) from_first_page(text, surname)
    if (is_pdf && !is.na(text)) {
      by_author <- papers %>%
        filter(same_surname(first_author, rep(surname, n()))) %>%
        distinct(id, .keep_all = TRUE)
      searched <- from_crossref_search(if (is.null(first_page)) substr(text, 1, 400) else first_page$title, year)
      # matched on the PDF's own title where its first page gives one, else on
      # the words of its first page
      pool <- bind_rows(by_author, searched) %>%
        filter(!is.na(title), !is_notice(title, doi), agrees(., .env$surname, .env$year)) %>%
        mutate(found = if (is.null(first_page)) {
          vapply(title, words_found, numeric(1), text = text)
        } else {
          vapply(title, title_overlap, numeric(1), b = first_page$title, of = max)
        }) %>%
        filter(found >= 0.8) %>%
        arrange(desc(found), is.na(doi)) %>%
        one_per_paper()
      if (nrow(pool) == 1) {
        return(bind_cols(select(r, database, citation, link),
                         as_result(pool, r$citation, "title in linked PDF")))
      }
    }
    # a readable linked PDF rules out the papers whose titles it doesn't hold
    if (is_pdf && !is.na(text) && nrow(local)) {
      local <- filter(local, vapply(title, words_found, numeric(1), text = text) >= 0.8)
    }
    # the papers of the citation's own year, where there are any; the one whose
    # citations share results with this citation's records, or else, of
    # several, the one sharing its sites. A single paper sharing no results is
    # left for manual review: a first author and year alone don't identify a
    # paper.
    if (any(local$year %in% year)) local <- filter(local, year %in% .env$year)
    chosen <- by_shared_results(local, r$results)
    if (nrow(chosen) == 1) {
      result <- as_result(chosen, r$citation, "same paper as another database's citation",
                          paste(nrow(local), if (nrow(local) == 1) "paper" else "papers",
                                "by this author and year; chosen by shared results"))
    } else if (nrow(local) > 1) {
      own <- unique(name_key(c(strsplit(r$countries, "\\+")[[1]], strsplit(r$sites, "\\|")[[1]])))
      shared <- vapply(local$places, function(p) {
        length(intersect(own, name_key(strsplit(p, "\\||\\+")[[1]])))
      }, numeric(1))
      if (sum(shared == max(shared)) == 1 && max(shared) > 0) {
        result <- as_result(local[which.max(shared), ], r$citation,
                            "same paper as another database's citation",
                            paste(nrow(local), "papers by this author and year; chosen by shared sites"))
      }
    }
    # a linked PDF whose first page gives a paper by the citation's first
    # author and year, holding none of the titles of the papers known by them,
    # is a paper of its own, with the details its first page gives
    if (is.null(result) && !is.null(first_page) && agrees(first_page, surname, year)) {
      known_titles <- papers %>%
        filter(same_surname(first_author, rep(surname, n()))) %>%
        pull(title)
      if (!any(vapply(known_titles, title_overlap, numeric(1), b = first_page$title) >= 0.8)) {
        result <- as_result(first_page, r$citation, "its own paper, read from its linked PDF",
                            paste("authors:", first_page$authors))
      }
    }
    if (is.null(result)) {
      found <- by_author_year(r$citation, r$countries, r$sites, text)
      result <- if (is.null(found)) {
        as_result(paper_row(), r$citation, "unresolved", "no candidate found")
      } else if (isTRUE(attr(found, "ambiguous"))) {
        as_result(paper_row(), r$citation, "ambiguous",
                  paste(paste0(found$first_author, " ", found$year, ": ", found$title,
                               " [", coalesce(found$doi, found$pmid), "]"),
                        collapse = " || "))
      } else {
        found
      }
    }
    bind_cols(select(r, database, citation, link), result)
  }
  later <- bind_rows(parallel::mclapply(seq_len(nrow(pending)), function(i) {
    r <- pending[i, ]
    tryCatch(resolve_pending(r), error = function(e) {
      bind_cols(select(r, database, citation, link),
                as_result(paper_row(), r$citation, "unresolved",
                          paste("error:", conditionMessage(e))))
    })
  }, mc.cores = 4))

  # a citation whose link leads to a different paper than the one found another
  # way is a conflict, for manual review, with both papers
  found <- bind_rows(resolved, later)
  # nothing left to resolve: no rows, but the table's columns
  if (!nrow(found)) found <- done[0, ]
  conflicts <- parallel::mclapply(seq_len(nrow(found)), function(i) {
    f <- found[i, ]
    if (f$method %in% c("unresolved", "ambiguous") || grepl("^(doi|pubmed id|pmc id)", f$method)) return(NULL)
    p <- tryCatch(by_identifier_any(f$citation, f$link), error = function(e) NULL)
    if (is.null(p) || p$check %in% "ok" || is.na(p$title)) return(NULL)
    if (identical(coalesce(p$doi, p$pmid), coalesce(f$doi, f$pmid))) return(NULL)
    paste(paste0(c(f$first_author, p$first_author), " ", c(f$year, p$year), ": ",
                 c(f$title, p$title), " [", c(coalesce(f$doi, f$pmid), coalesce(p$doi, p$pmid)), "]"),
          collapse = " || ")
  }, mc.cores = 4)
  for (i in which(!vapply(conflicts, is.null, logical(1)))) {
    found[i, c("doi", "pmid", "title", "first_author", "year", "journal", "check")] <- NA
    found$method[i] <- "link conflict"
    found$note[i] <- conflicts[[i]]
  }

  table <- done %>%
    semi_join(seen, by = c("database", "citation", "link")) %>%
    anti_join(found, by = c("database", "citation", "link")) %>%
    bind_rows(mutate(found, across(everything(), as.character))) %>%
    # (a paper left unidentified by manual review keeps its "label:" id)
    mutate(paper = coalesce(paper_id(doi, pmid, first_author, year, title), paper))

  # one paper found under two identifiers (e.g. the Journal of Medical
  # Entomology's old and new DOIs): the same first author and year, and titles
  # sharing nine tenths of their words, take the identifier most citations have
  ids <- table %>%
    filter(!is.na(paper), !is.na(title)) %>%
    count(paper, first_author, year, title, sort = TRUE)
  alias <- setNames(ids$paper, ids$paper)
  for (i in seq_len(nrow(ids))) {
    if (alias[[ids$paper[i]]] != ids$paper[i]) next
    same <- which(seq_len(nrow(ids)) > i &
                    ids$year == ids$year[i] &
                    same_surname(ids$first_author, rep(ids$first_author[i], nrow(ids))) &
                    vapply(ids$title, title_overlap, numeric(1), b = ids$title[i]) >= 0.9)
    alias[ids$paper[same]] <- ids$paper[i]
  }
  table %>%
    mutate(paper = coalesce(unname(alias[paper]), paper)) %>%
    arrange(database, citation) %>%
    write_excel_csv(table_file, na = "")
}
