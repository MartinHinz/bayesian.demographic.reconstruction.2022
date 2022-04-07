filter_by_element <- function(x, element) {
  x <- lapply(x, function(y) {
    select <- y$taxon.list$variable.element == element
    y$taxon.list <- y$taxon.list[select,]
    y$counts <- y$counts[,colnames(y$counts) %in% y$taxon.list$taxon.name]
    return(y)
  })
  class(x) <- append(class(x), "download_list")
  return(x)
}

filter_by_unit <- function(x, unit) {
  x <- lapply(x, function(y) {
    select <- y$taxon.list$variable.unit == unit
    y$taxon.list <- y$taxon.list[select,]
    y$counts <- y$counts[,colnames(y$counts) %in% y$taxon.list$taxon.name]
    return(y)
  })
  class(x) <- append(class(x), "download_list")
  return(x)
}

transform_to_percent <- function(x)
{
  x <- lapply(x, function(y) {
    y$counts <- analogue::tran(x = y$counts, method ="percent")
    return(y)
  })
  class(x) <- append(class(x), "download_list")
  return(x)
}

aggregate_by_joint_name <- function(x, common_names) {

  out <- pollen_data_comp %>% lapply(function(y) {
    y$taxon.list <- merge(y$taxon.list, common_names, by.x = "taxon.name", by.y = "name_orig")
    counts.new <- data.frame(t(y$counts))
    counts.new$name_joint <- NA
    counts.new[y$taxon.list$taxon.name,]$name_joint <- y$taxon.list$name_joint
    counts.collected <- counts.new  %>%
      as_tibble() %>%
      group_by(name_joint) %>%
      summarise(across(where(is.numeric), sum))
    counts.new <- data.frame(t(counts.collected[,2:ncol(counts.collected)]))

    colnames(counts.new) <- counts.collected$name_joint
    rownames(counts.new) <- rownames(y$counts)

    y$counts <- counts.new
    return(y)
  })

  class(out) <- append(class(out), "download_list")
  return(out)
}
