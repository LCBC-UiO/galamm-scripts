formula <- value_z ~ 0 + itemedu_father + itemedu_mother + itemedu_self +
  itemincome_father + itemincome_mother + itemincome_self + 
  itemgrouphippocampus:(siteousAvanto + siteousPrisma + siteousSkyra + icv_z + sexmale) + 
  s(age_z, by = itemgrouphippocampus, k = 15, bs = "cr")

random <- list(id = pdDiag(~ 0 + weight + itemgrouphippocampus))
weights <- varIdent(form =~ 1 | itemgroup)

modcall <- function(formula, random, dat, weights, loadings){
  dat$weight <- loadings[dat$item] + dat$itemgrouphippocampus * 
    loadings[["interaction"]] * dat$age_z
  gamm(
    formula = formula,
    random = random,
    data = dat,
    weights = weights,
    method = "ML"
  )
}

lambda_names <- c("edu_father", "edu_mother", "edu_self", 
                  "income_father", "income_mother", 
                  "income_self", "hippocampus", "interaction")

indsfun <- function(case){
  if(case == "parents_equal"){
    fixed_inds <- c(1L, 2L)
    equal_inds <- c(4L, 5L)
  } else if(case == "itemgroups_equal"){
    fixed_inds <- c(1L, 2L, 3L)
    equal_inds <- c(4L, 5L, 6L)
  } else if(case == "free_loadings"){
    fixed_inds <- 1L
    equal_inds <- integer()
  }
  list(fixed_inds = fixed_inds, equal_inds = equal_inds)
}

constrain_loadings <- function(loadings, equal_inds, zero_inds){
  if(length(equal_inds) > 1){
    loadings[equal_inds[-1]] <- loadings[[equal_inds[[1]]]]
  }
  if(length(zero_inds) > 0){
    loadings[zero_inds] <- 0
  }      
  loadings
}
