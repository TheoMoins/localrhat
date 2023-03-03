# r_folder <- "../R/"
# source(paste(r_folder, "import/monitornew.R", sep=""))


trad_ess <- function(chains) {
    return(ess_rfun(chains))
}


local_ess <- function(x, chains) {
    res <- trad_ess(chains <= x)
    return(res)
}


all_local_ess <- function(chains, max_nb_points = 500) {
    grid <- grid_for_R(chains, max_nb_points)
    tab_rhat <- sapply(grid, function(x) local_ess(x, chains))
    return(tab_rhat[!is.na(tab_rhat)])
}

ess_infinity <- function(chains, dir = NULL, max_nb_points = 500) {
    val_ess = all_local_ess(chains, max_nb_points)
    return(max(val_ess[is.finite(val_ess)]))
}
