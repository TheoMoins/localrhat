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
    tab_rhat <- c()
    grid <- grid_for_R(chains, max_nb_points)
    for (x in grid) {
        if (!is.na(local_ess(x, chains))){
            ess_x <- local_ess(x, chains)
        }
        tab_rhat <- c(tab_rhat, ess_x)
    }
    return(tab_rhat)
}

ess_infinity <- function(chains, dir = NULL, max_nb_points = 500) {
    val_ess = all_local_ess(chains, max_nb_points)
    return(max(val_ess[is.finite(val_ess)]))
}
