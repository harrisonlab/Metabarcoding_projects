########################################################################
#
# This is a modified (for speed) version of the CRAN package cooccur
# All options remain the same
#
########################################################################

cooccur2 <- 
function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE,
    true_rand_classifier = 0.1, prob = "hyper", site_mask = NULL,
    only_effects = FALSE, eff_standard = TRUE, eff_matrix = FALSE)
{
    if (type == "spp_site") {
        spp_site_mat <- mat
    }
    if (type == "site_spp") {
        spp_site_mat <- t(mat)
    }
    if (spp_names == TRUE) {
        spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
    }
    if (!is.null(site_mask)) {
        if (nrow(site_mask) == nrow(spp_site_mat) & ncol(site_mask) ==
            ncol(spp_site_mat)) {
            N_matrix <- create.N.matrix(site_mask)
        }
        else {
            stop("Incorrect dimensions for site_mask, aborting.")
        }
    }
    else {
        site_mask <- matrix(data = 1, nrow = nrow(spp_site_mat),
            ncol = ncol(spp_site_mat))
        N_matrix <- matrix(data = ncol(spp_site_mat), nrow = nrow(spp_site_mat),
            ncol = nrow(spp_site_mat))
    }
    spp_site_mat[spp_site_mat > 0] <- 1
    tsites <- ncol(spp_site_mat)
    nspp <- nrow(spp_site_mat)
    spp_pairs <- choose(nspp, 2)
    incidence <- prob_occur <- obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs,
        ncol = 3)
    incidence <- prob_occur <- matrix(nrow = nrow(N_matrix),
        ncol = ncol(N_matrix))

### hacked bit - produces incidence and a matrix of the pair data for each species
    mat_mask <- as.matrix(mat*site_mask)
    incidence <- t(apply(mat_mask,1,function(m) site_mask%*%m))
    pairs <- t(apply(mat_mask,1,function(m) mat_mask%*%m))
    #	pairs[!upper.tri(pairs)]<-NA
    diag(incidence) <- NA
    prob_occur <- incidence/N_matrix
###    
    row <- 0

	sapply(seq(1,(nspp-1)),function(spp) {
	## This bit has been kept as is (just replaced the for loops with a couple of applies),
  ## but, along with the for loop below could be optimised
		sapply(seq((spp+1),nspp),function(spp_next) {
			row <<- row+1
			obs_cooccur[row, 1] <<- spp
			obs_cooccur[row, 2] <<- spp_next
			obs_cooccur[row, 3] <<- pairs[spp,spp_next]
			prob_cooccur[row, 1] <<- spp
			prob_cooccur[row, 2] <<- spp_next
			prob_cooccur[row, 3] <<- prob_occur[spp, spp_next] * prob_occur[spp_next, spp]
			exp_cooccur[row, 1] <<- spp
			exp_cooccur[row, 2] <<- spp_next
			exp_cooccur[row, 3] <<- prob_cooccur[row, 3] *N_matrix[spp, spp_next]		
		})
	})

 pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), style = 3)

    if (thresh == TRUE) {
        n_pairs <- nrow(prob_cooccur)
        prob_cooccur <- prob_cooccur[exp_cooccur[, 3] >= 1, ]
        obs_cooccur <- obs_cooccur[exp_cooccur[, 3] >= 1, ]
        exp_cooccur <- exp_cooccur[exp_cooccur[, 3] >= 1, ]
        n_omitted <- n_pairs - nrow(prob_cooccur)
        pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)),
            style = 3)
    }
    output <- data.frame(matrix(nrow = 0, ncol = 9))
    colnames(output) <- c("sp1", "sp2", "sp1_inc", "sp2_inc",
        "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt",
        "p_gt")
    for (row in 1:nrow(obs_cooccur)) {
        sp1 <- obs_cooccur[row, 1]
        sp2 <- obs_cooccur[row, 2]
        sp1_inc <- incidence[sp1, sp2]
        sp2_inc <- incidence[sp2, sp1]
        max_inc <- max(sp1_inc, sp2_inc)
        min_inc <- min(sp1_inc, sp2_inc)
        nsite <- N_matrix[sp1, sp2]
        psite <- as.numeric(nsite + 1)
        prob_share_site <- rep(x = 0, times = psite)
        if (prob == "hyper") {
            if (only_effects == FALSE) {
                all.probs <- phyper(0:min_inc, min_inc, nsite -
                  min_inc, max_inc)
                prob_share_site[1] <- all.probs[1]
                for (j in 2:length(all.probs)) {
                  prob_share_site[j] <- all.probs[j] - all.probs[j -
                    1]
                }
            }
            else {
                for (j in 0:nsite) {
                  if ((sp1_inc + sp2_inc) <= (nsite + j)) {
                    if (j <= min_inc) {
                      prob_share_site[(j + 1)] <- 1
                    }
                  }
                }
            }
        }
        if (prob == "comb") {
            if (only_effects == FALSE) {
                for (j in 0:nsite) {
                  if ((sp1_inc + sp2_inc) <= (nsite + j)) {
                    if (j <= min_inc) {
                      prob_share_site[(j + 1)] <- coprob(max_inc = max_inc,
                        j = j, min_inc = min_inc, nsite = nsite)
                    }
                  }
                }
            }
            else {
                for (j in 0:nsite) {
                  if ((sp1_inc + sp2_inc) <= (nsite + j)) {
                    if (j <= min_inc) {
                      prob_share_site[(j + 1)] <- 1
                    }
                  }
                }
            }
        }
        p_lt <- 0
        p_gt <- 0
        for (j in 0:nsite) {
            if (j <= obs_cooccur[row, 3]) {
                p_lt <- prob_share_site[(j + 1)] + p_lt
            }
            if (j >= obs_cooccur[row, 3]) {
                p_gt <- prob_share_site[(j + 1)] + p_gt
            }
            if (j == obs_cooccur[row, 3]) {
                p_exactly_obs <- prob_share_site[(j + 1)]
            }
        }
        p_lt <- round(p_lt, 5)
        p_gt <- round(p_gt, 5)
        p_exactly_obs <- round(p_exactly_obs, 5)
        prob_cooccur[row, 3] <- round(prob_cooccur[row, 3], 3)
        exp_cooccur[row, 3] <- round(exp_cooccur[row, 3], 1)
        output[row, ] <- c(sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row,
            3], prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt,
            p_gt)
        setTxtProgressBar(pb, nspp + row)
    }
    close(pb)
    if (spp_names == TRUE) {
        sp1_name <- merge(x = data.frame(order = 1:length(output$sp1),
            sp1 = output$sp1), y = spp_key, by.x = "sp1", by.y = "num",
            all.x = T, sort = FALSE)
        sp2_name <- merge(x = data.frame(order = 1:length(output$sp2),
            sp2 = output$sp2), y = spp_key, by.x = "sp2", by.y = "num",
            all.x = T, sort = FALSE)
        output$sp1_name <- sp1_name[with(sp1_name, order(order)),
            "spp"]
        output$sp2_name <- sp2_name[with(sp2_name, order(order)),
            "spp"]
    }
    true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >=
        0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <=
        (tsites * true_rand_classifier)), ]))
    output_list <- list(call = match.call(), results = output,
        positive = nrow(output[output$p_gt < 0.05, ]), negative = nrow(output[output$p_lt <
            0.05, ]), co_occurrences = (nrow(output[output$p_gt <
            0.05 | output$p_lt < 0.05, ])), pairs = nrow(output),
        random = true_rand, unclassifiable = nrow(output) - (true_rand +
            nrow(output[output$p_gt < 0.05, ]) + nrow(output[output$p_lt <
            0.05, ])), sites = N_matrix, species = nspp, percent_sig = (((nrow(output[output$p_gt <
            0.05 | output$p_lt < 0.05, ])))/(nrow(output))) *
            100, true_rand_classifier = true_rand_classifier)
    if (spp_names == TRUE) {
        output_list$spp_key <- spp_key
        output_list$spp.names = row.names(spp_site_mat)
    }
    else {
        output_list$spp.names = c(1:nrow(spp_site_mat))
    }
    if (thresh == TRUE) {
        output_list$omitted <- n_omitted
        output_list$pot_pairs <- n_pairs
    }
    class(output_list) <- "cooccur"
    if (only_effects == F) {
        output_list
    }
    else {
        effect.sizes(mod = output_list, standardized = eff_standard,
            matrix = eff_matrix)
    }
}
