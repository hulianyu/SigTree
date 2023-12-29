# library(cubt) # The package 'cubt'
###############            
# Tree Growth #
############### 
library(partitions) # Load the partitions package
CUBT_Cat_growth <- function (B, minsplit = 20, mindev = 0.001, minsize = 10, lp = 7, 
                             method = "cubt", critopt = "entropy") 
{
  if (!is.matrix(B)) 
    stop("Data must be within a matrix")
  if (mindev > 1 | mindev < 0) 
    stop("mindev should be within [0,1]")
  Terms <- colnames(B)
  nbmod = length(unique(as.vector(B)))
  if (critopt == "entropy" & nbmod > 32) 
    stop("Too much different values in data")
  N = nrow(B)
  individus = 1:N
  who = list()
  who[["1"]] = individus
  decoupe = rep(TRUE, 2^(lp + 1) - 1)
  decoupe[(2^lp):(2^(lp + 1) - 1)] = FALSE
  prof = 0
  res = list()
  while (prof < lp + 1) {
    nnodes = 2^prof
    nodes = nnodes:(2 * nnodes - 1)
    for (node in nodes) {
      nodec = as.character(node)
      nodecL = as.character(2 * node)
      nodecR = as.character(2 * node + 1)
      Bact = B[who[[nodec]], ]
      nr = nrow(Bact)
      if (is.null(nr) & (length(Bact) > 0)) {
        nr = 1
        Bact = t(as.matrix(Bact))
      }
      devp = crit(Bact, N, critopt)
      if (!is.nan(devp)) 
        if (decoupe[node]) {
          term = FALSE
          if (nr > minsplit && devp >= 1e-10) {
            zzzz = dim(Bact)
            cat(node, "-", zzzz, "\n", file = "essai.txt", 
                append = T)
            aa = split.node(Bact, N, critopt)
            if (gregexpr("in", aa$splits.cutleft) != 
                -1) 
              left = Bact[, aa$coord] %in% aa$seuil
            else left = Bact[, aa$coord] < aa$seuil
            devl = crit(Bact[left, ], N, critopt)
            devr = crit(Bact[!left, ], N, critopt)
            delta = (devp - devl - devr)
            if (prof == 0) 
              delta0 = delta
            gain = delta/delta0
            if (gain > mindev & min(sum(left), sum(!left)) >= 
                minsize) {
              res[[nodec]] = aa
              res[[nodec]]$dev = devp
              who[[nodecL]] = who[[nodec]][left]
              who[[nodecR]] = who[[nodec]][!left]
            }
            else term = TRUE
          }
          else term = TRUE
          if (term) {
            who[[nodecL]] = NULL
            who[[nodecR]] = NULL
            decoupe[2 * node + 1] = FALSE
            decoupe[2 * node] = FALSE
            res[[nodec]] = list(var = "<leaf>", n = nr, 
                                dev = devp, yval = "", splits.cutleft = "", 
                                splits.cutright = "", coord = 0, seuil = 0)
          }
        }
      else {
        if (!is.null(who[[nodec]])) 
          res[[nodec]] = list(var = "<leaf>", n = length(who[[nodec]]), 
                              dev = devp, yval = "", splits.cutleft = "", 
                              splits.cutright = "", coord = 0, seuil = 0)
        decoupe[2 * node + 1] = FALSE
        decoupe[2 * node] = FALSE
        who[[nodecL]] = NULL
        who[[nodecR]] = NULL
      }
    }
    prof = prof + 1
  }
  out = list(res = res, who = who, decoupe = decoupe)
  out$frame = const.frame(out$res)
  out$frame$label = as.numeric(row.names(out$frame))
  out$call = match.call()
  out$terms = Terms
  class(out) = "cubt"
  out
}

# Constructs the frame for a cubt object.
const.frame <- function (obj) 
{
  monres = do.call(rbind, obj)
  nbl = nrow(monres)
  monres = as.data.frame(monres)[, -c(7, 8)]
  monres$dev = as.numeric(monres$dev)
  monres$var = as.factor(as.character(monres$var))
  messplits = cbind(as.character(monres$splits.cutleft), as.character(monres$splits.cutright))
  colnames(messplits) = c("cutleft", "cutright")
  monres = monres[, -c(5, 6)]
  monres$splits = messplits
  monres$n = as.numeric(monres$n)
  monres$yval = monres$n
  ord = as.character(ordonne(as.numeric(row.names(monres))))
  monres = monres[ord, ]
  monres
}

# Computes the criterion used for the tree construction.
crit <- function (B, N, critopt = "entropy") 
{
  if (!is.loaded("cubt")) 
    library.dynam("cubt", "cubt", lib.loc = .libPaths())
  n1 = as.integer(ncol(B))
  n2 = as.integer(nrow(B))
  if (is.null(nrow(B)) || nrow(B) == 1) 
    return(0)
  if (critopt == "anova") 
    n3 = 0
  else {
    n3 = length(unique(as.integer(as.matrix(B))))
    B = as.numeric(as.factor(B))
  }
  res1 = as.double(numeric(1))
  res = .C("critR", as.double(as.matrix(B)), n2, n1, as.integer(n3), 
           as.integer(N), res1)
  res[[6]]
}

# Sorts a cubt frame as needed for the plot.cubt function.
ordonne <- function (val) 
{
  i = 1
  res = NULL
  while (length(res) <= length(val)) {
    pres = !is.na(match(i, val))
    if (pres) {
      res = c(res, i)
      i = 2 * i
    }
    else {
      i = i/2
      if (i%%2 == 0) 
        i = i + 1
      else if (i%%2 == 1) {
        j = i
        while (j%%2 == 1) j = j%/%2
        i = j + 1
      }
    }
  }
  res[-length(res)]
}

# Computes the optimal split for a node using a criterion defined by crit.
split.node <- function (xx, N, crit0 = "entropy") 
{
  xx = as.matrix(xx)
  if (mode(xx) != "character") {
    if (!is.loaded("cubt")) 
      library.dynam("cubt", "cubt", lib.loc = .libPaths())
    n1 = as.integer(ncol(xx))
    n2 = as.integer(nrow(xx))
    if (crit0 == "anova") 
      n3 = 0
    else n3 = length(unique(as.vector(xx)))
    res1 = as.double(numeric(1))
    res2 = as.integer(numeric(1))
    res3 = as.double(numeric(1))
    res = .C("split", as.double(as.matrix(xx)), n2, n1, 
             as.integer(n3), res1, res2, res3)
    res1 = res[[5]]
    res2 = res[[6]]
    res3 = res[[7]]
    out = list(var = paste("x", res2, sep = ""), n = n2, 
               dev = res3, yval = "", splits.cutleft = paste("<", 
                                                             signif(res1, 4), sep = ""), splits.cutright = paste(">", 
                                                                                                                 signif(res1, 4), sep = ""), coord = res2, seuil = res1)
  }
  else {
    critmax = 0
    p = ncol(xx)
    c00 = crit(xx, N, crit0)
    for (j in 1:p) {
      Xj = xx[, j]
      mod = sort(unique(Xj))
      nbmod = length(mod)
      if (nbmod == 1) {
        next
      }
      else {
        Pj = setparts(restrictedparts(nbmod, 2))
        for (k in 2:ncol(Pj)) {
          mod.l = mod[which(Pj[, k] == 1)]
          left = Xj %in% mod.l
          data.l = xx[left, ]
          data.r = xx[!left, ]
          val0 = c00 - crit(data.l, N, crit0) - crit(data.r, 
                                                     N, crit0)
          if (val0 > critmax) {
            var = paste("x", j, sep = "")
            cut.l = paste(" in {", paste(mod.l, collapse = " "), 
                          "}", sep = "")
            cut.r = paste(" not in {", paste(mod.l, 
                                             collapse = " "), "}", sep = "")
            critmax = val0
            seuil = mod.l
          }
        }
      }
    }
    out = list(var = var, n = nrow(xx), dev = critmax, yval = "", 
               splits.cutleft = cut.l, splits.cutright = cut.r, 
               coord = as.numeric(substr(var, 2, 10)), seuil = seuil)
  }
  out
}


##############            
# Tree Prune #
##############
CUBT_Cat_prune <- function (trclust, Dt, mindist = NULL, qdist = 0.8, alph = 0.3, 
                            nleaves = NULL, dist = FALSE, ...) 
{
  trout = trclust
  struct = trout$res
  who = trclust$who
  frame = trclust$frame
  node <- row.names(frame)
  depth <- tree.depth(as.numeric(node))
  depth <- split(seq(as.numeric(node)), depth)
  lp = length(depth)
  vals = list()
  jj = 1
  size = length(leaf(trclust))
  if (size == 1) {
    warning("Root node tree, cannot prune it...")
    return(trclust)
  }
  if (qdist < 0 | qdist > 1) {
    warning("Wrong value for qdist, i'll take the median distance for mindist...")
    qdist = 0.5
  }
  for (prof in lp:2) {
    ress = numeric()
    nodes = node[depth[[prof]]]
    while (length(nodes) > 0) {
      ress = c(ress, distAB(Dt, who[[nodes[1]]], who[[nodes[2]]], 
                            alpha = alph, ...))
      nodes = nodes[-c(1, 2)]
    }
    vals[[prof]] = ress
  }
  vals = vals[-1]
  if (dist) 
    return(vals)
  else {
    arret = FALSE
    if (is.null(mindist)) 
      mindist = quantile(unlist(vals), qdist, na.rm = T)
    if (!is.null(mindist) | !is.null(nleaves)) {
      for (prof in lp:2) {
        nodes = node[depth[[prof]]]
        dd = vals[[prof - 1]]
        for (jj in 1:length(dd)) {
          elag = FALSE
          if (is.null(struct[[as.character(as.numeric(nodes[1])%/%2)]])) {
            struct[[as.character(nodes[1])]] = NULL
            who[[as.character(nodes[1])]] = NULL
            struct[[as.character(nodes[2])]] = NULL
            who[[as.character(nodes[2])]] = NULL
          }
          else if (dd[jj] < mindist) {
            feuilleg = struct[[as.character(nodes[1])]]$var == 
              "<leaf>"
            feuilled = struct[[as.character(nodes[2])]]$var == 
              "<leaf>"
            if (feuilleg & feuilled) {
              struct[[as.character(nodes[1])]] = NULL
              who[[as.character(nodes[1])]] = NULL
              struct[[as.character(nodes[2])]] = NULL
              who[[as.character(nodes[2])]] = NULL
              pere = struct[[as.character(as.numeric(nodes[1])%/%2)]]
              struct[[as.character(as.numeric(nodes[1])%/%2)]] = list(var = "<leaf>", 
                                                                      n = pere$n, dev = pere$dev, yval = "", 
                                                                      splits.cutleft = "", splits.cutright = "", 
                                                                      coord = 0, seuil = 0)
              elag = TRUE
            }
          }
          nodes = nodes[-c(1, 2)]
          if (!is.null(nleaves)) {
            if (size <= nleaves) {
              arret = TRUE
              break
            }
            if (elag) 
              size = size - 1
          }
        }
        if (arret) 
          break
      }
    }
    trout$res = struct
    trout$who = who
    trout$frame = const.frame(trout$res)
    trout$frame$label = as.numeric(row.names(trout$frame))
  }
  trout
}

# Computes the depth of a tree from the list of nodes numbers.
tree.depth <- function (nodes) 
{
  depth <- floor(log(nodes, base = 2) + 1e-07)
  as.vector(depth - min(depth))
}

# Retrieve leaves names from a cubt tree.
leaf <- function (arb) 
{
  as.numeric(row.names(arb$frame)[arb$frame$var == "<leaf>"])
}

# Computes distance between nodes as defined in the original reference.
distAB <- function (A, lft, rgt, alpha = 0.1, met = "im", p0 = 2, 
                    typ = 7, prec = 1e-10, ncores = 10) 
{
  ll = length(lft)
  lr = length(rgt)
  n = ll + lr
  if (met == "im") {
    data = A[c(lft, rgt), ]
    madst = im(as.matrix(t(data)))
    AB = madst[1:ll, (ll + 1):n]
    rm(madst, data)
  }
  else if (met == "hamming") {
    data = A[c(lft, rgt), ]
    madst = hamming.distance(as.matrix(data))
    AB = madst[1:ll, (ll + 1):n]
    rm(madst, data)
  }
  res = apply(AB, 1, min)
  qB = quantile(res, alpha, type = typ)
  dB = mean(res[(res - qB) <= prec])
  res = apply(AB, 2, min)
  qA = quantile(res, alpha, type = typ)
  dA = mean(res[(res - qA) <= prec])
  rm(AB)
  gc()
  max(dA, dB)
}

# Computes hamming distance between 2 matrices.
hamming.distance <- function (x, y) 
{
  z <- NULL
  if (is.vector(x) && is.vector(y)) {
    z <- sum(x != y)
  }
  else {
    z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    for (k in 1:(nrow(x) - 1)) {
      for (l in (k + 1):nrow(x)) {
        z[k, l] <- hamming.distance(x[k, ], x[l, ])
        z[l, k] <- z[k, l]
      }
    }
    dimnames(z) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
  }
  z
}

# Computes Mutual information between columns of a matrix.
im <- function (xx) 
{
  n1 = as.integer(ncol(xx))
  n2 = as.integer(nrow(xx))
  for (j in 1:n1) xx[, j] = as.numeric(as.factor(as.character(xx[, 
                                                                 j])))
  n3 = as.integer(length(unique(as.numeric(xx))))
  res = as.double(numeric(n1 * n1))
  res = .C("im", as.double(xx), n1, n2, n3, res)
  res = res[[5]]
  dim(res) = c(n1, n1)
  res
}


#############            
# Tree Join #
#############
CUBT_Cat_join <- function (obj, Don, pas = NULL, nclass = 4, verb = FALSE, crit0 = "entropy") 
{
  obj2 = obj
  feuilles = which(obj$frame$var == "<leaf>")
  lf = leaf(obj)
  ou = where(obj)
  N = nrow(Don)
  nbf = length(unique(ou))
  if (nbf == 1) {
    obj2$classe = 1
    obj2$classopt = 1
    obj2$frame$yval[feuilles] = 1
    return(obj2)
  }
  matdist = matrix(NA, nbf, nbf)
  classe = 1:nbf
  rownames(matdist) = colnames(matdist) = classe
  if (is.null(pas)) 
    pas = nbf
  if (nclass < nbf) {
    for (j in 1:(nbf - 1)) for (k in (j + 1):nbf) {
      cjk = crit(Don[which(ou == j | ou == k), ], N, critopt = crit0)
      cj = crit(Don[which(ou == j), ], N, critopt = crit0)
      ck = crit(Don[which(ou == k), ], N, critopt = crit0)
      matdist[j, k] = (cjk - cj - ck)/cjk
    }
    if (verb) 
      print(c("Number of leaves", nbf))
    if (verb) 
      print(c("Number of steps", pas))
    if (verb) 
      print(c("Number of classes", nclass))
    hierar = matrix(NA, nrow = pas, ncol = nbf)
    ncl = numeric(pas)
    n.pas = pas
    kk = 1
    groupes = as.list(1:nbf)
    names(groupes) = 1:nbf
    for (i in 1:(pas - 2)) {
      clstar = mat.min(matdist)
      cl1 = clstar[1]
      cl2 = clstar[2]
      clrest = rownames(matdist)
      etiq1 = clrest[cl1]
      etiq2 = clrest[cl2]
      cop1 = which(classe == etiq1)
      cop2 = which(classe == etiq2)
      newcl = paste("CL", kk, sep = "")
      if (verb) 
        print(c(etiq1, "+", etiq2, "=>", newcl))
      newgrp = c(groupes[[etiq1]], groupes[[etiq2]])
      groupes[[etiq1]] = groupes[[etiq2]] = NULL
      classe[cop1] = classe[cop2] = paste("CL", kk, sep = "")
      ncl[i] = nbcl = length(unique(classe))
      matdist = matdist[-clstar, -clstar]
      newdist = rep(NA, nbcl - 1)
      ouj = !is.na(match(ou, newgrp))
      cj = crit(Don[ouj, ], N, critopt = crit0)
      for (k in 1:(nbcl - 1)) {
        oujk = !is.na(match(ou, c(newgrp, groupes[[k]])))
        cjk = crit(Don[oujk, ], N, critopt = crit0)
        ouk = !is.na(match(ou, groupes[[k]]))
        ck = crit(Don[ouk, ], N, critopt = crit0)
        newdist[k] = (cjk - cj - ck)/cjk
      }
      matdist = rbind(cbind(matdist, newdist), NA)
      rownames(matdist)[nbcl] = colnames(matdist)[nbcl] = newcl
      groupes[[newcl]] = newgrp
      if (verb) 
        print(matdist)
      hierar[i, ] = classe
      kk = kk + 1
    }
    if (verb) 
      print("distance matrix after join")
    if (verb) 
      print(matdist)
    if (verb) 
      print("number of classes through join")
    if (verb) 
      print(ncl)
    if (verb) 
      print("hierar")
    if (verb) 
      print(hierar)
    laquelle = match(nclass, ncl)
    if (is.na(laquelle)) 
      laquelle = n.pas
    if (is.null(dim(hierar))) 
      classe = hierar
    else classe = hierar[laquelle, ]
    if (verb) 
      print("selected classes")
    if (verb) 
      print(classe)
    obj2$hierar = hierar
    obj2$nbclasses = ncl
    obj2$classopt = classe
    obj2$frame$yval[feuilles] = classe
    obj2$frame$label = obj2$frame$yval
  }
  else {
    print("Nothing to do")
    obj2 = obj
  }
  out = clean(obj2)
  out
}

# Computes the cluster of each observation used in the tree construction.
where <- function (arb) 
{
  leaves = row.names(arb$frame)[arb$frame$var == "<leaf>"]
  nbcl = length(leaves)
  who = arb$who[leaves]
  long = unlist(lapply(who, length))
  who = unlist(who)
  classe = rep(1:nbcl, long)
  classe[order(who)]
}

# Aggregates sibling leaves having the same label in a joined cubt tree.
clean <- function (ar) 
{
  fr = ar$frame
  leaf = fr$var == "<leaf>"
  labs = fr$label
  lableaf = unique(labs[leaf])
  nl = length(lableaf)
  newlab = paste("C", 1:nl, sep = "")
  for (i in 1:nl) {
    fr$label[labs == lableaf[i]] = newlab[i]
    fr$yval[labs == lableaf[i]] = newlab[i]
  }
  ar$frame = fr
  who = ar$who
  done = F
  while (!done) {
    labs = fr$label
    nds = as.numeric(rownames(fr)[fr$var == "<leaf>"])
    pleaf = nds[nds%%2 == 0]
    pleafn = match(pleaf, as.numeric(rownames(fr)))
    pleafnp1 = match(pleaf + 1, as.numeric(rownames(fr)))
    plabs = cbind(fr[as.character(pleaf), "label"], fr[as.character(pleaf + 
                                                                      1), "label"])
    del = NULL
    nf = 0
    for (i in 1:nrow(plabs)) {
      if ((fr[as.character(pleaf[i] + 1), "var"] == "<leaf>") & 
          (plabs[i, 1] == plabs[i, 2])) {
        nf = nf + 1
        del = c(del, pleafn[i], pleafnp1[i])
        pere = as.character(pleaf[i]%/%2)
        fr[pere, "label"] = plabs[i, 1]
        fr[pere, "yval"] = plabs[i, 1]
        fr[pere, "var"] = "<leaf>"
        fr[pere, "splits"] = ""
        who[[as.character(pleaf[i])]] = NULL
        who[[as.character(pleaf[i] + 1)]] = NULL
      }
    }
    if (nf == 0) 
      done = T
    else fr = fr[-del, ]
  }
  ar$frame = fr
  ar$who = who
  ar
}

# Computes the column number where the minimum is reached within each row of a matrix.
mat.min <- function (x) 
{
  d <- dim(x)
  A <- cbind(rep(1:d[1], d[2]), sort(rep(1:d[2], d[1])))
  out <- A[which.min(x), ]
  return(out)
}


#############            
# Plot Tree #
#############
Plot.CUBT <- function (x, y = NULL, type = c("proportional", "uniform"), 
                       ...) 
{
  if (inherits(x, "singlenode")) 
    stop("cannot plot singlenode tree")
  if (!inherits(x, "cubt")) 
    stop("not legitimate tree")
  type <- match.arg(type)
  uniform <- type == "uniform"
  dev <- dev.cur()
  if (dev == 1L) 
    dev <- 2L
  assign(paste(".Tree.unif", dev, sep = "."), uniform, envir = .GlobalEnv)
  invisible(treepl(treeco(x), node = as.numeric(row.names(x$frame)), 
                   ...))
}

# Preparing for the tree plot.
treepl <- function (xy, node, erase = FALSE, ...) 
{
  x <- xy$x
  y <- xy$y
  parent <- match((node%/%2), node)
  sibling <- match(ifelse(node%%2, node - 1L, node + 1L), 
                   node)
  xx <- rbind(x, x, x[sibling], x[sibling], NA)
  yy <- rbind(y, y[parent], y[parent], y[sibling], NA)
  if (any(erase)) {
    lines(c(xx[, erase]), c(yy[, erase]), col = par("bg"))
    return(x = x[!erase], y = y[!erase])
  }
  plot(range(x), range(y), type = "n", axes = FALSE, xlab = "", 
       ylab = "")
  text(x[1L], y[1L], "|", ...)
  lines(c(xx[, -1]), c(yy[, -1]), ...)
  list(x = x, y = y)
}

# Computes coordinates for tree plot.
treeco <- function (tree, uniform = paste(".Tree.unif", dev.cur(), sep = ".")) 
{
  frame <- tree$frame
  node <- as.numeric(row.names(frame))
  depth <- tree.depth(node)
  x <- -depth
  if (exists(uniform)) 
    uniform <- get(uniform)
  else uniform <- 0
  if (uniform) 
    y <- x
  else {
    y <- dev <- frame$dev
    depth <- split(seq(node), depth)
    parent <- match(node%/%2, node)
    sibling <- match(ifelse(node%%2, node - 1L, node + 1L), 
                     node)
    for (i in depth[-1L]) y[i] <- y[parent[i]] - dev[parent[i]] + 
      dev[i] + dev[sibling[i]]
  }
  depth <- -x
  leaves <- frame$var == "<leaf>"
  x[leaves] <- seq(sum(leaves))
  depth <- split(seq(node)[!leaves], depth[!leaves])
  left.child <- match(node * 2L, node)
  right.child <- match(node * 2 + 1L, node)
  for (i in rev(depth)) x[i] <- 0.5 * (x[left.child[i]] + 
                                         x[right.child[i]])
  list(x = x, y = y)
}

# Put labels on the tree plot.
Text.CUBT <- function (x, splits = TRUE, label = "yval", all = FALSE, pretty = NULL, 
                       digits = getOption("digits") - 3, adj = par("adj"), xpd = TRUE, 
                       ...) 
{
  oldxpd <- par(xpd = xpd)
  on.exit(par(oldxpd))
  if (inherits(x, "singlenode")) 
    stop("cannot plot singlenode tree")
  if (!inherits(x, "cubt")) 
    stop("not legitimate tree")
  frame <- x$frame
  column <- names(frame)
  if (!is.null(ylevels <- attr(x, "ylevels"))) 
    column <- c(column, ylevels)
  if (!is.null(label) && is.na(match(label, column))) 
    stop("label must be a column label of the frame component of the tree")
  charht <- par("cxy")[2L]
  if (!is.null(srt <- list(...)$srt) && srt == 90) {
    if (missing(adj)) 
      adj <- 0
    ladj <- 1 - adj
  }
  else ladj <- adj
  xy <- treeco(x)
  if (splits) {
    node <- as.numeric(row.names(frame))
    left.child <- match(2 * node, node)
    rows <- labels.cubt(x, pretty = pretty)[left.child]
    ind <- !is.na(rows)
    text(xy$x[ind], xy$y[ind] + 0.5 * charht, rows[ind], 
         adj = adj, ...)
  }
  if (!is.null(label)) {
    leaves <- if (all) 
      rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    if (label == "yval" & !is.null(ylevels)) 
      stat <- as.character(frame$yval[leaves])
    else if (!is.null(ylevels) && !is.na(lev <- match(label, 
                                                      ylevels))) 
      stat <- format(signif(frame$yprob[leaves, lev], 
                            digits = digits))
    else stat <- format(frame[leaves, label])
    if (!is.null(dim(stat)) && dim(stat)[2L] > 1) {
      if (length(dimnames(stat)[[2L]])) 
        stat[1, ] <- paste(sep = ":", dimnames(stat)[[2L]], 
                           stat[1, ])
      stat <- do.call("paste", c(list(sep = "\n"), split(stat, 
                                                         col(stat))))
    }
    text(xy$x[leaves], xy$y[leaves] - 0.5 * charht, labels = stat, 
         adj = ladj, ...)
  }
  invisible()
}

# Internal function for package CUBT used for labelling the tree plot.
labels.cubt <- function (object, pretty = TRUE, collapse = TRUE, ...) 
{
  if (!inherits(object, "cubt")) 
    stop("not legitimate tree")
  frame <- object$frame
  xlevels <- attr(object, "xlevels")
  var <- as.character(frame$var)
  splits <- matrix(sub("^>", " > ", sub("^<", " < ", frame$splits)), 
                   ,2)
  if (!is.null(pretty)) {
    if (pretty) 
      xlevels <- lapply(xlevels, abbreviate, minlength = pretty)
    for (i in grep("^:", splits[, 1])) for (j in 1L:2L) {
      sh <- splits[i, j]
      nc <- nchar(sh)
      sh <- substring(sh, 2L:nc, 2L:nc)
      xl <- xlevels[[var[i]]][match(sh, letters)]
      splits[i, j] <- paste(": ", paste(as.vector(xl), 
                                        collapse = ","), sep = "")
    }
  }
  if (!collapse) 
    return(array(paste(var, splits, sep = ""), dim(splits)))
  node <- as.numeric(row.names(frame))
  parent <- match((node%/%2), node)
  odd <- as.logical(node%%2)
  node[odd] <- paste(var[parent[odd]], splits[parent[odd], 
                                              2L], sep = "")
  node[!odd] <- paste(var[parent[!odd]], splits[parent[!odd], 
                                                1L], sep = "")
  node[1L] <- "root"
  node
}