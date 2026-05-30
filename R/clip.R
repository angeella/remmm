#' @title Cluster wise sign-flipping score test
#' @description This function implements a *cluster wise sign-flipping score-based test* for
#' one or multiple outcomes under within-cluster dependence (i.e., non-independent observations).
#' @usage clip(formula, data = NULL, cluster = NULL,
#' n_flips = 5000, alternative= "two.sided", seed = 1234)
#' @param formula One of:
#'   \enumerate{
#'     \item A \code{formula} object. For multivariate outcomes, \code{formula}
#'   may refer to a matrix response in \code{data} (as in \code{\link[stats]{lm}})
#'   or to a column in \code{data} (e.g. \code{cbind(y1,y2) ~ x})
#'     \item A list of fitted \code{glm} objects
#'     \item A list of \code{formula} objects
#'   }
#' @param data A \code{data.frame} containing the variables referenced in \code{formula} (used in cases 1 and 3).
#' @param n_flips Integer. Number of flips to generate. Default is 5000.
#' @param cluster Cluster identifier. \code{vector} or \code{formula} (e.g. \code{~ cluster}) naming the clustering column in \code{data}.
#' @param alternative Character string specifying the alternative hypothesis used to test the fixed effect
#' coefficients. One of \code{"two.sided"}, \code{"greater"}, \code{"less"}.
#' Default is \code{"two.sided"}.
#' @param seed Optional integer seed for reproducibility. Default \code{NULL}.
#' @param ... not implemented, yet.
#' @return A \code{remm} object, i.e., a list containing the following objects:
#' \describe{
#'   \item{Tspace}{\code{data.frame} where rows represents the sign-flipping transformed (plus the identity one) test and columns the variables.}
#'   \item{summary_table}{\code{data.frame} containing for each model the estimated parameter(s), score(s), std error(s), test(s), partial correlation(s) and p-value(s).}
#'   \item{call}{The matched call.}
#' }
#'
#' @seealso \code{\link[flipscores]{flipscores}}
#'
#' @importFrom flipscores flipscores
#' @importFrom flipscores combine_tests
#' @importFrom flipscores combine_contrasts
#' @import stats
#' @author Angela Andreella, Livio Finos
#' @examples
#' ##############################################
#' set.seed(1)
#' n=50
#' D=data.frame(matrix(rnorm(n*5),n,5))
#' Y=matrix(rnorm(n*2),n,2)
#' names(D)[4:5]=c("Y1","Y2")
#' cluster=rep(1:10,5)
#' res=clip(cbind(Y1,Y2)~X1+X2+X3,data=D,cluster=cluster)
#' summary(res)
#'
#' res=clip(list(Y1~X1+X2+X3,Y2~X1+X2+X3),data=D,cluster=cluster)
#' summary(res)
#'
#' mods=list(lm(Y1~X1+X2+X3,data=D),lm(Y2~X1+X2+X3,data=D))
#' res=clip(mods,data=D,cluster=cluster)
#' summary(res)
#' mods=list(list(formula=Y1~X1+X2+X3,data=D,cluster=cluster),
#'           list(formula=Y2~X1+X2+X3,data=D,cluster=cluster))
#' res=clip(mods)
#' summary(res)
#' summary(combine_tests(res))
#' summary(combine_tests(res, by="model"))
#' summary(combine_tests(res, by="coefficient"))
#' # to be implemented in flipscores: summary(p.adjust(res))
#' @export



clip <- function(formula,
                 data             = NULL,
                 cluster          = NULL,
                 n_flips          = 5000,
                 alternative      = "two.sided",
                 seed             = NULL,
                 tested_coeffs    = NULL,
                 flips            = NULL,
                 ...){

  original_call <- match.call()

  if(is.null(cluster)){
    if(is.list(formula)){
      cltrs=sapply(formula,function(md) {
        if(is(md$cluster,"formula"))
          cluster <- model.matrix(md$cluster, data = md$data)
        else
            cluster= md$cluster

        unique(cluster)})
      cluster_names=unique(as.vector(cltrs))
    }
  } else {
    if(is(cluster,"formula"))
      cluster <- model.matrix(cluster, data = data)
    cluster_names=unique(cluster)
  }
  n_obs=length(cluster_names)



  if(!is.null(seed)) set.seed(seed)

  if(is.null(flips)){
    flips=make_flips(n_obs=n_obs,n_flips=n_flips)
  }#else{
  #flips = flips
  #}


  ##############################################################
  # CASE 2: list of glm, lm or list objects
  ##############################################################
  if (is.list(formula) && all(sapply(formula, function(x)
    inherits(x, "glm")||inherits(x, "lm")))) {

    for(i in 1:length(formula))
      formula[[i]]$data=eval(formula[[i]]$call$data, parent.frame())

    out_list=lapply(formula, function(frm)
      .clip (formula(frm), frm$data, cluster,flips,alternative,
             cluster_names=cluster_names,tested_coeffs=tested_coeffs)
    )
    out=.make_output_from_list_Tspace_summary_table(out_list,original_call)
    return(out)

  }
  ## questa utile se vuoi mettere cluster specifici per ogni modello:
  if (is.list(formula) && all(sapply(formula, function(x)
    inherits(x, "list")))) {

    for(i in 1:length(formula))
      formula[[i]]$data=eval(formula[[i]]$data, parent.frame())

    out_list=lapply(formula, function(frm){
      if(is(frm$cluster,"formula"))
        cluster=model.matrix(frm$cluster, data = frm$data) else
          cluster= frm$cluster

      .clip(frm$formula, frm$data,
             cluster,flips,alternative,
             cluster_names=cluster_names,tested_coeffs=tested_coeffs)
    })
    out=.make_output_from_list_Tspace_summary_table(out_list,original_call)
    return(out)

  }

  ##############################################################
  # CASE 3: list of formulas
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "formula"))) {
    #message("flipscores: list of formulas detected -> converting to glms")

    out_list=lapply(formula, function(frm)
      .clip (frm, data, cluster,flips,alternative,
             cluster_names=cluster_names,tested_coeffs=tested_coeffs)
    )

    out=.make_output_from_list_Tspace_summary_table(out_list,original_call)
    return(out)
  }


  ##############################################################
  # CASE 1: standard formula -> original flipscores engine
  ##############################################################
  out=.clip(formula, data, cluster,flips,alternative,
            cluster_names=cluster_names,tested_coeffs=tested_coeffs)

  out$call <- original_call
  class(out) <- c("remmm", class(out))
  class(out) <- c("joint_flipscores", class(out))
  return(out)
}

###################################
.clip <- function(formula, data, cluster,flips,alternative,
                  cluster_names,tested_coeffs=NULL){
  D <- formula_to_matrices(formula, data = data)
  names_X=colnames(D$X)
  if(is.null(tested_coeffs)) tested_coeffs=names_X
  scores=lapply(tested_coeffs,function(i).get_scores(X=D$X[,i,drop=FALSE],
                                                     Y=D$Y,
                                                     Z=D$X[,setdiff(names_X,i),drop=FALSE],
                                                     cluster=cluster,
                                                     cluster_names=cluster_names))
  Tspace=lapply(scores,.flip_test,
                flips=flips)
  names(Tspace) <- names(scores) <- tested_coeffs


  summary_table=lapply(names(scores),function(i)
    cbind(coefficient=i,.make_summary_table(scores[[i]],Tspace[[i]],alternative)))

  Tspace=do.call(cbind,Tspace)
  summary_table=do.call(rbind,summary_table)
  summary_table=summary_table[,c(2,1,3:ncol(summary_table))]
  rownames(summary_table)=NULL
  list(Tspace=Tspace,summary_table=summary_table)
}



# for standardized (see in flipscores):
.score_std=function(flp,scores_objs) {
  # scr_eff # un vettore
  numerator=crossprod(flp,scores_objs$scores) #t(scr_eff)%*%flp
  if (all(sign(flp)==1)|(all(sign(flp)==-1))){
    denominator = 1
  } else {
    denominator = 1 - sum((colSums(scores_objs$vars_objs$A[flp==1,,drop=FALSE])
                           -colSums(scores_objs$vars_objs$A[flp==-1,,drop=FALSE]))^2)
  }
  as.vector(numerator/((denominator)**0.5))
}

######################
#X solo colonna
.get_scores<- function(X,Y,Z,cluster,cluster_names){
  Yr <- .get_IH(Z)%*%Y
  Q=qr.Q(qr(Z))
  Xr=crossprod(diag(nrow(Z))-tcrossprod(Q),X)
  m = sum(Xr^2)
  # we divide it by sqrt(m) which is the sd scaling factor of the observed test stat (i.e. effective and standardized have the same observed test stat)
  A=Xr[,]*Q/sqrt(m)
  scores=Xr[,]*Yr

  # raggruppa per cluster
  if(is.null(cluster)) cluster=1:nrow(scores)
  A <- rowsum(A, group = cluster)
  scores <- rowsum(scores, group = cluster)

  temp=fill_scores_by_cluster(list(scores=scores,A=A),cluster_names)
  list(scores=temp$scores, vars_objs=list(A=temp$A))#,Xr=Xr))
}

###################

.flip_test<- function(scores,
                      flips=NULL,
                      n_flips=NULL,
                      seed=NULL,
                      ...){


  if(!is.null(flips)){
    Tspace=sapply(1:nrow(flips),
                  function(i).score_std(flips[i,],scores))
    if(is.matrix(Tspace)) Tspace=t(Tspace)
    if(is.vector(Tspace)) Tspace=matrix(Tspace)


  } else {
    set.seed(seed)
    n_obs=nrow(scores)
    Tspace=rbind(.score_std(rep(1,n_obs),scores),t((replicate(n_flips-1,{
      .score_std(sample(c(-1,1),n_obs, replace = T),scores)
    }))))
  }

  return(Tspace)
}
