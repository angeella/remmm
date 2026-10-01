
mix_perm <- function(db,mod,rsubj_int,
                     rsubj_slop,
                     rstim_int, rstim_slop){
  mm=getME(mod,'X')

  y = db$y
  x = db$x
  subj = db$subj
  items = db$items

  mmis=model.matrix(~0+factor(subj))
  mmxs=model.matrix(~0+x:factor(subj))
  mmii=model.matrix(~0+factor(items))
  mmxi=model.matrix(~0+x:factor(items))


  y_=y-(mmxs%*%rsubj_slop+
          mmis%*%rsubj_int+
          mmii%*%rstim_int+
          mmxi%*%rstim_slop)

  rs=flip(y_,mm[,2])
  rs@res[,"p-value"]
}