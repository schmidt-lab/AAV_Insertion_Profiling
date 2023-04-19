#helper functions for Hoffmann_2023_main.R

#helper functions
'%!in%' = function(x,y)!('%in%'(x,y))

#write values to B-factor column of PDB file for plotting in pymol
wcpdb = function(
    file_in = 'aav_map_target.pdb',
    data_source,
    pos_col = 'pos_7kfr',
    data_col,
    na_val = 999999,
    file_out
)
{
  s = Rpdb::read.pdb(file_in)
  idx = match(s$atoms$resid, eval(as.symbol(data_source))[,pos_col])
  s$atoms$temp = as.numeric(eval(as.symbol(data_source))[idx,data_col])
  s$atoms$temp[is.na(s$atoms$temp)] = na_val
  Rpdb::write.pdb(s, file = file_out)
  
}

#pairwise Kolmogorov-Smirnov test
pairwise_ks_test = function(value, group, n_min = 50, warning = 0, alternative = "two.sided" ){
  
  lev = unique(group)
  
  lst = lapply( seq_along(lev), function(i) value[group == lev[i]] )
  names(lst)=lev
  
  if (sum(lengths(lst)< n_min)) {
    lst = lst [-which(lengths(lst)< n_min)]}
  
  f = function(x, y){ 
    w = getOption("warn") 
    options(warn = warning)  # ignore warnings 
    p = ks.test(x, y, alternative = alternative, exact = 
                  F)$p.value 
    options(warn = w) 
    return(p) 
  } 
  
  res = lapply(lst, function(x) lapply(lst, function(y) f(x, y))) 
  
  res = unlist(res)
  res = matrix(res, nrow = length(lst), ncol = length(lst), byrow = T)
  row.names(res) = colnames(res) = names(lst)
  cat("Pairwise Kolmogorov-Smirnov Test p-value Matrix","\n","\n")
  return(res)
}

#end helper functions
