options(stringsAsFactors = FALSE)

tax = read.csv("ICTV.csv",
  as.is=TRUE,strip.white=TRUE)

n = c( "Species", paste("sp_",tax$Species,sep="") )
tnv = data.frame(Species=tax$Species)

TM = cbind( tnv, diag( length(tnv$Species) ) )
names(TM) = n

txc = function(prefix,data){
  data = unique(data)
  data = data[ data != "Unassigned" ]
  data = data[ data != "" ]
  data = data[ !is.na(data) ] 	
  paste(prefix,(data),sep="")
}

txc = c( txc("B_", tax$Baltimore),  txc("O_", tax$Order),
  txc("F_", tax$Family), txc("Sf_", tax$Subfamily), txc("G_", tax$Genus) )

##txc = chartr(" ","_",txc)
##txc = chartr("?","Q",txc)



txct = rep("Order",length(txc) )
txct = ifelse( substr(txc,1,1) =="F", "Family",txct)
txct = ifelse( substr(txc,1,1) =="S", "Subfamily",txct)
txct = ifelse( substr(txc,1,1) =="G", "Genus",txct)
txct = ifelse( substr(txc,1,1) =="B", "Baltimore",txct)


for(i in 1:length(txc) ){
  cat(". ")
  a = tax[[ txct[i] ]]
  b = substr( txc[i], 3, nchar(txc[i]) )
  d = 0 + (a==b)
  d = ifelse( is.na(d), 0, d)
  eval( parse(text=paste("TM$\"",txc[i],"\" = d",sep="")) )
}

write.table(TM,"FullPseudoTraitMatrixHP3virus.csv",sep=",",row.names=FALSE)


