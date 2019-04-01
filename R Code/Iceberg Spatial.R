

library(SpRanger)


Names = intersect(FHN, names(MammalRanges))
N = length(Names)

t1 = Sys.time()

RangesList <- lapply(1:N, function(a){
  
  print(Names[a])
  
  lapply(1:N, function(b){
    
    if(a<b){
      
      if(RangeAdj[Names[a],Names[b]]>0){
        
        raster::intersect(MammalRanges[[Names[a]]], MammalRanges[[Names[b]]])
        
      }
      
    }
    
  })
  
})

Areas = sapply(unlist(RangesList), function(a) raster::freq(a))

t2 = Sys.time()


m1 = PairsWisely(MammalRanges)



t1 - t2

