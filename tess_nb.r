tess_nb<-function(coords,nb.type="queen",buffer_dist=0.1){
  # Voronoi tesselation and neighbour calculation
  # issue: st_voronoi changes the order of the geometries (points/polys in this case)
  # see https://github.com/r-spatial/sf/issues/824
  
  pp<-coords %>% as.data.frame %>% st_as_sf(coords = c(1,2)) %>% st_geometry
  dd<-st_geometry(st_multipoint(coords))
  hull<-st_geometry(st_buffer(dd,buffer_dist))
  v <- st_voronoi(do.call(c,dd), envelope = hull)
  inte<-st_collection_extract(v)
  inte<-st_intersection(st_cast(inte), hull)
  inte<-inte[unlist(st_intersects(pp,inte))] # this recovers the original order but is slow
  
  # neighbours calculation
  st_rook = function(a, b = a) st_relate(a, b, pattern = "F***1****")
  st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")
  if (nb.type=="queen"){
    queen_neighbours<-TRUE
  }else{
    queen_neighbours<-FALSE
  }
  if (queen_neighbours){
    nb<-st_queen(inte) 
  }else{
    nb<-st_rook(inte)
  }
  nb.matrix<-as.matrix(nb)+0 # added zero to convert to numeric
  A<-nb.matrix # adjacency matrix
 
  res<-list(A=A,tess=inte)
  return(res)
   
}
  