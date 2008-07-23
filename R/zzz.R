.onLoad <- function(lib,pkgname){
	require(cluster) || stop("cluster package could not be found") 
	require(Biobase) || stop("Biobase package could not be found")
}

.onUnload <- function( libpath ){
	library.dynam.unload( "hopach", libpath )		  
}

.First.lib<-function(libname,pkgname){
	library.dynam("hopach",pkgname,libname)
   require(cluster) || stop("can't load without cluster package")
	require(Biobase) || stop("can't load without Biobase package")
}


