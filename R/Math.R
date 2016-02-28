make.pairs<-function(cl1="",cl2=NULL){
  if(is.null(cl2)){
    m<-combn(cl1,m=2)
    data.frame(Var1=m[1,],Var2=m[2,])
  }
  else{
    expand.grid(1:3,4:5)
  }
}

