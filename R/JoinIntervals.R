join.intervals<-function(s1,e1,s2,e2){
  ## check args
  if(length(s1)!=length(e1))
    stop("unequal length of first intervals")
  if(length(s2)!=length(e2))
    stop("unequal lenght of the second intervals")
  if(any(s1>e1))
    stop("problem with chronology of the first set of intervals")
  if(any(s2>e2))
    stop("problem with chronology of the first set of intervals")
  
  l1=length(s1)
  l2=length(s2)
  ## for results
  s3=integer()
  e3=integer()
  index=1
  for (i in 1:l1)
    for (j in 1:l2)
    {
      # check if there is overlap between the two intervals
      if ((s2[j]>=s1[i] & s2[j] < e1[i]) |
          (s1[i]>=s2[j] & s1[i] < e2[j]))
      {
        # find the beginning of overlap, largest start
        if (s1[i]>= s2[j])
          s3[index]=s1[i]
        else
          s3[index]=s2[j]
        # find the end of overlap, smallest end
        if (e1[i]<= e2[j])
          e3[index]=e1[i]
        else
          e3[index]=e2[j]
        index=index+1
      }
    }
  data.frame(s=s3,e=e3)
}
