detectUps<-function(x, threshold=10000)
{  # check if x[i+1]-x[i] > threshold
  .Call("detect_ttl_ups_cwrap",
          x,
          length(x),
          threshold)
}
detectDowns<-function(x, threshold=10000)
{  # check if x[i+1]-x[i] < 0-threshold
  .Call("detect_ttl_downs_cwrap",
        x,
        length(x),
        threshold)
}
