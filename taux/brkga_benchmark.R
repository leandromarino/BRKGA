https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745


Data = as.matrix(dist(Data_capitals)),
Fo = Fobj_TSP,
Dc = Decoder_TSP,
n = nrow(Data_capitals),
rc = 0.7,
pe = 0.2,
pm = 0.2,
p = 100,
ng = 2000,
ngw = 500,
MaxTime = 3600,
MAX = FALSE,
Exa1 = NULL,
Exa2 = NULL,







Data = as.matrix(dist(Data_capitals))
Fo = Fobj_TSP
Dc = Decoder_TSP
n = nrow(Data_capitals)
rc = 0.7
pe = 0.2
pm = 0.2 
p = 100
ng = 2000
ngw = 500
MaxTime = 3600
MAX = FALSE
Exa1 = NULL
Exa2 = NULL





brkga3 <-function(Data, Fo, Dc, rc = 0.7, pe = 0.2,
                  pm = 0.2, n, p = 100, ng = 2000, 
                  ngw = 500,
                  MaxTime = 3600, MAX = FALSE, 
                  Exa1 = NULL, Exa2 = NULL)
{
  cpu_time<-proc.time()
  time_iter<-0
  
  # starting parallel clustering
  cl <- parallel::makeCluster(parallel::detectCores())
  
  # RK Gen
  f<-popgen(n,p)
  
  # Decoder
  # g chave aleatória decodificada
  if(is.null(Exa1)==TRUE)
  {
    g<-t(parallel::parApply(cl, f,1,function(x) Dc(x)))
  } else {
    g<-t(parallel::parApply(cl, f,1,function(x) Dc(x,Exa1)))
  }
  if(nrow(g)==1) {g<-t(g)}
  
  # Objetive Function
  # ft valores da função objetivo nas solções de g
  if(is.null(Exa2)==TRUE){
    ft<-parallel::parApply(cl, g,1,function(x) Fo(Data,x))
  }  else {
    ft<-parallel::parApply(cl, g,1,function(x) Fo(Data,x,Exa2))}
  
  
  # finishing parallel clustering
  parallel::stopCluster(cl)
  
  
  #setting variables
  fbest <- ifelse(MAX == TRUE, -Inf,Inf)
  i <- 0
  pelite <- round(pe*p)
  pmutant <- round(pm*p)
  ngwb <- 0
  
  
  while( (i<ng) & (ngwb<=ngw) & (time_iter<MaxTime) )
  {
    i<-i+1
    ngwb<-ngwb+1
    pq <- order(ft, decreasing = MAX)
    f<-f[pq,] #Sorting by Fitness
    g<-g[pq,]
    if (is.matrix(g)==FALSE){g<-as.matrix(g)}
    fmin<-ft[pq[1]]
    ft<-ft[pq]
    if (MAX==FALSE)    {
      if (fmin<fbest)      {
        fbest<-fmin
        gbest<-g[1,]
        
        solution_best<-c(fbest,gbest)
        cat("Best Solution Generation ",i," = ",fbest,"\n")
        flush.console()
        ibest<-i
        ngwb<-0
      }
    } else
    {
      if (fmin>fbest)
      {fbest<-fmin
      gbest<-g[1,]
      solution_best<-c(fbest,gbest)
      cat("Best Solution Generation ",i," = ",fbest,"\n")
      flush.console()
      ibest<-i
      }
      
    }
    felite<-f[1:pelite,]
    fnonelite<-f[(pelite+1):p,] #Non-Elite
    fmutant<-popgen(n,pmutant)
    fnovos<-crossover(felite,fnonelite,rc,p,pelite,pmutant,n)
    fnew<-rbind(fmutant,fnovos)
    if(is.null(Exa1)==TRUE)
    { gnew<-t(apply(fnew,1,function(x) Dc(x)))}
    else {gnew<-t(apply(fnew,1,function(x) Dc(x,Exa1)))}
    if (nrow(gnew)==1) {gnew<-t(gnew)}
    
    if (ncol(g)==1)
    {g<-rbind(as.matrix(g[1:pelite,]),gnew)}
    else {g<-rbind(g[1:pelite,],gnew)}
    
    f<-rbind(felite,fnew)
    glk<-g[(pelite+1):p,]
    if (is.matrix(glk)==FALSE) {glk<-as.matrix(glk)}
    if (is.null(Exa2)==TRUE)
    {ftk<-apply(glk,1,function(x) Fo(Data,x))}
    else
    {ftk<-apply(glk,1,function(x) Fo(Data,x,Exa2))}
    
    ft<-c(ft[1:pelite],ftk)
    time_iter<-(proc.time()-cpu_time)[3]
  }
  

  
  cpu_time<-(proc.time()-cpu_time)[3]
  
  return(list(fbest=fbest,gbest=gbest,cpu_time=cpu_time,ibest=ibest))
}

















brkga2 <- compiler::cmpfun(brkga)


microbenchmark(
  BRKGA::brkga(Data = as.matrix(dist(Data_capitals)),
               Fo = Fobj_TSP,
               Dc = Decoder_TSP,
               n = nrow(Data_capitals),
               rc = 0.7,
               pe = 0.2,
               pm = 0.2,
               p = 20,
               ng = 50,
               ngw = 500,
               MaxTime = 3600,
               MAX = FALSE,
               Exa1 = NULL,
               Exa2 = NULL),
  brkga2(Data = as.matrix(dist(Data_capitals)),
         Fo = Fobj_TSP,
         Dc = Decoder_TSP,
         n = nrow(Data_capitals),
         rc = 0.7,
         pe = 0.2,
         pm = 0.2,
         p = 50,
         ng = 50,
         ngw = 500,
         MaxTime = 3600,
         MAX = FALSE,
         Exa1 = NULL,
         Exa2 = NULL),
  
  brkga3(Data = as.matrix(dist(BRKGA::Data_capitals)),
         Fo = BRKGA::Fobj_TSP,
         Dc = BRKGA::Decoder_TSP,
         n = nrow(BRKGA::Data_capitals),
         rc = 0.7,
         pe = 0.2,
         pm = 0.2,
         p = 50,
         ng = 50,
         ngw = 500,
         MaxTime = 3600,
         MAX = FALSE,
         Exa1 = NULL,
         Exa2 = NULL)
)




system.time(
  brkga(Data = as.matrix(dist(BRKGA::Data_capitals)),
        Fo = BRKGA::Fobj_TSP,
        Dc = BRKGA::Decoder_TSP,
        n = nrow(BRKGA::Data_capitals),
        rc = 0.7,
        pe = 0.2,
        pm = 0.2,
        p = 5000,
        ng = 50,
        ngw = 500,
        MaxTime = 3600,
        MAX = FALSE,
        Exa1 = NULL,
        Exa2 = NULL))

system.time(
  brkga3(Data = as.matrix(dist(BRKGA::Data_capitals)),
         Fo = BRKGA::Fobj_TSP,
         Dc = BRKGA::Decoder_TSP,
         n = nrow(BRKGA::Data_capitals),
         rc = 0.7,
         pe = 0.2,
         pm = 0.2,
         p = 5000,
         ng = 50,
         ngw = 500,
         MaxTime = 3600,
         MAX = FALSE,
         Exa1 = NULL,
         Exa2 = NULL)
)
