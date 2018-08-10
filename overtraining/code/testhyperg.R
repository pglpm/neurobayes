library('gsl')

for(n in seq(0.1, 110, 10)){
    for(k in seq(0.1, 2*n-0.1, 10)){
        x <- 2*n-k
        ifhyperg_2F1(k,2*n,k+1,-1)/k
    }}

resu1 <- unlist(sapply(seq(0.1, 110, 2),
       function(n){sapply(seq(0.1, 2*n-0.1, 2),
                          function(k){hyperg_2F1(k,2*n,k+1,-1)/k})}))

resu2 <- unlist(sapply(seq(0.1, 110, 2),
       function(n){sapply(seq(0.1, 2*n-0.1, 2),
                          function(k){x <- 2*n-k; hyperg_2F1(2*n,x,x+1,-1)/x})}))

resub2 <- unlist(sapply(seq(0.1, 110, 2.2),
       function(n){sapply(seq(0.1, 2*n-0.1, 2.1),
                          function(k){qqb2(k,n)})}))

write.table(resub2,'resub2.csv',sep=',',row.names=F,col.names=F,na='Null')

resuo <- unlist(sapply(seq(0.1, 110, 2.2),
       function(n){sapply(seq(0.1, 2*n-0.1, 2.1),
                          function(k){qqo(k,n)})}))

write.table(resuo,'resuo.csv',sep=',',row.names=F,col.names=F,na='Null')

lims <- function(xd,nt,ng){-(log(qqb2((nt+ng)*(2+xd),2*(ng+nt))/qqb2((nt)*(2+xd),2*(nt))))/(2*(ng))}



unlist(sapply(3:5,function(n){sapply(1:n,function(k){paste0(n,k)})}))

test <- (sapply(seq(0.1, 110, 10),
       function(n){sapply(seq(0.1, 2*n-0.1, 10),
                          function(k){qqb2(k,n)})}))

ttt <- unlist(sapply(seq(0.1, 110, 20),
       function(n){sapply(seq(0.1, 2*n-0.1, 20),
                          function(k){paste0(n,' ',k)})}))
