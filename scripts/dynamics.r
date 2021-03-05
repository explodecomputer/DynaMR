###########################

#solving an ode

 

require(deSolve)

 

 

expgrowth=function(Time,State,Params){

                with(as.list(c(State,Params)),{

                                                                dy=a*y

                                                                dx=b*y

                                                                res<-list(c(dy,dx))

                                                                return(res)

                                                })

}

 

#use documentation for plotting odes here:

#http://www.inside-r.org/packages/cran/deSolve/docs/hist.deSolve

 

out<-ode(func=expgrowth,y=c(y=5,x=0),parms=c(a=1.5,b=-1),times=seq(0, 5, by = .1))

 

###