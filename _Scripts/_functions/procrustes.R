#### ROTATION ####

#for loop version of the rotation function
rad_rot <- function(x,y,trace.x,trace.y,rot.plt,print.plt){
  trace_center <- c(mean(trace.x),mean(trace.y))
  cart_trace <- data.frame(x=trace.x,y=trace.y) %>%
    #once you check the implementaiton of this code this transformation to cart_trace
    #can be removed. I was using it as a sanity check
    mutate(x=x-trace_center[1],y=y-trace_center[2]) %>%
    mutate(
      #this could be written as a function cart2pol (cartesian to polar coordinates)
      #input would be two vectors x&y
      #output would be two vectors radius&radians
      radius=sqrt(x^2+y^2), #find radius
      radian = atan(y/x), #find relative radians
      radian = ifelse(x<0 & y<0, radian+pi,
                      ifelse(x<0 & y>0, radian+pi,
                             ifelse(x>0 & y<0, radian+2*pi,
                                    radian))), #convert to actual radian
      
      quad = ifelse(radian<pi/2,1,
                    ifelse(radian<pi&radian>pi/2,2,
                           ifelse(radian<3*pi/2&radian>pi,3,
                                  4))),# assign quadrants based on actual radian
      
      #this could be written as a function pol2cart (polar to cartesian coordinates)
      #input would be two vectors radius & actual radian
      #output would be two vectors x&y
      rel_radian = ifelse(quad==2 | quad==3,radian-(pi),
                          ifelse(quad==4,radian-(2*pi),
                                 radian)), #recalculate relative radians
      ret_x=radius*cos(rel_radian), #obtain x
      ret_y=radius*sin(rel_radian), #obtain y
      ret_x=ifelse(quad == 2 | quad==3,-ret_x,ret_x), #reasign direction based on quadrant
      ret_y=ifelse(quad == 2 | quad==3,-ret_y,ret_y), #reasign direction based on quadrant
      ret_x=(ret_x)+trace_center[1], #transpose x back to original position
      ret_y=(ret_y)+trace_center[2] #transpose y back to original position
    )
  
  #range of rotation is 2pi radians or one full circle by pi/180 radians (one degree)
  rotation_range <- seq(from = 0, to = 2*pi, by = pi/180)
  
  #preallocate output (col1 is rotation factor, col 2 is error)
  res <- matrix(0, ncol=3,nrow=length(rotation_range))
  for (i in 1:length(rotation_range)) {
    #rotate the data and convert back to cartesian from polar
    rot_trace <- cart_trace %>% # call output from before
      mutate(rot_radian=radian+rotation_range[i], #rotate shape
        rot_radian = ifelse(rot_radian>(2*pi),rot_radian-(2*pi),rot_radian),
        
        new_quad = ifelse(radian<pi/2,1,
                   ifelse(radian<pi&radian>pi/2,2,
                   ifelse(radian<3*pi/2&radian>pi,3,
                   4))),# assign quadrants based on actual radian
       
        #this could be written as a function pol2cart (polar to cartesian coordinates)
        #input would be two vectors radius & actual radian
        #output would be two vectors x&y
        rel_radian = ifelse(new_quad==2,rot_radian-(pi),
                     ifelse(new_quad==3,rot_radian-(pi),
                     ifelse(new_quad==4,rot_radian-(2*pi),
                     rot_radian))), #recalculate relative radians
        ret_x=radius*cos(rel_radian), #obtain x
        ret_y=radius*sin(rel_radian), #obtain y
        ret_x=ifelse(quad == 2 | quad==3,-ret_x,ret_x), #reasign direction based on quadrant
        ret_y=ifelse(quad == 2 | quad==3,-ret_y,ret_y), #reasign direction based on quadrant
        ret_x=(ret_x)+trace_center[1], #transpose x back to original position
        ret_y=(ret_y)+trace_center[2] #transpose y back to original position
      ) %>%
      select(c(1:4,rot_radian,quad,new_quad,6:8))
    
    #plot rotation
    if(rot.plt==TRUE){
      print(plot_tracing_err(trace.x,trace.y,
                             rot_trace$ret_x,rot_trace$ret_y)+
              ggtitle(paste0(rotation_range[i]*(180/pi)," degree(s) of rotation.")))
      invisible(readline(prompt="Press [enter] to continue"))
    }
    res[i,1]=rotation_range[i]
    res[i,2]=ifelse(min(rot_trace$ret_x)< 0 |
                    max(rot_trace$ret_x)>screen_res[1] |
                    min(rot_trace$ret_y)< 0 |
                    max(rot_trace$ret_y)>screen_res[2], 0, 1) 
    #the line above flags shapes that have rotated off the screen
    res[i,3]=mean(line_length(x, y, rot_trace$ret_x,rot_trace$ret_y))
  }
  plaus <- res[res[,2]==1,]
  min <- which.min(plaus[,3])
  best_rot <- plaus[min,1]
  
  best_rot_trace <- cart_trace %>% # call output from before
    mutate(rot_radian=radian+best_rot, #rotate shape
      rot_radian = ifelse(rot_radian>(2*pi),rot_radian-(2*pi),rot_radian),
           
      new_quad = ifelse(radian<pi/2,1,
                 ifelse(radian<pi&radian>pi/2,2,
                 ifelse(radian<3*pi/2&radian>pi,3,
                 4))),# assign quadrants based on actual radian
           
      #this could be written as a function pol2cart (polar to cartesian coordinates)
      #input would be two vectors radius & actual radian
      #output would be two vectors x&y
      rel_radian = ifelse(new_quad==2,rot_radian-(pi),
                   ifelse(new_quad==3,rot_radian-(pi),
                   ifelse(new_quad==4,rot_radian-(2*pi),
                   rot_radian))), #recalculate relative radians
      ret_x=radius*cos(rel_radian), #obtain x
      ret_y=radius*sin(rel_radian), #obtain y
      ret_x=ifelse(quad == 2 | quad==3,-ret_x,ret_x), #reasign direction based on quadrant
      ret_y=ifelse(quad == 2 | quad==3,-ret_y,ret_y), #reasign direction based on quadrant
      ret_x=(ret_x)+trace_center[1], #transpose x back to original position
      ret_y=(ret_y)+trace_center[2] #transpose y back to original position
    ) %>%
    select(c(1:4,rot_radian,quad,new_quad,6:8))

  
  #plot best result
  if(print.plt==TRUE){
    plot_tracing_err(x,y, best_rot_trace$ret_x,best_rot_trace$ret_y)+
      labs(title = paste0("The attained error post rotation. 
                          Rotated =", as.character(plaus[min,1]*(180/pi)), " degree(s)",
                          " Error = ", as.character(plaus[min,3]),"."))
    
    #plot unrotated for comparison
    #plot_tracing_err(x,y,trace.x,trace.y)+
    # labs(title = paste0("The attained error post rotation. 
    #                     Rotated = 0 degree(s)",
    #                     " Error = ", as.character(plaus[1,3]),"."))
  }
  return(c(plaus[min,1]*(180/pi),plaus[min,3]))
}

#this implements the logic that we discussed and only computes a subset of the full space
# of scaling my only worry is that im not immediately sure if this function will find the
#global minimum of the full space. It might get stuck on a local minmum as the logic only
#compares the current vs the previous trial.
part_rad_rot <- function(x,y,trace.x,trace.y,rot.plt,print.plt){
  trace_center <- c(mean(trace.x),mean(trace.y))
  cart_trace <- data.frame(x=trace.x,y=trace.y) %>%
    #once you check the implementaiton of this code this transformation to cart_trace
    #can be removed. I was using it as a sanity check
    mutate(x=x-trace_center[1],y=y-trace_center[2]) %>%
    mutate(
      #this could be written as a function cart2pol (cartesian to polar coordinates)
      #input would be two vectors x&y
      #output would be two vectors radius&radians
      radius=sqrt(x^2+y^2), #find radius
      radian = atan(y/x), #find relative radians
      radian = ifelse(x<0 & y<0, radian+pi,
                      ifelse(x<0 & y>0, radian+pi,
                             ifelse(x>0 & y<0, radian+2*pi,
                                    radian))), #convert to actual radian
      
      quad = ifelse(radian<pi/2,1,
                    ifelse(radian<pi&radian>pi/2,2,
                           ifelse(radian<3*pi/2&radian>pi,3,
                                  4))),# assign quadrants based on actual radian
      
      #this could be written as a function pol2cart (polar to cartesian coordinates)
      #input would be two vectors radius & actual radian
      #output would be two vectors x&y
      rel_radian = ifelse(quad==2 | quad==3,radian-(pi),
                          ifelse(quad==4,radian-(2*pi),
                                 radian)), #recalculate relative radians
      ret_x=radius*cos(rel_radian), #obtain x
      ret_y=radius*sin(rel_radian), #obtain y
      ret_x=ifelse(quad == 2 | quad==3,-ret_x,ret_x), #reasign direction based on quadrant
      ret_y=ifelse(quad == 2 | quad==3,-ret_y,ret_y), #reasign direction based on quadrant
      ret_x=(ret_x)+trace_center[1], #transpose x back to original position
      ret_y=(ret_y)+trace_center[2] #transpose y back to original position
    )
  
  #range of rotation is 2pi radians or one full circle by pi/180 radians (one degree)
  rotation_range <- seq(from = 0, to = 2*pi, by = pi/180)
  
  #preallocate output (col1 is rotation factor, col 2 is error)
  res <- matrix(0, ncol=3,nrow=length(rotation_range))
  reverse=0
  minimum=0
  i=1
  while (minimum==0){
    if (reverse==0){
      this_rot=rotation_range[i]
    }
    else {
      this_rot=max(rotation_range[!rotation_range %in% res[,1]])
    }
    
    rot_trace <- cart_trace %>% # call output from before
      mutate(rot_radian=radian+this_rot, #rotate shape
             rot_radian = ifelse(rot_radian>(2*pi),rot_radian-(2*pi),rot_radian),
             
             new_quad = ifelse(radian<pi/2,1,
                               ifelse(radian<pi&radian>pi/2,2,
                                      ifelse(radian<3*pi/2&radian>pi,3,
                                             4))),# assign quadrants based on actual radian
             
             #this could be written as a function pol2cart (polar to cartesian coordinates)
             #input would be two vectors radius & actual radian
             #output would be two vectors x&y
             rel_radian = ifelse(new_quad==2,rot_radian-(pi),
                                 ifelse(new_quad==3,rot_radian-(pi),
                                        ifelse(new_quad==4,rot_radian-(2*pi),
                                               rot_radian))), #recalculate relative radians
             ret_x=radius*cos(rel_radian), #obtain x
             ret_y=radius*sin(rel_radian), #obtain y
             ret_x=ifelse(quad == 2 | quad==3,-ret_x,ret_x), #reasign direction based on quadrant
             ret_y=ifelse(quad == 2 | quad==3,-ret_y,ret_y), #reasign direction based on quadrant
             ret_x=(ret_x)+trace_center[1], #transpose x back to original position
             ret_y=(ret_y)+trace_center[2] #transpose y back to original position
      ) %>%
      select(c(1:4,rot_radian,quad,new_quad,6:8))
    
    #plot rotation
    if(rot.plt==TRUE){
      print(plot_tracing_err(trace.x,trace.y,
                             rot_trace$ret_x,rot_trace$ret_y)+
              ggtitle(paste0(this_rot*(180/pi)," degree(s) of rotation.")))
      invisible(readline(prompt="Press [enter] to continue"))
    }
    res[i,1]=this_rot
    res[i,2]=ifelse(min(rot_trace$ret_x)< 0 |
                      max(rot_trace$ret_x)>screen_res[1] |
                      min(rot_trace$ret_y)< 0 |
                      max(rot_trace$ret_y)>screen_res[2], 0, 1)
    #the line above flags shapes that have rotated off the screen
    res[i,3]=mean(line_length(x, y, rot_trace$ret_x,rot_trace$ret_y))
    
    if (i==1){
    }
    else if (i==2 & res[i,3]>res[i-1,3]) {
      reverse=1
    }
    else if (res[i,3]>res[i-1,3]){
      minimum=1
    }
    i=i+1
  }
  plaus <- res[res[,2]==1,]
  min <- ifelse(length(plaus)==3,1,which.min(plaus[,3]))
  best_rot <- ifelse(length(plaus)==3,plaus[1],plaus[min,1])
  best_rot_trace <- cart_trace %>% # call output from before
    mutate(rot_radian=radian+best_rot, #rotate shape
           rot_radian = ifelse(rot_radian>(2*pi),rot_radian-(2*pi),rot_radian),
           
           new_quad = ifelse(radian<pi/2,1,
                             ifelse(radian<pi&radian>pi/2,2,
                                    ifelse(radian<3*pi/2&radian>pi,3,
                                           4))),# assign quadrants based on actual radian
           
           #this could be written as a function pol2cart (polar to cartesian coordinates)
           #input would be two vectors radius & actual radian
           #output would be two vectors x&y
           rel_radian = ifelse(new_quad==2,rot_radian-(pi),
                               ifelse(new_quad==3,rot_radian-(pi),
                                      ifelse(new_quad==4,rot_radian-(2*pi),
                                             rot_radian))), #recalculate relative radians
           ret_x=radius*cos(rel_radian), #obtain x
           ret_y=radius*sin(rel_radian), #obtain y
           ret_x=ifelse(quad == 2 | quad==3,-ret_x,ret_x), #reasign direction based on quadrant
           ret_y=ifelse(quad == 2 | quad==3,-ret_y,ret_y), #reasign direction based on quadrant
           ret_x=(ret_x)+trace_center[1], #transpose x back to original position
           ret_y=(ret_y)+trace_center[2] #transpose y back to original position
    ) %>%
    select(c(1:4,rot_radian,quad,new_quad,6:8))
  
  #plot best result
  if(print.plt==TRUE){
    plot_tracing_err(x,y, best_rot_trace$ret_x,best_rot_trace$ret_y)+
      labs(title = paste0("The attained error post rotation. 
                          Rotated =", as.character(plaus[min,1]*(180/pi)), " degree(s)",
                          " Error = ", as.character(plaus[min,3]),"."))
    
    #plot unrotated for comparison
    #plot_tracing_err(x,y,trace.x,trace.y)+
    # labs(title = paste0("The attained error post rotation. 
    #                     Rotated = 0 degree(s)",
    #                     " Error = ", as.character(plaus[1,3]),"."))
  }
  return(data.frame(x,y,best_rot_trace$ret_x,best_rot_trace$ret_y))
}

#### SCALING ####

#vectorized version of the scaling function that removes the implementaiton of line length
fullVec.scale <- function(x,y,trace.x,trace.y,scale_res,print.plt){
  x_scale <- screen_res[1]/(range(trace.x)[2]-range(trace.x)[1])
  y_scale <- screen_res[2]/(range(trace.y)[2]-range(trace.y)[1])
  scale_range <- seq(from = 0.01, to = round(min(x_scale,y_scale),2), by = scale_res)
  trace.xad = trace.x %*% t(scale_range)
  trace.yad = trace.y %*% t(scale_range)
  xad.err = trace.xad-x
  yad.err = trace.yad-y
  err = sqrt(xad.err**2+yad.err**2)
  scale_err = colMeans(err)
  min <- which.min(scale_err)
  if(print.plt==TRUE){
    thisy <- (trace.yad[,min])
    thisx <- (trace.xad[,min])
    plot_tracing_err(x,y,thisx,thisy)+
      labs(title = paste0("The attained error post scaling with vectorized. 
                          Scale factor =", as.character(min*0.01),
                          " Error = ", as.character(scale_err[min]),".")
      )
  }
  return(c(min*0.01, scale_err[min]))
}

#for loop version of the scaling function
fullvec.scale <- function(x,y,trace.x,trace.y,scale_res,print.plt){
  x_scale <- screen_res[1]/(range(trace.x)[2]-range(trace.x)[1])
  y_scale <- screen_res[2]/(range(trace.y)[2]-range(trace.y)[1])
  scale_range <- seq(from = 0.01, to = round(min(x_scale,y_scale),2), by = scale_res)
  res <- matrix(0, ncol=2,nrow=length(scale_range))
  for (i in 1:length(scale_range)) {
    trace.xad <- (trace.x*scale_range[i])
    trace.yad <- (trace.y*scale_range[i])
    res[i,1]=scale_range[i]
    res[i,2]=mean(line_length(x, y, trace.xad,trace.yad))
  }
  min <- which.min(res[,2])
  if(print.plt==TRUE){
    thisy <- (trace.y*res[min,1])
    thisx <- (trace.x*res[min,1])
    plot_tracing_err(x,y,thisx,thisy)+
      labs(title = paste0("The attained error post scaling. 
                          Scale factor =", as.character(res[min,1]),
                          " Error = ", as.character(res[min,2]),".")
      )
  }
  return(res[min,])
}

#this implements the logic that we discussed and only computes a subset of the full space
# of scaling
partvec.scale <- function(x,y,trace.x,trace.y,scale_res,print.plt){
  x_scale <- screen_res[1]/(range(trace.x)[2]-range(trace.x)[1])
  y_scale <- screen_res[2]/(range(trace.y)[2]-range(trace.y)[1])
  scale_range <- seq(from = 0.01, to = round(min(x_scale,y_scale),2), by = scale_res)
  res <- matrix(NA, ncol=2,nrow=length(scale_range))
  reverse=0
  minimum=0
  i=1
  while (minimum==0){
    if (reverse==0){
      this_scale=scale_range[which(scale_range==1)+(i-1)]
      if (this_scale==max(scale_range)){
        reverse=1
      }
    }
    else {
      this_scale=max(scale_range[which(scale_range<this_scale)])
    }
    trace.xad <- (trace.x*this_scale)
    trace.yad <- (trace.y*this_scale)
    res[i,1]=this_scale
    res[i,2]=mean(line_length(x, y, trace.xad,trace.yad))
    if (i==1){
    }
    else if (i==2&is.na(res[i,2]>res[i-1,2])){
      browser()
    }
    else if (i==2 & res[i,2]>res[i-1,2]) {
        reverse=1
      }
      else if (res[i,2]>res[i-1,2]){
        minimum=1
      }
    i=i+1
  }
  min <- which.min(res[,2])
  scaley <- (trace.y*res[min,1])
  scalex <- (trace.x*res[min,1])
  if(print.plt==TRUE){
    plot_tracing_err(x,y,scalex,scaley)+
      labs(title = paste0("The attained error post scaling. 
                          Scale factor =", as.character(res[min,1]),
                          " Error = ", as.character(res[min,2]),".")
      )
  }
  return(data.frame(x,y,scalex,scaley))
}

#### TRANSLATION ####

fullgrid.trans <- function(x,y,trace.x,trace.y,print.plt){
  rangex <- ((0-round(min(trace.x)))+1):((screen_res[1]-round(max(trace.x)))-1)
  rangey <- ((0-round(min(trace.y)))+1):((screen_res[2]-round(max(trace.y)))-1)
  res <- matrix(NA,ncol=3,nrow=length(rangex)*length(rangey))
  for (i in 1:length(rangex)) {
    trace.xad <- (trace.x+rangex[i])
    for (j in 1:length(rangey)) {
      thisind <- (i-1)*length(rangey)+j
      trace.yad <- (trace.y+rangey[j])
      res[thisind,1]=rangex[i]
      res[thisind,2]=rangey[j]
      res[thisind,3]=mean(line_length(x, y, trace.xad,trace.yad))
    }
  }
  min <- which.min(res[,3])
  thisy <- (trace.y+res[min,2])
  thisx <- (trace.x+res[min,1])
  if(print.plt==T){
    plot_tracing_err(x,y,thisx,thisy)+
      labs(title = paste0("The attained error post transform ",
                          "xTrans =", as.character(res[min,1]),
                          " Ytrans = ", as.character(res[min,2]),
                          " Error = ", as.character(res[min,3]),".")
      )
  }
  return(c(res[min,1],res[min,2],res[min,3]))
}

fullgridvec.trans <- function(x,y,trace.x,trace.y,print.plt){
  rangex <- ((0-round(min(trace.x)))+1):((screen_res[1]-round(max(trace.x)))-1)
  rangey <- ((0-round(min(trace.y)))+1):((screen_res[2]-round(max(trace.y)))-1)
  perms <- cbind(rep(rangex, each = length(rangey)), rep(rangey, length(rangex)))
  trace.xad = rep(1,length(perms[,1])) %*% t(trace.x)
  trace.yad = rep(1,length(perms[,1])) %*% t(trace.y)
  trace.xad = trace.xad + perms[,1]
  trace.yad = trace.yad + perms[,2]
  trace.xad = t(trace.xad)-x
  trace.yad = t(trace.yad)-y
  err = sqrt(trace.xad**2+trace.yad**2)
  err = colMeans(err)
  min <- which.min(err)
  tranx <- (trace.x+perms[min,1])
  trany <- (trace.y+perms[min,2])
  if(print.plt==TRUE){
    plot_tracing_err(x,y,tranx,trany)+
      labs(title = paste0("The attained error post vectorized translation.", 
                          "xTrans =", as.character(perms[min,1]),
                          " Ytrans = ", as.character(perms[min,2]),
                          " Error = ", as.character(err[min]),".")
      )
  }
  return(data.frame(x,y,tranx,trany))
}

partgrid.trans <- function(x,y,trace.x,trace.y,print.plt){
  rangex <- ((0-round(min(trace.x)))+1):((screen_res[1]-round(max(trace.x)))-1)
  rangey <- ((0-round(min(trace.y)))+1):((screen_res[2]-round(max(trace.y)))-1)
  res <- matrix(NA,ncol=3,nrow=length(rangex)*length(rangey))
  reversex=0
  reversey=0
  minimumx=0
  minimumy=0
  i=1
  j=1
  thisind=1
  while (minimumx==0){
    if (reversex==0){
      this_tranx=rangex[which(rangex==0)+(i-1)]
    }
    else {
      this_tranx=max(rangex[which(rangex<this_tranx)])
    }
    trace.xad <- (trace.x+this_tranx)
    while (minimumy==0){
      if (reversey==0){
        this_trany=rangey[which(rangey==0)+(j-1)]
      }
      else {
        this_trany=max(rangey[which(rangey<this_trany)])
      }
      trace.yad <- (trace.y+this_trany)
      res[thisind,1]=this_tranx
      res[thisind,2]=this_trany
      res[thisind,3]=mean(line_length(x, y, trace.xad,trace.yad))
      if (j==1){
      }
      else if (j==2 & res[thisind,3]>res[thisind-1,3]) {
        reversey=1
      }
      else if (res[thisind,3]>res[thisind-1,3]){
        minimumy=1
        thisind=thisind-1
      }
      j=j+1
      thisind=thisind+1
    }
    res[thisind,1]=this_tranx
    res[thisind,2]=this_trany
    res[thisind,3]=mean(line_length(x, y, trace.xad,trace.yad))
    if (i==1 & res[thisind,3]>res[thisind-1,3]) {
      reversex=1
    }
    else if (res[thisind,3]>res[thisind-1,3]){
      minimumx=1
    }
    i=i+1
    thisind=thisind+1
  }
  min <- which.min(res[,3])
  thisy <- (trace.y+res[min,2])
  thisx <- (trace.x+res[min,1])
  if(print.plt==T){
    plot_tracing_err(x,y,thisx,thisy)+
      labs(title = paste0("The attained error post transform ",
                          "xTrans =", as.character(res[min,1]),
                          " Ytrans = ", as.character(res[min,2]),
                          " Error = ", as.character(res[min,3]),".")
      )
  }
  return(c(res[min,1],res[min,2],res[min,3]))
}
