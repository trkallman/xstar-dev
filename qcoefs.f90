!-----------------------------------------------------------------------
       subroutine qcoefs(x1,x2,y1,y2,xx,y)
       rm=(y2-y1)/(x2-x1)
       y=y1+rm*(xx-x1)
       return
       end
