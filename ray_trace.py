# coding:utf8
'''    Mini ray tracing                 k
  =============================        /
                                      /
       Le repère par défaut :        O------- i
                                     |
       L'observateur est en O        |
       et regarde dans la            j
       direction de k (qui s'"enfonce" dans l'écran) .
       Pour tourner  la tête vers la droite       : jrot( a > 0 )
                                  la gauche       : jrot( a < 0 )
                                  le haut         : irot( a > 0 )
                                  le bas          : irot( a < 0 )
       Pour incliner la tête dans le sens horaire : krot( a > 0 )
                                  l'autre sens    : krot( a < 0 )
       Pour se déplacer le long d'un cercle centré sur un point observé : hrot et vrot
'''

import math

epsilon  = 1e-7
infini   = 1e99
degrad   = math.pi/180
sqrt     = lambda x  : x**0.5
cosin    = lambda a  : (math.cos(a*degrad),math.sin(a*degrad)) # a in degrees
vsum     = lambda it : sum( it,Vector())

def precomputed ( value,func,*args ):
    if value is not None : return value
    res = func( *args )
    if res is not None :
        return res[1]

###################################################################################################

class Vector (list):
    def __init__ (self,*e) :
        for a in ('[0][0]','[0]','') :
            try    : self[:] = (map( float,eval('e'+a ))+[.0]*3 )[:3] ; break
            except : pass
    for o in '+add -sub ^xor |or *mul /div'.split():
        exec('def __i%s__ (S,P) : S[:] = S%sP ; return S'%(o[1:],o[0]))
    __getattr__= lambda S,n : S['xyz'.find(n.lower())%3]
    __str__    = lambda S   : '['+','.join('%g'%x for x in S )+']'
    __repr__   = lambda S   : 'Vector('+','.join('%g'%x for x in S )+')'
    __invert__ = lambda S   : (lambda n : S/(n**0.5) if n else S*0)( S.norm2())         # ~ == normalization
    __abs__    = lambda S   : S.norm2() ** 0.5
    norm2      = lambda S   : sum( x**2 for x in S )
    __and__    = lambda S,P : sum( x*y for x,y in zip( S,P ))                           # & == scalar product
    __xor__    = lambda S,P : Vector( S[i]*P[(i+1)%3]-S[(i+1)%3]*P[i] for i in (1,2,0)) # ^ == cross product
    __or__     = lambda S,P : Vector( x*y for x,y in zip( S,P ))
    __add__    = lambda S,P : Vector( x+y for x,y in zip( S,P ))
    __sub__    = lambda S,P : Vector( x-y for x,y in zip( S,P ))
    __mul__    = lambda S,a : Vector( x*a for x in S )
    __rmul__   = lambda S,a : Vector( x*a for x in S )
    __div__    = lambda S,a : Vector( x/a for x in S )
    __neg__    = lambda S   : Vector(  -x for x in S )
    min        = lambda S,a : Vector( min( x,a ) for x in S )
    map        = lambda S,f : tuple( map( f,S ))
    def iter (self) :
        for z in xrange( int( self.z )):
            for y in xrange( int( self.y )):
                    for x in xrange( int( self.x )):
                        yield Vector( x,y,z )

class Repere :
    def __init__ (S,*m) : S.O,S.i,S.j,S.k = map( Vector,( 0,1,( 0,1 ),( 0,0,1 ))) ; S.move(*m)
    def move  (S,*m) : S.O += Vector( m ) ; return S
    def irot  (S,a)  : c,s = cosin(a) ; S.j,S.k     = S.j*c + S.k*s , S.k*c - S.j*s ; return S
    def jrot  (S,a)  : c,s = cosin(a) ; S.k,S.i     = S.k*c + S.i*s , S.i*c - S.k*s ; return S
    def krot  (S,a)  : c,s = cosin(a) ; S.i,S.j     = S.i*c + S.j*s , S.j*c - S.i*s ; return S
    def vrot  (S,a,d): c,s = cosin(a) ; S.j,S.k,S.O = S.j*c + S.k*s , S.k*c - S.j*s , S.O + (S.k*(1-c) + S.j*s)*d ; return S
    def hrot  (S,a,d): c,s = cosin(a) ; S.i,S.k,S.O = S.i*c + S.k*s , S.k*c - S.i*s , S.O + (S.k*(1-c) + S.i*s)*d ; return S
    def get   (S,c=1): return S.O,S.i*c,S.j*c,S.k*c
    def __str__ (S)  : return '(%s,%s,%s,%s)'%(S.O,S.i,S.j,S.k)

###################################################################################################

class Ray :
    def __init__ ( self,base,dir,color=(1,1,1),inside=None,previous=None ) :
        self.base,self.dir,self.inside,self.color,self.total = Vector( base ),~Vector( dir ),inside or False,Vector(color),0
        if previous :
            self.color = previous.color if color is None else (previous.color|Vector(color))
            self.total = abs( self.base - previous.base ) + previous.total
            if inside is None : self.inside = previous.inside
    def colorize (self,col)  :
        self.color |= col ; return self
    coeff    = lambda self   : 1/(self.total/200+1)
    __call__ = lambda self,t : self.base + self.dir*t
    __repr__ = lambda self   : '(%s,%s)'%(self.base,self.dir)

###################################################################################################

class Light :
    def __init__(self,center,color=(1,1,1),alpha=200):
        self.center,self.color,self.alpha = Vector( center ),Vector( color ),alpha

    def coeff (self,ray) :
        BC,U = self.center-ray.base,ray.dir
        t = BC & U
        #if t < 0 : return 0
        r = abs( BC - t*U ) if t > 0 else abs( BC )
        return 1/( 1 + r*r/self.alpha )  # ==  self.alpha / ( self.alpha + r*r )

###################################################################################################

class Plane :
    def __init__(self,O,I,J,color=(1,1,1),indice=1,opacity=1):
        assert indice >= 1 and opacity >= -1
        self.O,I,J,self.indice,self.opacity = Vector( O ),Vector( I ),Vector( J ),indice,opacity
        self.color = color if callable( color ) else lambda *x : color
        self.I,self.J,self.K = I/I.norm2(),J/J.norm2(),~(I^J)

    def intersect (self,ray,dmin=infini) :
        O,K,B,U = self.O,self.K,ray.base,ray.dir
        t = U&K
        if abs( t ) > epsilon :
            d = ((O-B)&K)/t
            if dmin > d > epsilon :
                I   = ray( d )
                OI  = I-O
                col = self.color( OI&self.I , OI&self.J )
                if col is not None :
                    return d,I,cmp(t,0)*K,Vector( col )

###################################################################################################

class Sphere :
    def __init__(self,center,radius,color=(1,1,1),indice=1,opacity=1):
        assert indice >= 1 and opacity >= -1
        self.center,self.radius,self.indice,self.opacity = Vector( center ),radius,indice,opacity
        self.color = color if callable( color ) else lambda *x : color

    def intersect (self,ray,dmin=infini) : # returns (dist,point,normal,color) or None
        O,r = self.center,self.radius      # where dist   : distance( ray.base , point )
        B,U = ray.base,ray.dir             #       point  : intersection point
        OB  = O-B                          #       normal : normal vector to the sphere at point
        b2  = r*r - ( OB ^ U ).norm2()     #       color  : color of the sphere at point
        if b2 > 0 :
            a,b = OB&U,sqrt( b2 )
            d   = a-b if a > b+epsilon else a+b
            if dmin > d > epsilon :
                I   = ray( d )
                N   = ~(O-I)
                col = self.color( N )
                if col is not None :
                    return d,I,N,Vector( col )

###################################################################################################

def reflect (object,ray,I,N,col) : # I = intersection point , N = normal at this point
    U       = ray.dir              # returns tuple of rays
    V       = (N&U)*N
    refl    = Ray( I,U-2*V,previous=ray )
    if object.opacity >= 1 : # fully opaque : no transmission at all
        return (refl.colorize(col),)
    T       = ~(U-V)
    i,s     = (object.indice,-1) if ray.inside else (1./object.indice,+1)
    sina    = T&U
    sinb    = i*sina
    if abs( sinb ) >= 1 :
        return (refl.colorize(col),)
    a       = max( sina,sinb )
    if object.opacity > 0 : a = (1-a)*object.opacity + a
    else                  : a = (1+object.opacity) * a # opacity == -1 --> fully transparent
    if a > 1-epsilon :
        return (refl.colorize(col),)
    cosb    = s*sqrt( 1 - sinb**2 )
    R       = N*cosb + T*sinb
    trans   = Ray( I,R,inside=not ray.inside,previous=ray )
    if a < epsilon : # fully transparent : no reflection at all
        return (trans.colorize( col ),)
    return (refl.colorize(col*a),trans.colorize(col*(1-a)))

def trace1ray ( objects,lights,ray,deep ):
    imin,omin,col = (infini,0),None,ray.color
    for obj in objects :
        i = obj.intersect( ray,imin[0] ) # i = (dist,point,normal,color)
        if i and epsilon < i[0] < imin[0] :
            imin,omin = i,obj
    if not omin :
        return vsum( l.color*l.coeff( ray ) for l in lights )|col
    rays = reflect( omin,ray,*imin[1:] )
    if deep == 0 :
        return vsum( l.color*l.coeff( r ) for l in lights for r in rays )|col
    if omin.opacity > 1+epsilon :
        deep = 1
    return vsum( trace1ray( objects,lights,r,deep-1 ) for r in rays )|col

def viewbox ( size,obs,aperture_deg ):
    a = math.tan( aperture_deg*degrad/2 )
    I = obs.i*a
    J = obs.j*a*size[1]/size[0]
    O = obs.k - I - J
    return O,2*I/max(1,size[0]-1),2*J/max(1,size[1]-1)

def TraceRays ( size,scene,obs,aperture,plot,deepmax=8 ):
    global debug
    lights  = filter( lambda s : hasattr( s,'coeff'),scene )
    objects = filter( lambda s : hasattr( s,'intersect'),scene )
    O,I,J   = viewbox( size,obs,aperture )
    for y in xrange( size[1] ):
        for x in xrange( size[0] ):
            debug = x == 95 and y == size[1]/2
            ray = Ray( obs.O,O + x*I + y*J )
            col = trace1ray( objects,lights,ray,deepmax )
            if not plot((x,y),col ) : return

###################################################################################################

class PlotterTk :
    import Tkinter,thread,time
    def __init__(self,size,colorize,text='',zoom=5,**args):
        self._pos = lambda (x,y) : (x,y)
        self._col = lambda c     : '#%02x%02x%02x'%colorize(c)
        if zoom > 1 :
            self._pos = lambda (x,y) : (x*zoom,y*zoom)
            self._col = lambda c     : (lambda s : (' {'+s*zoom+' }')*zoom)(' #%02x%02x%02x'%colorize(c))
        self.root = self.Tkinter.Tk()
        self.img  = self.Tkinter.PhotoImage( width=size[0]*zoom,height=size[1]*zoom )
        self.lock = self.thread.allocate_lock()
        self.t0   = self.time.time()
        self.root.title('%s (press ESCAPE to quit)'%( text or 'running'))
        self.Tkinter.Label( self.root,image=self.img ).pack()
        self.root.bind('<Escape>',lambda e : self.root.quit())
        self.root.geometry('+600+400')
        self.thread.start_new_thread( lambda : (self.lock.acquire(),self.root.mainloop(),self.lock.release()),() )
    __call__ = lambda self,pos,color : self.lock.locked() and [self.img.put( self._col( color ),self._pos( pos ))]
    _time    = lambda self           : self.root.title('%ds : press ESCAPE to quit'%(self.time.time()-self.t0))
    finish   = lambda self           : self.root.after( 1,self._time ) is self.lock.acquire()

class TimeCounter :
    import timeit
    def __init__(self):
        self.t0 = self.timeit.default_timer()
    def __call__ ( self,ratio ) :
        def timedelta ( seconds ):
            s = int( seconds )
            m,h,j = s/60,s/3600,s/(3600*24)
            if j : return '%dj%02dh'%( j,h%24 )
            if h : return '%dh%02dm'%( h,m%60 )
            if m : return '%dm%02ds'%( m,s%60 )
            return '%d.%02ds'%( s%60,(seconds-s)*100 )
        if ratio == 0 : return '?','?','?'
        t,r = self.timeit.default_timer() - self.t0 , 1/ratio
        return timedelta( t ),timedelta( t*( r-1 )),timedelta( t*r )

class PlotterPIL :
    from PIL import Image
    import datetime
    def __init__(self,size,colorize,filename='',text='',**args ):
        self.col  = colorize
        self.name = filename or str( self.datetime.datetime.now())[:19].replace(':','-').replace(' ','_')+('_' if text else '')+text
        self.img  = self.Image.new('RGB',size)
        self.pct  = [-1,0,size[0]*size[1]]
        self.tc   = TimeCounter()
    def __call__ ( self,pos,color ):
        self.pct[1] += 1
        p = (1000*self.pct[1])/self.pct[2]
        if p != self.pct[0] :
            print '%.1f%%    %s + %s = %s  \r'%(( p/10.0,)+self.tc( 1.0*self.pct[1]/self.pct[2] )), 
            self.pct[0] = p
        self.img.putpixel( pos,self.col( color ))
        return 1
    def finish (self):
        self.img.save( self.name+'.png' )

###################################################################################################

def test_lut ( nb=60 , maxx=1800 ):
    print 
    maxy     = max( 100,maxx/nb )
    colorize = lambda v   : (int(v[0]*255),int(v[1]*255),int(v[2]*255))
    rnd      = lambda i,j : 1111111111111111111/(1+i*3+j)
    comp     = lambda i,j : (1001 + rnd(i,j)%2000)/3000.0
    lut      = [Vector( comp(i,1),comp(i,2),comp(i,3)) for i in xrange( nb )]
    plot     = PlotterTk( (maxx,maxy),colorize,'test the LUT')
    for x,y in ( (x,y) for y in xrange( maxy ) for x in xrange( maxx )):
        plot( (x,y),lut[ x*nb/maxx ] )
    plot.finish()
    
if 0 :
    test_lut( 200 )
    quit() 

def test ():
    globals().update({'C%x%x%x'%tuple(v):v/9 for v in Vector( 10,10,10 ).iter()})
    f2byte   = lambda x         : min( 255,max( 0,int( x*255 )))
    colorize = lambda c=1       : lambda v : tuple( map( f2byte,(v*c)[:3]))
    checker  = lambda x,y,d=1   : min( 1 , int((x/d)%2) ^ int((y/d)%2) )
    nchecker = lambda x,y,n,d=1 : ( int((x/d)%n) ^ int((y/d)%n))%n
    rnd      = lambda i,j       : 1111111111111111111/(1+i*i*3+j)
    comp     = lambda i,j       : (1001 + rnd(i,j)%2000)/3000.0
    make_lut = lambda nb        : [Vector( comp(i,1),comp(i,2),comp(i,3)) for i in xrange( nb )]
    I,J,K    = Vector(1),Vector(0,1),Vector(0,0,1)
    def D ( x,*y ) :
        print ' '.join( map( str,y ))
        return x

    class Mandel :
        def __init__(self,pos=(0,0),coeff=1,max=5,start=0):
            f = lambda i,j : (1001 + (i*(j*100+73))%2000)/3000.0
            self.m,self.a,self.b = max,coeff*1.,complex( *pos )
            self.lut = [None]*start+[Vector( f(i,1),f(i,2),f(i,3)) for i in xrange( max )] + [Vector()]
        def __call__( self,x,y ):
            p,z = self.a*( x*1j+y )+self.b,0
            for n in xrange( self.m ):
                if abs( z ) > 2 :
                    return self.lut[n]
                z = z*z + p
            return self.lut[-1]

    def lens_z ( center,r,e,c,i,o=0 ): # lens along z axis : r = lens radius ; e = lens thickness
        R = (r*r + e*e/4.)/e
        u = 1 - e/(2*R)
        z = Vector(0,0,u*R)
        try   : c1,c2,x = c[0],c[1],c[0][0]
        except: c1,c2   = c,c
        scol1 = lambda v : c1 if v.z >  u else None
        scol2 = lambda v : c2 if v.z < -u else None
        return [Sphere( center+z,R,scol1,indice=i,opacity=o ),
                Sphere( center-z,R,scol2,indice=i,opacity=o )]

    def image10 ( size,Plotter ):
        plot  = Plotter( size,colorize())
        dcol  = lambda x,y : (C966,C669)[checker(x,y,3)] if -8<x<8 and -8<y<8 else None
        scol1 = lambda v   : C999 if v.z >  0.5 else None
        scol2 = lambda v   : C999 if v.z < -0.5 else None
        scene = [Light((0,0,42),alpha=300),
                 Sphere((0,0,20),5,scol1,1.1 ),
                 Sphere((0,0,15),5,scol2,1.1 ),
                 Plane((0,0,30),1,(0,1),dcol,opacity=1 )]
        TraceRays( size,scene,Repere().hrot(5,20).vrot(2,20),45,plot,3 )
        plot.finish()

    def image11 ( size,Plotter ):
        plot  = Plotter( size,colorize())
        dcol1 = lambda a   : lambda x,y : a if 0<x<8 and 0<y<8 else None
        dcol2 = lambda a,b : lambda x,y : (a,b)[checker(x,y,2)] if 0<x<8 and 0<y<8 else None
        scene = [Light((4,-4,-4),alpha=100),
                 Sphere((4,-4,-4),1,C999,indice=1.1,opacity=1 ),
                 #Sphere((0,0,20),5,scol1,1.1 ),
                 #Sphere((0,0,15),5,scol2,1.1 ),
                 #Sphere(0,1,C666,1 ),Sphere(8,1,C900,1 ),Sphere((0,-8),1,C090,1 ),Sphere((0,0,-8),1,C009,1 ),
                 Plane(0,(1, 0,0),(0,-1,0),dcol1(C977),opacity=1 ),
                 Plane(0,(0,-1,0),(0,0,-1),dcol1(C799),opacity=1 ),
                 Plane(0,(1, 0,0),(0,0,-1),dcol2(C666,C444),opacity=-1 )]
        TraceRays( size,scene,Repere(0,-2,-20).hrot(35,20).vrot(-15,20).jrot(3),50,plot,5 )
        plot.finish()

    def image12 ( size,Plotter,d=3 ):
        plot  = Plotter( size,colorize(2))#,'zz_loupe_%d'%int(d*1000))
        m     = 5
        rnd   = lambda n : ((1111111121/(n+1))%101+50)/151.0
        lut   = [None]*3 + [Vector(rnd(i),rnd(i+m),rnd(i+m*2)) for i in xrange( m )]
        dcol  = lambda a,b : lambda x,y : (a,b)[checker(x,y)] if -4<x<4 and -4<y<4 else None
        miro  = lambda x,y : C977 if -3<x<3 and -3<y<3 else None
        scol1 = lambda v   : C999 if v.z >  0.5 else None
        scol2 = lambda v   : C999 if v.z < -0.5 else None
        def mandel ( x,y ):
            p,z = complex( y-.75,x ),0
            for n in xrange( m ):
                if abs( z ) > 2 : return lut[n]
                z = z*z + p
            return Vector()
        scene = [Light((0,0,10),alpha=300),
                 Plane((0,0,-10),(1,0,-1),(0,1),miro,opacity=1 ),
                 Plane((0,0,.1),1,(0,1),dcol(C666,C444),opacity=-1 ),
                 Plane( 0,1,(0,1),mandel,opacity=1 ),
                 Plane( 4,1,(0,1),mandel,opacity=1 ),
                 Sphere((0,0,-d)  ,2,scol1,indice=1.27,opacity=0 ),
                 Sphere((0,0,-d-2),2,scol2,indice=1.27,opacity=0 )]
        TraceRays( size,scene,Repere(0,0,-30).hrot(87,20).jrot(15),40,plot,3 )
        #TraceRays( size,scene,Repere(0,0,-10),50,plot,3 )
        plot.finish()

    def image13 ( size,Plotter ):
        plot  = Plotter( size,colorize(2))
        dcol  = lambda x,y : (C777,C333)[checker(x,y)]
        scene = [Light((0,0,10),alpha=300),Plane(0,J,K,dcol,opacity=-1 )]
        for x in xrange(5) :
            scene += lens_z( Vector(0,0,(x-2)*2),3,1+x/4.,C933,1,-1 )
        TraceRays( size,scene,Repere(0,0,-90).hrot(90,90),7,plot,3 )
        plot.finish()

    def image14 ( size,Plotter ):
        plot  = Plotter( size,colorize(2))
        dcol  = lambda a,b : lambda x,y : (a,b)[checker(x,y,.05)] if -8<x<8 and -8<y<8 else None
        scene = [Light(K*10,alpha=300),
                 Plane(K*5,I,J,dcol(C777,C333))] #+ lens_z( 0,0,0,3,2,C977,1.1 )
        for pos in Vector( 5,3,1 ).iter() :
            scene += lens_z( pos*2-(4,2,0),0.9,0.5+pos.x/4.,(C966,C699),1.1+pos.y/9.,-0.8 )
        TraceRays( size,scene,Repere(0,0,-90),7,plot,5 )
        plot.finish()

    def image15 ( size,Plotter ):
        plot  = Plotter( size,colorize(1))
        dcol  = lambda a,b : lambda x,y : (a,b)[checker(x,y,.05)] if -8<x<8 and -8<y<8 else None
        scene = [Light(K*16,alpha=300),
                 Plane(K*5,I,J,Mandel((-1.7111,0),0.01,99))]
        scene += lens_z( 0,0.9,0.5,C999,1.2,-0.8 )
        TraceRays( size,scene,Repere(0,0,-9),15,plot,2 )
        plot.finish()

    def image16 ( size,Plotter ):
        plot  = Plotter( size,colorize(1))
        dcol  = lambda *cols : lambda x,y : cols[nchecker(x,y,len( cols ),.3)] if -8<x<8 and -8<y<8 else None
        miro  = lambda  d    : lambda x,y : C889 if -d<x<d and -d<y<d else None
        scene = [Light(K*16,alpha=300),
                 Plane(-1*K,(I-K)/1.41,J,miro( 3 )),
                 #Plane(K*5,I,J,dcol(C333,C399,C939)),
                 Plane(K*5,I,J,Mandel((-1.7111,0),0.01,99)),
                 Plane(K*5,I,J,dcol(C333,C399,C939,C993,C339,C393,C933,C999,C963)),
                 ]
        if 0 :
            for i in xrange( 1,5 ) :
                scene += [Sphere((i,0,0),0.1,C900 ),Sphere((0,i,0),0.1,C090 ),Sphere((0,0,i),0.1,C009 )]
        if 1 :
            scene += lens_z( 0,0.9,0.5,C999,1.2,-0.8 )
            for pos in Vector( 5,5,1 ).iter() :
                scene += [Sphere( pos-(2,2,4),0.1,C191 )]
        TraceRays( size,scene,Repere(0,0,-20).hrot(80,20).vrot(3,20),30,plot,6 )
        plot.finish()

    def image17 ( size,Plotter ):
        plot  = Plotter( size,colorize(1))
        scene = [Light(K*6,alpha=300),Sphere( -2*I,1.5,C977,opacity=0 ),Sphere( 2*I,1.5,C779,opacity=0 )]
        TraceRays( size,scene,Repere(-9*K),25,plot,4 )
        plot.finish()

    def image18 ( size,Plotter ):
        dcol   = lambda x,y : (C333,C999)[checker(x,y,.5)] if -4<x<4 and -4<y<4 else None
        plot   = Plotter( size,colorize(0.51))
        scene  = [Light(I*3,C999,4),Light(J*3,C999,4),Light(-I*3,C999,4)]
        scene  = [Light(I*3,C900,4),Light(J*3,C090,4),Light(-I*3,C009,4)]
        scene += [Plane( 5*K,I,J,dcol,opacity=2 )]
        #scene += [Sphere(-2*I,1.5,C977,indice=1.2,opacity=2 ),Sphere( 2*I,1.5,C797,indice=1.2,opacity=1 )]
        #scene += [Sphere(-2*J,1.5,C779,indice=1.2,opacity=0 ),Sphere( 2*J,1.5,C886,indice=1.2,opacity=-0.5 )]
        TraceRays( size,scene,Repere(-9*K),55,plot,4 )
        plot.finish()

    def image19 ( size,Plotter ):
        a = 0.1
        plot  = Plotter( size,colorize(90))
        scene = [Light(-K*500,C111,1000),Light(I*3,C900,a),Light(J*3,C090,a),Light(-I*3,C009,a),Sphere(K*3,1,C555,opacity=2 )]
        TraceRays( size,scene,Repere(),45,plot,1 )
        plot.finish()

    def image20 ( size,Plotter ):
        a = 10
        plot  = Plotter( size,colorize(1))
        dcol  = lambda x,y : (C666,C999)[checker(x,y,.2)] if -1<x<1 and -1<y<1 else None
        scene = [Light(I*3,C900,a),Light(J*3,C090,a),Light(-I*3,C009,a),Plane( 5*K,I,J,dcol,opacity=1.1 )]
        TraceRays( size,scene,Repere(),45,plot,1 )
        plot.finish()

    def image21 ( size,Plotter ):
        dcol  = lambda a,b : lambda x,y : (a,b)[checker(x,y,.5)] if -4<x<4 and -4<y<4 else None
        plot  = Plotter( size,colorize(1.5))
        scene = [Light(-2*K+4*I,alpha=200),
                 Plane( 6*K,I,J,dcol(C222,C777),opacity=1 ),
                 Plane( 5*K,I,J,Mandel((-0.5,0),0.4,999,5)),
                 ]+lens_z( (1,0),1.5,0.5,(0.95,1,1),2,-0.7 ) # lens along z axis : r = lens radius ; e = lens thickness
        TraceRays( size,scene,Repere(-12*K),45,plot,5 )
        plot.finish()

    def image21a ( size,Plotter,text='',indice=1.5,dist=0.5,repere=Repere(-10*K)):
        dcol  = lambda a,b : lambda x,y   : (a,b)[checker(x,y,.5)] if -4<x<4 and -4<y<4 else None
        icol  = lambda a,b : lambda x,y   : (a,b)[checker(x,y,.5)]
        scol  = lambda a,b : lambda v     : (a,b)[checker(abs(v.x+2)**2.9,abs(v.y+2)**2.9,1)]
        plot  = Plotter( size,colorize(1.5),text=text )
        scene = [Light( 100*K+40*I,alpha=600),
                 #Light( -10*K-20*J,alpha=100),
                 Light( I,alpha=100),
                 #Plane( 6*K,I,J,dcol(C222,C777),opacity=1 ),
                 #Plane( 10*K,I,J,icol(C581,C518),opacity=1 ),
                 Sphere( 100*K,45,scol(C586,C568),opacity=1 ),
                 #Sphere( 4.9*K+0.5*I-0.55*J,0.2,C911,opacity=1 ),
                 Plane( 5*K,I,J,Mandel((-0.5,0),0.4,999,5)),
                 ]+lens_z( dist*K+0.5*I-0.55*J,1.5,0.5,(0.99,1,1),indice,-0.7 ) # center,radius,thickness,color,refr.index,transparency
        TraceRays( size,scene,repere,45,plot,5 )
        plot.finish()

    #image20( (240,150),PlotterTk )
    #image21( (1920,1200),PlotterPIL )
    #image21a( (192,120),PlotterTk )
    #image21a( (1920,1200),PlotterPIL )
    for i in xrange(15):
        print '------ repere =',i,'      '
        image21a( (192,120),PlotterPIL,'repere=%d'%i,repere=Repere( (-10-15+i)*K ))
    return
    for i in xrange(-10,11):
        print '------ dist =',i/20.0
        image21a( (192,120),PlotterPIL,'dist=%d'%(i+10),dist=i/20.0 )

if __name__ == '__main__' :
    test()
