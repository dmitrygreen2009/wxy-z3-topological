function infinite_pairs(family,w;ordering="matter_first",cell_slices=2)
    @assert cell_slices>=1
    @assert family in ["zigzag","armchair"]
    n=9w*cell_slices
    function site(t,j,k)
        jj=mod(j,w)
        if ordering=="star"
            return k<=3 ? 9w*t+9jj+k : k<=6 ? 9w*t+9jj+k+2 :
                k==8 ? 9w*(t-1)+9jj+9 : 9w*t+9jj+(k==7 ? 4 : 5)
        end
        @assert ordering=="matter_first"
        return k<=6 ? 9w*t+6jj+k : k==8 ? 9w*(t-1)+6w+3jj+3 :
            9w*t+6w+3jj+(k==7 ? 1 : 2)
    end
    pairs=[]
    for t=0:cell_slices-1,j=0:w-1,s=0:1,a=1:3,i=1:3
        m=site(t,j,3s+a)
        gt,gj=s==0 ? (t,j) : i==1 ? (t,j) : i==2 ? (t+1,j) :
            family=="zigzag" ? (t,j+1) : (t-1,j+1)
        g=site(gt,gj,6+i)
        shift=-fld(min(m,g)-1,n)*n
        push!(pairs,(m=m+shift,g=g+shift,a=a,i=i))
    end
    pairs
end

function infinite_opsum(family,w;periodic=false,ordering="matter_first",cell_slices=2)
    @assert cell_slices>=1
    @assert family in ["zigzag","armchair"]
    n=9w*cell_slices;os=OpSum()
    for p in infinite_pairs(family,w;ordering,cell_slices)
        m=periodic ? mod1(p.m,n) : p.m;g=periodic ? mod1(p.g,n) : p.g
        os += -W[p.a,p.i],"S-",m,"S+",g
        os += -conj(W[p.a,p.i]),"S+",m,"S-",g
    end
    if !periodic
        for j=1:n;os += 0.0,"I",j;end
    end
    os
end
