function [im_ds,kk,L,kt,im]=down_samp(k,im2,from,to)
L=1;
kt=k;
if (round(k)~=k)
    kp=k;
    t=0;
    while (round(kp)~=kp)
        kp=kp*10;
        t=t+1;
    end
    g_c_d=gcd(kp,t*10);
    L=t*10/g_c_d;
    M=kp/g_c_d;
    k=M;
end

im=im2(from:to,from:to);
if (ndims(im)==3)
    im=rgb2gray(im);
end

kk=k;

[lx,ly]=size(im);

if (round(kt)~=kt)
    im_n=zeros(round(lx*L/M),round(ly*L/M));
    im_new=zeros(round(lx*L/M),round(ly*L/M));
    for it=1:lx
        jj=1;
        j=1;
        while jj<=ly
            ii=1;
            tt=1;
            i=1;
            if (j*M<=ly)
                im_vec=im(it,(j-1)*M+1:j*M);
            else
                im_vec=im(it,(j-1)*M+1:ly);
            end
            jj=jj+M;
            while ii<L
                while i<=length(im_vec)
                    im_vec_s(ii)=im_vec(i);
                    i=i+2;
                    ii=ii+1;
                end
                i=2;
                im_vec_ss(tt)=im_vec(i);
                i=i+2;
                tt=tt+1;
                ii=ii+1;
            end
            ii=1;
            i=1;
            while ii<=L
                im_vec_sn(ii)=im_vec_s(ii);
                if i<=length(im_vec_ss)
                    im_vec_sn(ii+1)=im_vec_ss(i);
                    ii=ii+1;
                    i=i+1;
                end
                ii=ii+1;
            end
            im_n(it,(j-1)*L+1:j*L)=im_vec_sn;
            j=j+1;
        end
    end
    
    [xx,yy]=size(im_n);
    for jt=1:yy
        ii=1;
        i=1;
        while ii<=lx
            jj=1;
            tt=1;
            j=1;
            if (i*M<=lx)
                im_vec=im_n((i-1)*M+1:i*M,jt);
            else
                im_vec=im_n((i-1)*M+1:lx,jt);
            end
            ii=ii+M;
            while jj<L
                while j<=length(im_vec)
                    im_vec_s(jj)=im_vec(j);
                    j=j+2;
                    jj=jj+1;
                end
                j=2;
                im_vec_ss(tt)=im_vec(j);
                j=j+2;
                tt=tt+1;
                jj=jj+1;
            end
            jj=1;
            j=1;
            while jj<=L
                im_vec_sn(jj)=im_vec_s(jj);
                if j<=length(im_vec_ss)
                    im_vec_sn(jj+1)=im_vec_ss(j);
                    jj=jj+1;
                    j=j+1;
                end
                jj=jj+1;
            end
            im_new((i-1)*L+1:i*L,jt)=im_vec_sn.';
            i=i+1;
        end
    end
    
    im_ds=im_new(1:round(lx*L/M),1:round(ly*L/M));
else
    im_ds=downsample(im,kt);
    im_ds=downsample(im_ds.',kt);
    im_ds=im_ds.';
end
end

