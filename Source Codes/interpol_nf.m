function psnr1=interpol_nf(k,im2,method,from,to)

[im_ds,kk,L,kt,imm]=down_samp(k,im2,from,to);

[~,ly]=size(im_ds);
im_interpol=upsample(im_ds,kk);
im_interpol=upsample(im_interpol.',kk);
im_interpol=im_interpol.';
im_down1=im_interpol;
[xx,yy]=size(im_down1);

for i=1:kk:xx
    for j=1:ly
        for k=1:kk-1
            if j==1
                fxi_1=double(im_down1(i,(j-1)*kk+1));
            else
                fxi_1=double(im_down1(i,(j-2)*kk+1));
            end
            if (j==ly)
                fxi1=double(im_down1(i,(j-1)*kk+1));
                fxi2=double(im_down1(i,(j-1)*kk+1));
            elseif (j==ly-1)
                fxi1=double(im_down1(i,j*kk+1));
                fxi2=double(im_down1(i,j*kk+1));
            else
                fxi1=double(im_down1(i,j*kk+1));
                fxi2=double(im_down1(i,(j+1)*kk+1));
            end
            s=(kk-k)/kk;
            fxi=double(im_down1(i,(j-1)*kk+1));
            if (method==1)
                im_interpol(i,(j-1)*kk+k+1)=(1-s)*fxi+(s)*fxi1;
            elseif method==2
                im_interpol(i,(j-1)*kk+k+1)= fxi_1*(-0.5*(abs(1+s))^3+(5/2)*(1+s)^2-4*abs(1+s)+2)+fxi*(1.5*(abs(s))^3-2.5*s^2+1)+fxi1*(1.5*(abs(s-1))^3-2.5*(s-1)^2+1)+fxi2*(-0.5*(abs(s-2))^3+2.5*(s-2)^2-4*abs(s-2)+2);
            elseif method==3
                im_interpol(i,(j-1)*kk+k+1)=(1/6)*fxi_1*(-(abs(1+s))^3+6*(1+s)^2-12*abs(1+s)+8)+(1/6)*fxi*(3*(abs(s))^3-6*s^2+4)+(1/6)*fxi1*(3*(abs(s-1))^3-6*(s-1)^2+4)+(1/6)*fxi2*(-(abs(s-2))^3+6*(s-2)^2-12*abs(s-2)+8);
            end
            pix_int=im_interpol(i,(j-1)*kk+k+1);
            
        end
    end
end
for j=1:yy
    for i=1:kk:xx-kk+1
        for k=1:kk-1
            s=(kk-k)/kk;
            fxi=double(im_interpol(i,j));
            if i==1
                fxi_1=double(im_interpol(i,j));
            else
                fxi_1=double(im_interpol(i-kk,j));
            end
            if (i==xx-kk+1)
                fxi1=double(im_interpol(i,j));
                fxi2=double(im_interpol(i,j));
            elseif i==xx-2*kk+1
                fxi1=double(im_interpol(i+kk,j));
                fxi2=double(im_interpol(i+kk,j));
            else
                fxi1=double(im_interpol(i+kk,j));
                fxi2=double(im_interpol(i+2*kk,j));
            end
            if (method==1)
                im_interpol(i+k,j)=(1-s)*fxi+(s)*fxi1;
            elseif method==2
                im_interpol(i+k,j)= fxi_1*(-0.5*(abs(1+s))^3+(5/2)*(1+s)^2-4*abs(1+s)+2)+fxi*(1.5*(abs(s))^3-2.5*s^2+1)+fxi1*(1.5*(abs(s-1))^3-2.5*(s-1)^2+1)+fxi2*(-0.5*(abs(s-2))^3+2.5*(s-2)^2-4*abs(s-2)+2);
            elseif method==3
                im_interpol(i+k,j)=(1/6)*fxi_1*(-(abs(1+s))^3+6*(1+s)^2-12*abs(1+s)+8)+(1/6)*fxi*(3*(abs(s))^3-6*s^2+4)+(1/6)*fxi1*(3*(abs(s-1))^3-6*(s-1)^2+4)+(1/6)*fxi2*(-(abs(s-2))^3+6*(s-2)^2-12*abs(s-2)+8);
            end
            im_interpol(i+k,j)=(1-s)*fxi+(s)*fxi1;
            pix_int=im_interpol(i+k,j);
            
        end
    end
end
if (round(kt)==kt)
    figure;
    imshow(im_interpol,[]);
    title('Interpolate Without LEAD');
    if size(im_interpol)==size(imm)
        mse_eror=mse(im_interpol(2:xx-1,2:yy-1)-double(imm(2:xx-1,2:yy-1)));
    else
        mse_eror=mse(im_interpol(2:xx-2,2:yy-2)-double(imm(2:xx-2,2:yy-2)));
    end
end
if (round(kt)~=kt)
    im_d_int=downsample(im_interpol,L);
    im_d_int=downsample(im_d_int.',L);
    im_d_int=im_d_int.';
    figure;
    imshow(im_d_int,[]);
    if size(im_d_int)==size(imm)
        mse_eror=mse(im_d_int(2:xx-1,2:yy-1)-double(imm(2:xx-1,2:yy-1)));
    else
        mse_eror=mse(im_d_int(2:xx-2,2:yy-2)-double(imm(2:xx-2,2:yy-2)));
    end
end
psnr1=10*log10((255^2)/mse_eror);