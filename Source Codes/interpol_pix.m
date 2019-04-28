function [fit_nes,out_taw1,ftx]=interpol_pix(im_inp,MM,kk,chr_inp,out_taw,ci,ci1,method)

fxi=double(im_inp(2));
fxi1=double(im_inp(3));
fxi_1=double(im_inp(1));
fxi2=double(im_inp(4));

if (out_taw~=0)
    out_taw1=out_taw;
else
    out_taw1=fit_fuzzy(chr_inp,ci,ci1,MM);
end
fit_nes=0;

if (fxi==255)
    x1=5;
elseif (fxi==0)
    x1=-5;
else
    x1=(-1/5)*log((255-fxi)/fxi);
end
if fxi1==255
    x2=5;
elseif fxi1==0
    x2=-5;
else
    x2=(-1/5)*log((255-fxi1)/fxi1);
end

err=0;
ftx=zeros(1,kk-1);
if method==1
    pt=8;
elseif method==2
    pt=3;
else
    pt=5;
end
if (out_taw == 0)
    if chr_inp(2*MM-1)==0
        for p=1:pt-1
            s=(pt-p)/pt;
            [s_prim]=fit_eval(out_taw1,s,ci,ci1);
            s=s_prim;
            if (method==1)
                im_interpol=(1-s)*fxi+(s)*fxi1;
            elseif method==2
                im_interpol= fxi_1*(-0.5*(abs(1+s))^3+(5/2)*(1+s)^2-4*abs(1+s)+2)+fxi*(1.5*(abs(s))^3-2.5*s^2+1)+fxi1*(1.5*(abs(s-1))^3-2.5*(s-1)^2+1)+fxi2*(-0.5*(abs(s-2))^3+2.5*(s-2)^2-4*abs(s-2)+2);
            elseif method==3
                im_interpol=(1/6)*fxi_1*(-(abs(1+s))^3+6*(1+s)^2-12*abs(1+s)+8)+(1/6)*fxi*(3*(abs(s))^3-6*s^2+4)+(1/6)*fxi1*(3*(abs(s-1))^3-6*(s-1)^2+4)+(1/6)*fxi2*(-(abs(s-2))^3+6*(s-2)^2-12*abs(s-2)+8);
            end
            x=p*(x2-x1)/pt+x1;
            fx_sigt=255/(1+exp(-5*x));
            ftxt=im_interpol;
            mse_er=mse(fx_sigt-ftxt);
            err=mse_er+err;
        end
        if err==0
            err=eps;
        end
        fit_nes=1/err;
    else
        fit_nes=chr_inp(2*MM-1);
    end
end

for k=1:kk-1
    s=(kk-k)/kk;
    [s_prim]=fit_eval(out_taw1,s,ci,ci1);
    s=s_prim;
    if (method==1)
        im_interpol=(1-s)*fxi+(s)*fxi1;
    elseif method==2
        im_interpol= fxi_1*(-0.5*(abs(1+s))^3+(5/2)*(1+s)^2-4*abs(1+s)+2)+fxi*(1.5*(abs(s))^3-2.5*s^2+1)+fxi1*(1.5*(abs(s-1))^3-2.5*(s-1)^2+1)+fxi2*(-0.5*(abs(s-2))^3+2.5*(s-2)^2-4*abs(s-2)+2);
    elseif method==3
        im_interpol=(1/6)*fxi_1*(-(abs(1+s))^3+6*(1+s)^2-12*abs(1+s)+8)+(1/6)*fxi*(3*(abs(s))^3-6*s^2+4)+(1/6)*fxi1*(3*(abs(s-1))^3-6*(s-1)^2+4)+(1/6)*fxi2*(-(abs(s-2))^3+6*(s-2)^2-12*abs(s-2)+8);
    end
    ftx(k)=im_interpol;
end
end