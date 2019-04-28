function psnr=main(k,method,im2,use_table,from,to)

if method==1
    M_M=3;
elseif method==2
    M_M=5;
else
    M_M=3;
end
[im_ds,kk,L,kt,imm]=down_samp(k,im2,from,to); %down sampling with k factor


im_test=im_ds;

[xx,yy]=size(im_test);
im_new=zeros(kk*xx,kk*yy);
if use_table==1 % using lookup tables or not
    if method==1
        table_t=load('table2.mat');
        table_t=table_t.table2;
    elseif method==2
        table_t=load('table3.mat');
        table_t=table_t.table3;
    else
        table_t=load('table4.mat');
        table_t=table_t.table4;
    end
else
    table_t=zeros(256,256);
end
cc=0;
tt=0;
for i=1:xx   % do inteprolaion in horizontal dorection
    for j=1:yy
        if(j==1)
            fxi_1=double(im_test(i,j));
        else
            fxi_1=double(im_test(i,j-1));
        end
        if(j==yy-1)
            fxi2=double(im_test(i,j+1));
            fxi1=double(im_test(i,j+1));
        elseif (j==yy)
            fxi2=double(im_test(i,j));
            fxi1=double(im_test(i,j));
        else
            fxi2=double(im_test(i,j+2));
            fxi1=double(im_test(i,j+1));
        end
        
        im_inp=[fxi_1,im_test(i,j),fxi1,fxi2];
        fxi=im_inp(2);
        
        ci=round(abs(fxi1-fxi_1));
        ci1=round(abs(fxi2-fxi));
        if table_t(ci+1,ci1+1)~=0
            out_taw=table_t(ci+1,ci1+1);
            [~,~,ftx]=interpol_pix(im_inp,M_M,kk,0,out_taw,ci,ci1,method);
        else
            [pop,I,M_M,table_t]=Genetic(im_inp,kk,table_t,ci,ci1,method,M_M);
            chr_inp=pop(I,:);
            [~,~,ftx]=interpol_pix(im_inp,M_M,kk,chr_inp,0,ci,ci1,method);
            cc=cc+1;
        end
        im_new(kk*(i-1)+1,(j-1)*kk+1)=im_test(i,j);
        im_new(kk*(i-1)+1,(j-1)*kk+2:j*kk)=ftx;
    end
end
[xx,yy]=size(im_new);
for j=1:yy   %vertical direct
    for i=1:kk:xx-kk+1
        if(i==1)
            fxi_1=double(im_new(i,j));
        else
            fxi_1=double(im_new((i-kk),j));
        end
        if(i==xx+1-2*kk)
            fxi2=double(im_new(i+kk,j));
            fxi1=double(im_new(i+kk,j));
        elseif (i==xx+1-kk)
            fxi2=double(im_new(i,j));
            fxi1=double(im_new(i,j));
        else
            fxi2=double(im_new(i+2*kk,j));
            fxi1=double(im_new(i+kk,j));
        end
        im_inp=[fxi_1,im_new(i,j),fxi1,fxi2];
        fxi=double(im_inp(2));
        
        ci=round(abs(fxi1-fxi_1));
        ci1=round(abs(fxi2-fxi));
        if table_t(ci+1,ci1+1)~=0
            out_taw=table_t(ci+1,ci1+1);
            [~,~,ftx]=interpol_pix(im_inp,M_M,kk,0,out_taw,ci,ci1,method);
        else
            [pop,I,M_M,table_t]=Genetic(im_inp,kk,table_t,ci,ci1,method,M_M);
            chr_inp=pop(I,:);
            [~,~,ftx]=interpol_pix(im_inp,M_M,kk,chr_inp,0,ci,ci1,method);
            tt=tt+1;
        end
        im_new(i+1:i+kk-1,j)=ftx.';
    end
end

if round(kt)==kt
    figure;
    imshow(im_new,[]);
    title('Interpolate with LEAD');
    if size(im_new)==size(imm)
        mse_er=mse(im_new(2:xx-1,2:yy-1)-double(imm(2:xx-1,2:yy-1)));
    else
        mse_er=mse(im_new(2:xx-2,2:yy-2)-double(imm(2:xx-2,2:yy-2)));
    end
else
    im_d_int1=downsample(im_new,L);
    im_d_int1=downsample(im_d_int1.',L);
    im_d_int1=im_d_int1.';
    figure;
    imshow(im_d_int1,[]);
    if size(im_d_int1)==size(imm)
        mse_er=mse(im_d_int1(2:xx-1,2:yy-1)-double(imm(2:xx-1,2:yy-1)));
    else
        mse_er=mse(im_d_int1(2:xx-2,2:yy-2)-double(imm(2:xx-2,2:yy-2)));
    end
end

psnr=10*log10((255^2)/mse_er);

if use_table==1
    if method==1
        table2=table_t;
        save('table2.mat','table2');   
    elseif method==2
        table3=table_t;
        save('table3.mat','table3');
    else
        table4=table_t;
        save('table4.mat','table4');
    end
end

