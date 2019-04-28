function GUI_MV()
%% Grafical User Interface----vahid asbaghi- machine vision project
%
%
%% ************************************************************************
close all;
clear all;
clc;

%Controls Colors Define
panel_color=[0.55 0.75 0.65];
entryField_color=[1 1 1];
button_color=[0.4 0.4 0.9];
res_color=[0.2 0.2 0.3];

%% Handles and Controls Define
hFigure=figure(...
    'Units','Pixels',...
    'Position',[100 100 400 400],...
    'Toolbar','none',...
    'MenuBar','none',...
    'NumberTitle','off',...
    'Color',[.5 .5 .5],...
    'Name','Image Enhancement Genetic Algorithm');

hPanel=uipanel(...
    'Parent', hFigure,...
    'Units','Pixels',...
    'Position',[0 0 400 400],...
    'BackgroundColor',panel_color);

himlist=uicontrol(...
    'Style','listbox',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[10 305 110 70],...
    'String','Fishing_Boat|Lenna|Baboon',...
    'BackgroundColor',entryField_color);

himlistl=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[120 290 100 100],...
    'HorizontalAlignment','center',...
    'String','Select Test Image. Press "Show image" if Want to See Image.',...
    'BackgroundColor',panel_color);

hButton1=uicontrol(...
    'Style','pushbutton',...
    'Parent',hPanel,...
    'Units','Pixels',...
    'Position',[10 375 110 20],...
    'String','Show image',...
    'BackgroundColor',button_color,...
    'Callback',@Showimage_callback);
himlist1=uicontrol(...
    'Style','listbox',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[240 305 110 70],...
    'String','Bilinear|Bicubic|Cubic Spline',...
    'BackgroundColor',entryField_color);
hlistl1=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[240 375 100 17],...
    'String','Select Method',...
    'BackgroundColor',panel_color);
hButton=uicontrol(...
    'Style','pushbutton',...
    'Parent',hPanel,...
    'Units','Pixels',...
    'Position',[5 250 100 40],...
    'String','With EALD',...
    'BackgroundColor',button_color,...
    'Callback',@Start_callback);
hButton2=uicontrol(...
    'Style','pushbutton',...
    'Parent',hPanel,...
    'Units','Pixels',...
    'Position',[110 250 100 40],...
    'String','Without EALD',...
    'BackgroundColor',button_color,...
    'Callback',@Start_callback2);
hRa1=uicontrol(...
    'Style','Radio',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[10 210 100 17],...
    'String','Using Table',...
    'Value',1,...
    'BackgroundColor',panel_color);
hselpix11=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[210 220 120 25],...
    'String','Specify Expriment Pixels From----------------To',...
    'BackgroundColor',panel_color);
hfromE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[185 190 80 25],...
    'String','150',...
    'BackgroundColor',entryField_color);
htoE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[270 190 80 25],...
    'String','450',...
    'BackgroundColor',entryField_color);
hkl=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[5 170 50 25],...
    'String','K value',...
    'BackgroundColor',panel_color);
hkE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[60 175 70 25],...
    'String','2',...
    'BackgroundColor',entryField_color);
hselpix11=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[5 140 400 25],...
    'String','******************************** Results Section **************************************',...
    'BackgroundColor',panel_color);
hPSNRl=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[5 115 50 25],...
    'String','PSNR',...
    'BackgroundColor',panel_color);
hPSNRE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[55 120 70 25],...
    'String','0',...
    'BackgroundColor',entryField_color);
hMl=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[5 85 50 25],...
    'String','M',...
    'BackgroundColor',panel_color);
hME=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[55 90 70 25],...
    'String','0',...
    'BackgroundColor',entryField_color);
hPtl=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[5 55 50 25],...
    'String','Pt',...
    'BackgroundColor',panel_color);
hPtE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[55 60 70 25],...
    'String','0',...
    'BackgroundColor',entryField_color);
htawl=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[150 100 30 25],...
    'String','Taw',...
    'BackgroundColor',panel_color);
htawE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[185 105 70 25],...
    'String','0',...
    'BackgroundColor',entryField_color);
hBetal=uicontrol(...
    'Style','Text',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[150 70 30 25],...
    'String','beta',...
    'BackgroundColor',panel_color);
hBetaE=uicontrol(...
    'Style','Edit',...
    'Parent',hPanel,...
    'Units','Pixel',...
    'Position',[185 75 70 25],...
    'String','0',...
    'BackgroundColor',entryField_color);

handles=[hFigure,hPanel,hfromE,hRa1,...
    hButton,hButton1,himlist,himlist1,hButton2];

set(handles,...
    'Units','Normalized'); %normalized position units: from (0,0) to (1,1)

%% Callbacks for Buttons

%callback for button 'With EALD'
    function Start_callback(hObject,eventdata)        
        %reading parameters
        imname=get(himlist,'Value'); %get selected value number from listbox
        Meth=get(himlist1,'Value');
        use_table = get(hRa1,'Value');
        from=str2double(get(hfromE,'String'));
        to=str2double(get(htoE,'String'));
        k=str2double(get(hkE,'String'));
        
        switch Meth
            case 1
                method=1;
            case 2
                method=2;
            case 3
                method=3;
        end
        
        %read selected image
        switch imname
            case 1
                imtest=imread('fishingboat.tif');
            case 2
                imtest=imread('lena_gray_512.tif');
            case 3
                imtest=imread('baboon_512.tif');
        end
        
        % send required value to 'main' function for evaluate 'PSNR'
        % and show image
        psnr=main(k,method,imtest,use_table,from,to);
        if method==1
            taw=2;
            beta=1.2;
            pt=8;
            M=3;
        elseif method==2
            taw=0.4;
            beta=0.1;
            pt=3;
            M=5;
        else
            taw=2;
            beta=0.5;
            pt=5;
            M=3;
        end
        
        set(hBetaE,'String',num2str(beta)); %return evaluated values and show
        set(hPSNRE,'String',num2str(psnr));
        set(hME,'String',num2str(M));
        set(hPtE,'String',num2str(pt));
        set(htawE,'String',num2str(taw));
    end

%callback for button 'Without EALD'
    function Start_callback2(hObject,eventdata)        
        %reading parameters
        imname=get(himlist,'Value'); %get selected value number from listbox
        Meth=get(himlist1,'Value');
        from=str2double(get(hfromE,'String'));
        to=str2double(get(htoE,'String'));
        k=str2double(get(hkE,'String'));
        switch Meth
            case 1
                method=1;
            case 2
                method=2;
            case 3
                method=3;
        end
        
        %read selected image
        switch imname
            case 1
                imtest=imread('fishingboat.tif');
            case 2
                imtest=imread('lena_gray_512.tif');
            case 3
                imtest=imread('baboon_512.tif');
        end
        
        % send required value to 'main' function for evaluate 'PSNR'
        psnr1=interpol_nf(k,imtest,method,from,to);
        if method==1
            taw=0;
            beta=0;
            pt=0;
            M=0;
        elseif method==2
            taw=0;
            beta=0;
            pt=0;
            M=0;
        else
            taw=0;
            beta=0;
            pt=0;
            M=0;
        end
        
        set(hBetaE,'String',num2str(beta)); %return evaluated values and show
        set(hPSNRE,'String',num2str(psnr1));
        set(hME,'String',num2str(M));
        set(hPtE,'String',num2str(pt));
        set(htawE,'String',num2str(taw));
    end

% Callback for button 'Show Image'
% role: show selected image
    function Showimage_callback(hObject,eventdata)
        imname=get(himlist,'Value');
        switch imname
            case 1
                imtest=imread('fishingboat.tif');
            case 2
                imtest=imread('lena_gray_512.tif');
            case 3
                imtest=imread('baboon_512.tif');
        end
        figure('Position',[800 550 256 256]);
        imshow(imtest);
       
    end


%**************************************************************************
%***************************** END Function********************************
end
%**************************************************************************
%**************************************************************************