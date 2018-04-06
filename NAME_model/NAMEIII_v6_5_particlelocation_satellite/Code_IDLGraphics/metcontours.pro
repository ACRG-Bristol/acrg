pro metcontours,field,contours1,contourcolors,fieldunits

case field of

;'Topog':begin
;contours1=[-100,10,50,100,200,400,800,1600,3200,6400]
;fieldunits='m'
;end

'Topog':begin
contours1=[-100,10,50,100,200,300,400,500,600,$
           800,1000,1200,1400,1600,1800,2000,$
           2200,2400,2600,3000,3500,4000,4500,5000,6000,7000]
fieldunits='m'
end


'Zi':begin
contours1=[0,10,50,100,200,400,800,1600,3200]
fieldunits='m'
end

'Tot ppt':begin
contours1=[0.001,0.03,0.1,0.2,0.5,1.0,2.0,4.0,8.0,16.0,32.0]
fieldunits='mm/hr'
end

'Con ppt':begin
contours1=[0.001,0.03,0.1,0.2,0.5,1.0,2.0,4.0,8.0,16.0,32.0]
fieldunits='mm/hr'
end

'Dyn ppt':begin
contours1=[0.001,0.03,0.1,0.2,0.5,1.0,2.0,4.0,8.0,16.0,32.0]
fieldunits='mm/hr'
end

'Precipitation rate (mm/hr)':begin
contours1=[0.001,0.03,0.1,0.2,0.5,1.0,2.0,4.0,8.0,16.0,32.0]
fieldunits='mm/hr'
end

'Pmsl':begin
contours1=952+(findgen(25)*4.0)
fieldunits='hpa'
end

'Heat flux':begin
contours1=[-100,-50,0,50,100,200,400,800]
fieldunits='w/m2'
end

'U*':begin
contours1=[0.0,0.05,0.1,0.2,0.3,0.5,1.0,2.0]
fieldunits='m/s'
end

'Z0':begin
contours1=[0.0,0.0001,0.001,0.01,0.1,1,10]
fieldunits='m'
end

'Ts':begin
contours1=[-100,-50,-20,-10,0,10,20,30,40,50]
fieldunits='C'
end

'Tot cld':begin
contours1=[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.8,1.0,1.1]
fieldunits='Fraction'
end

'W':begin
contours1=[-1.0,-0.3,-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3,1]
fieldunits='m/s'
end

'P':begin
contours1=[600.0,650.0,700.0,750.0,800.0,850.0,900.0,950.0,1000.0,1050.0]
fieldunits='hpa'
end

'P*':begin
contours1=[600.0,700.0,800.0,850.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0]
fieldunits='hpa'
end

'Temp':begin
contours1=[-150,-100,-50,-20,-10,0,10,20,30,40,50]
fieldunits='C'
end

'Cldwat':begin
contours1=[1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3]
fieldunits='g/g'
end

'Q':begin
contours1=[1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1]
fieldunits='g/g'
end


'Land Use':begin
contours1=[1,2,3,4,5,6,7,8,9,10,11]
fieldunits=''
end

    'PARTCLAY':begin
        contours1=[-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1]
        fieldunits='Fraction'
    end

    'PARTVEG':begin
        contours1=[-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1]
        fieldunits='Fraction'
    end


    'PARTMREL1':begin
        contours1=[-1,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
        fieldunits='Fraction'
    end

    'PARTMREL2':begin
        contours1=[-1,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
        fieldunits='Fraction'
    end



    'PARTMREL3':begin
        contours1=[-1,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
        fieldunits='Fraction'
    end


    'PARTMREL4':begin
        contours1=[-1,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
        fieldunits='Fraction'
    end



    'PARTMREL5':begin
        contours1=[-1,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
        fieldunits='Fraction'
    end



    'PARTMREL6':begin
        contours1=[-1,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
        fieldunits='Fraction'
    end


   'PARTMOISTURE':begin
        contours1=[-1,0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0]
        fieldunits='Fraction'
    end

else:begin
return
end

endcase

case field of

'Land Use':begin
contourcolors=[15,17,3,7,12,2,8,21,1,32,99]
end

else:begin
ncol=n_elements(contours1)
dcol=(!d.n_colors-35)/(ncol+1)
contourcolors=!d.n_colors-(indgen(ncol)*dcol+dcol)
end

endcase

end
