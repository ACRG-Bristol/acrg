pro plot_obs_vs_mod,obs_datadir,all_mod_datadir,start_date,end_date,year,$
  field,level,grid,all_species,all_obs_loc,all_mod_loc,plot_type,$
  plot_length,output_dir

;-------------------------------------------------------------------------------------------
; procedure initialised by plot_obs_vs_mod.com in order to start the comparison of 
; modelled data with that obtained from the netcen website. 

; calls the routines:	convert_mod_data.pro
;			extract_obs_data.pro


;Version number		Date		Owner			Comment
;0.1			Dec 2003	Mark Harrison

;-------------------------------------------------------------------------------------------
print,'starting plot_obs_vs_mod'
n_mod_dataset=n_elements(all_mod_datadir)
print,'n_mod_dataset',n_mod_dataset
all_stats_filename=STRARR(n_mod_dataset)

For s=0,n_mod_dataset-1 do begin

  st=string(s)
 
  
  ;specifying the name and location of the stats file
   stats_file=output_dir+'/stats'+st+'.txt'
   stats_filename=strcompress(stats_file,/remove_all)

   all_stats_filename(s)=stats_filename

   get_lun,fout
   openw,fout,all_stats_filename(s)

  ;writing the header for the stats file
   printf,fout,'location','species','fac2','fb','r','rmse','nmse','stderr',$
   'meanerr',Format="(A30,A20,7A10)"

  ;closing the stats file
   close,fout
   free_lun,fout 
      
   print,'stats_filename',stats_filename
   
EndFor 

;calculates the number of locations specified 
;for the mod and obs in .com routine
 n_mod_locs=n_elements(all_mod_loc)
 n_obs_locs=n_elements(all_obs_loc)

;checks that the number of locations inputted for model
;is the same as that specified for the observations
 if(n_obs_locs eq  n_mod_locs) then begin
         
  no_of_days=end_date.julian-start_date.julian
  print,'no_of_days',no_of_days
  
  no_of_hrs=int(24*no_of_days)+1
  print,'no_of_hrs',no_of_hrs

;Do these routines for each location  

  For j=0,n_obs_locs-1 do begin
      mod_loc=all_mod_loc(j)
      obs_loc=all_obs_loc(j)  
  
      ;Do these routines for each species

      nspecies=n_elements(all_species)
   
      For i=0,nspecies-1 do begin
          species=all_species(i)

     ;Initialising the matrices with null values
 
          all_mod_data=replicate(-9.9,n_mod_dataset,no_of_hrs)
          all_mod_date=replicate(end_date,n_mod_dataset,no_of_hrs)
          all_mod_units=replicate('units',n_mod_dataset,no_of_hrs)

          ;Do this routine to get mod data for each mod_datadir specified in .com
      
          ymax=1.0
	  
          For s=0,n_mod_dataset-1 do begin
             mod_dataset=all_mod_datadir(s)

	     ;calling the convert_mod_data routine to obtain and convert 
	     ;the model output to the correct units.
	 
	      convert_mod_data,mod_dataset,field,level,grid,species,mod_loc,$
	        mod_date,mod_data,mod_units
	

	
	      ;need to make sure we take the relevant time period 
	
	      indmod=where(mod_date.julian ge start_date.julian and mod_date.julian le end_date.julian,nindmod)
	
	      If nindmod gt 0  then begin
	
	        mod_data=mod_data(indmod)
                ind=where(mod_data lt 0.0,nind)
		if nind gt 0 then mod_data(ind)=-9.9

info,all_mod_data(s,*)
info,mod_data
rows=n_elements(mod_data)
info,rows
		
	        all_mod_data(s,0:rows-1)=mod_data
	      
	        all_mod_date(s,0:rows-1)=mod_date(indmod)
	      	
	        all_mod_units(s,0:rows-1)=mod_units(indmod)
             
	        ymax=max([mod_data,ymax])
		print,'Ymax:',ymax

              EndIF
	     
	  EndFor     
           
	  ;calling the extract_obs_data routine to read the observed data
	 
	   extract_obs_data,obs_datadir,year,species,obs_loc,obs_date,obs_data,$
	     obs_units,end_date
           ind=where(obs_data lt 0.0,nind)
           if nind gt 0 then obs_data(ind)=-9.9

     info,obs_units
     obs_units='wrong'
     
           

          ;Now add the general plotting information
           
	   tek_color
	   
           loc_name=STRSUBST(obs_loc,' ','_')

	   fileplot=output_dir+'/'+loc_name+'_'+species+'.ps'
	   print,'fileplot',fileplot
		
	   set_plot,'ps'
  	   device,/portrait
	   device,/color
	   device,filename=fileplot
	   device, Ysize=23
	   device,yoffset=3
	   !p.font=-1
	   !p.multi=[0,1,2]
	   ;!y.margin=[4,2]
	   !P.Charsize=0.9
	   ;!x.margin=[10,4]
	   !X.style=1
	   !Y.Style=1

           indobs=where(obs_date.julian ge start_date.julian and obs_date.julian le end_date.julian,nindobs)
	 
           if (nindobs gt 0) then ymax=max([ymax,obs_data(indobs)])
	   
           reus=no_of_days / plot_length

           print,'reus',reus

           no_of_plots=NINT(reus)
	
           leftover=no_of_plots - reus

           If (reus - no_of_plots GE 0.5) THEN BEGIN
             no_of_plots=no_of_plots+1
           EndIf

           print,'no_of_plots',no_of_plots
      
           ;Do the plotting for each species
                
           plot_start=start_date.julian
           plot_end=plot_start+plot_length

          ;coding to organise the labelling of the x-axis specifically
        
           xtickv=[1,8,15,22,29]
           xticks=N_ELEMENTS(xtickv)-1
    
           k=n_elements(obs_data)

	   col_obs=2	
		   
           ;Need to specify that the -9 values are not to be included 
	   ;in the plotting
	 
	
          ;specify the type of plot that is required

           If plot_type eq 'overplot' then begin   
             ymin=0
             obs_plot=obs_data
             obs_label=0.95*ymax
           EndIf 
    
    
           If plot_type eq 'mirrorred' then begin   
             ymin=-ymax
             obs_plot=-obs_data
             obs_label=-ymax     
           EndIf

               mod_units=all_mod_units(0,0)
	       
               print,'mod_units',mod_units  

           For z=0,no_of_plots-1 DO BEGIN

	     ;now add information regarding the axes
	   
	     plot,obs_date,obs_data,DT_Range=[plot_start,plot_end],$
               Title = mod_loc+' '+species+' '+level,Yrange=[ymin,ymax],$
	       YTitle = species+' '+mod_units,/nodata


	     ;adding the obs to the plot
	   
	      oplotm,obs_date,obs_plot,color=col_obs,thick=2.5,missing=-9.9

             ;Plotting the mod data for each set specified in the .com file
	     ;This takes the information from the relevant column of the
	     ;matrix and reassigns it to a 1D array
	     
             For s=0,n_mod_dataset-1 do begin
               mod_data=all_mod_data(s,*)
               mod_data=reform(mod_data,no_of_hrs)

               mod_date=all_mod_date(s,*)
               mod_date=reform(mod_date,no_of_hrs)


	      oplotm,mod_date,mod_data,color=s+4,thick=2.5,missing=-9.9
	  
             
	     
	
	      If z eq 0 THEN BEGIN
	
         ;Calling the calc_stats subroutine
                calc_stats,mod_data,obs_data,start_date,end_date,mod_date,$
                obs_date,r,rmse,fac2,fb,nmse,stderr,meanerr,all_stats_filename(s),obs_loc,species       
	
	      EndIf
	     	  
	     EndFor     

            XYOUTS, 1, ymax, "modelled data",color=4
            XYOUTS, 1, obs_label, "observed data",color=2


            If n_mod_dataset eq 2 then begin
              XYOUTS, 1, 0.9*ymax, "second modelled data",color=5
            EndIf


            If n_mod_dataset eq 3 then begin
              XYOUTS, 1, 0.85*ymax, "third modelled data",color=6
            EndIf

            If n_mod_dataset eq 4 then begin
              XYOUTS, 1, 0.8*ymax, "fourth modelled data",color=7
            EndIf
	    
           plot_start=plot_end
           plot_end=plot_end+plot_length

           If z eq no_of_plots-2 THEN BEGIN

             If (reus - no_of_plots LT 0.5) THEN BEGIN
                plot_end = end_date.julian
             EndIf
           EndIf
         
    
        Endfor
     

        device,/close_file  
	 

       
 
     Endfor

     Endfor

   Endif

end

